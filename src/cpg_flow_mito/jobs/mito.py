"""
Hail Batch jobs needed to call mitochondrial SNVs
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.resource import PythonResult

from cpg_utils import Path, to_path, config, hail_batch
from cpg_utils.config import get_config, image_path, reference_path
from cpg_utils.hail_batch import (
    command,
    fasta_res_group,
)
from cpg_flow import filetypes, resources, targets
from cpg_flow_mito.utils import get_mito_references, get_control_region_intervals


def subset_cram_to_chrM(
    cram_path: filetypes.CramPath,
    job_attrs: dict | None = None,
) -> Job:
    """
    Extract read pairs with at least one read mapping to chrM

    Specific config for selection of which reads to extract may influence the impact of
    NUMTs. Because of this we stick to the exact command used by Broad until we have a
    good reason not to.

    Args:
        cram_path: Input cram to extract chrM reads from.

    Outputs:
        job.output_bam: a ResourceGroup containing bam of reads mapped to chrM

    Cmd from:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl#LL188-L244C2

    """
    batch_instance = hail_batch.get_batch()

    job_attrs = (job_attrs or {}) | dict(tool='gatk_PrintReads')
    j = batch_instance.new_job('subset_cram_to_chrM', job_attrs)
    j.image(image_path('gatk'))

    resources.STANDARD.set_resources(j, ncpu=2)

    reference = fasta_res_group(batch_instance)
    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bai',
        },
    )

    # We are only accessing a tiny fraction of the genome. Mounting is the best option.
    bucket = cram_path.path.drive
    bucket_mount_path = to_path('/bucket')
    j.cloudfuse(bucket, str(bucket_mount_path), read_only=True)
    mounted_cram_path = bucket_mount_path / '/'.join(cram_path.path.parts[2:])

    cmd = f"""
        gatk PrintReads \
            -R {reference.base} \
            -L chrM \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -I {mounted_cram_path} \
            -O {j.output_bam.bam}

    """

    j.command(command(cmd))
    return j


def mito_realign(
    sequencing_group_id: str,
    input_bam: hb.ResourceGroup,
    job_attrs: dict[str, str],
    shifted: bool = False,
) -> Job:
    """
    Re-align reads to mitochondrial genome

    Uses bazam to extract fastq from input bam then pipes this directly to bwa for
    mapping.

    Args:
        sequencing_group_id: CPG sequencing_group id for inclusion in RG header
        input_bam: Bam for realignment
        mito_ref: Resource group containing the mito genome and bwa index to align to
        job_attrs: Dictionary of job attributes to add to the job.

    Outputs:
        job.output_cram: A sorted cram

    Bwa command from:
    https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignmentPipeline.wdl#L59

    """
    batch_instance = hail_batch.get_batch()
    j = batch_instance.new_job('mito_realign', job_attrs | {'tool': 'bwa'})

    mito_ref = get_mito_references(shifted=shifted)

    # todo into config
    j.image(image_path('bwa'))

    nthreads = resources.STANDARD.set_resources(j, ncpu=4).get_nthreads()

    cmd = f"""\
        bazam -Xmx16g \
            -n{min(nthreads, 6)} -bam {input_bam.bam} | \
        bwa \
            mem -K 100000000 -p -v 3 -t 2 -Y {mito_ref.base} \
            -R '@RG\\tID:{sequencing_group_id}\\tSM:{sequencing_group_id}' \
            - | \
        samtools view -bSu -T {mito_ref.base} - | \
        samtools sort -o {j.output_cram}
        """
    j.command(command(cmd))

    return j


def collect_coverage_metrics(
    cram: hb.ResourceGroup,
    metrics: Path | None = None,
    theoretical_sensitivity: Path | None = None,
    read_length_for_optimization: int = 151,
    coverage_cap: int = 100000,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run CollectWgsMetrics

    Args:
        cram: Input Cram
        read_length_for_optimization:  Read length used for optimization only. If this is
            too small CollectWgsMetrics might fail, but the results are not affected
            by this number. [Default: 151] [default: 100000 (from wdl)].
        coverage_cap: Treat positions with coverage exceeding this value as if they had
            coverage at this value.

    Outputs:
        metrics: output file
        theoretical_sensitivity: Undocumented CollectWgsMetrics output?
    """
    reference = get_mito_references()
    batch_instance = hail_batch.get_batch()
    j = batch_instance.new_job('collect_coverage_metrics', job_attrs | {'tool': 'picard_CollectWgsMetrics'})

    # todo into config
    j.image(image_path('picard'))

    resources.STANDARD.set_resources(j, ncpu=2)

    cmd = f"""
        picard \
            CollectWgsMetrics \
            INPUT={cram.cram} \
            VALIDATION_STRINGENCY=SILENT \
            REFERENCE_SEQUENCE={reference.base} \
            OUTPUT={j.metrics} \
            USE_FAST_ALGORITHM=true \
            READ_LENGTH={read_length_for_optimization} \
            COVERAGE_CAP={coverage_cap} \
            INCLUDE_BQ_HISTOGRAM=true \
            THEORETICAL_SENSITIVITY_OUTPUT={j.theoretical_sensitivity}

    """

    j.command(command(cmd))
    if metrics:
        batch_instance.write_output(j.metrics, str(metrics))
    if theoretical_sensitivity:
        batch_instance.write_output(j.theoretical_sensitivity, str(theoretical_sensitivity))

    return j


def extract_coverage_mean(
    metrics: hb.ResourceFile,
    mean_path: Path | None = None,
    median_path: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Extract mean and median coverage values for sequencing group

    Args:
        metrics: CollectWgsMetrics metrics output to process.
        mean_path: Output path for mean file.
        median_path: Output path for median file.

    Outputs:
        job.mean_coverage: mean coverage of chrM
        job.median_coverage: median coverage of chrM
    """
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='R')
    j = batch_instance.new_job('extract_coverage_mean', job_attrs)
    j.image(image_path('peer'))

    resources.STANDARD.set_resources(j, ncpu=2)

    cmd = f"""

    R --vanilla <<CODE
        df = read.table(
            "{metrics}",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\\\\t',nrows=1
        )
        write.table(
            floor(df[,"MEAN_COVERAGE"]),
            "{j.mean_coverage}",
            quote=F, col.names=F, row.names=F)
        write.table(
            df[,"MEDIAN_COVERAGE"],
            "{j.median_coverage}",
            quote=F, col.names=F, row.names=F)
        CODE
    """

    j.command(command(cmd))
    if mean_path:
        batch_instance.write_output(j.mean_coverage, str(mean_path))
    if median_path:
        batch_instance.write_output(j.median_coverage, str(median_path))

    return j


def coverage_at_every_base(
    cram: hb.ResourceGroup,
    job_attrs: dict,
    shifted: bool = False,
) -> Job:
    """
    Run picard CollectHsMetrics to calculate read coverage at each base.

    Args:
        cram: Input cram

    Outputs:
        job.per_base_coverage
        job.hs_metrics_out
    """
    intervals_list = (
        get_control_region_intervals().control_region_shifted
        if shifted
        else get_control_region_intervals().non_control_region
    )
    reference = get_mito_references(shifted)
    batch_instance = hail_batch.get_batch()
    j = batch_instance.new_job('coverage_at_every_base', job_attrs | {'tool': 'picard_CollectHsMetrics'})
    j.image(image_path('picard'))

    resources.STANDARD.set_resources(j, ncpu=2)

    cmd = f"""
    picard CollectHsMetrics \
      I={cram.cram} \
      R={reference.base} \
      PER_BASE_COVERAGE={j.per_base_coverage} \
      O={j.hs_metrics_out} \
      TI={intervals_list} \
      BI={intervals_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    """

    j.command(command(cmd))

    return j


def merge_coverage(
    non_cr_coverage: hb.ResourceFile,
    shifted_cr_coverage: hb.ResourceFile,
    merged_coverage: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Merge per base coverage files from non-control region and shifted control region.

    Args:
        non_cr_coverage: Per-base coverage from CollectHsMetrics for all of
            chrM excluding the control region.
        shifted_cr_coverage: Per-base coverage from CollectHsMetrics for only the chrM
            control region in shifted coord space.
        merged_coverage: Path to write merged coverage tsv.

    Outputs:
        job.merged_coverage: Merged coverage tsv.
    """
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='R')
    j = batch_instance.new_job('merge_coverage', job_attrs)
    j.image(image_path('peer'))

    resources.STANDARD.set_resources(j, ncpu=2)

    cmd = f"""
    R --vanilla <<CODE
      shift_back = function(x) {{
        if (x < 8570) {{
          return(x + 8000)
        }} else {{
          return (x - 8569)
        }}
      }}

      control_region_shifted = read.table("{shifted_cr_coverage}", header=T)
      shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
      control_region_shifted[,"pos"] = shifted_back

      beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
      end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

      non_control_region = read.table("{non_cr_coverage}", header=T)
      combined_table = rbind(beginning, non_control_region, end)
      write.table(combined_table, "{j.merged_coverage}", row.names=F, col.names=T, quote=F, sep="\\\t")

    CODE
    """

    j.command(command(cmd))
    batch_instance.write_output(j.merged_coverage, str(merged_coverage))

    return j


def mito_mutect2(
    cram: hb.ResourceGroup,
    region: str,
    max_reads_per_alignment_start: int = 75,
    job_attrs: dict | None = None,
    shifted: bool = False,
) -> Job:
    """
    Call SNPs and indels in mitochondrial genome using Mutect2 in "mitochondria-mode"

    Args:
        cram: Cram to call variants in.
        region: Coordinate string restricting the region to call variants within.
        max_reads_per_alignment_start: Mutect argument. [Default: 75].
        job_attrs
        shifted: bool, determines which reference group is collected

    Output:
        job.output_vcf: resource group containing vcf, index AND statistics file.

    Cmd from:
    https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L417-L484
    """
    reference = get_mito_references(shifted)
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='Mutect2')
    j = batch_instance.new_job('mito_mutect2', job_attrs)
    j.image(image_path('gatk'))

    res = resources.STANDARD.set_resources(j, ncpu=4)

    j.declare_resource_group(
        output_vcf={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
            'vcf.gz.stats': '{root}.vcf.gz.stats',
        },
    )

    cmd = f"""
        gatk --java-options "{res.java_mem_options()}" Mutect2 \
            -R {reference.base} \
            -I {cram.cram} \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O {j.output_vcf['vcf.gz']} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {max_reads_per_alignment_start} \
            --max-mnp-distance 0 \
            -L {region}
    """

    j.command(command(cmd))

    return j


def liftover_and_combine_vcfs(
    vcf: hb.ResourceGroup,
    shifted_vcf: hb.ResourceGroup,
    job_attrs: dict | None = None,
) -> Job:
    """
    Lifts over shifted vcf and combines it with the rest of the chrM calls.

    Args:
        vcf: chrM variants mapped to wt chrM (exl control region).
        shifted_vcf: chrM control region variants mapped to shifted chrM.
        reference: resource group for wt chrM reference.
        shift_back_chain: chain file provided with shifted genome.

    Outputs:
        job.output_vcf: Merged vcf.gz in resources.STANDARD hg38 coordinate space.

    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL360-L415C2
    """
    shift_back_chain = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'shift_back_chain']))
    reference = get_mito_references()
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='picard_LiftoverVcf')
    j = batch_instance.new_job('liftover_and_combine_vcfs', job_attrs)
    j.image(image_path('picard'))

    resources.STANDARD.set_resources(j, ncpu=4)

    j.declare_resource_group(lifted_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    j.declare_resource_group(output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})

    cmd = f"""
        picard LiftoverVcf \
            I={shifted_vcf['vcf.gz']} \
            O={j.lifted_vcf['vcf.gz']} \
            R={reference.base} \
            CHAIN={shift_back_chain} \
            REJECT={j.rejected_vcf}.vcf.gz

        picard MergeVcfs \
            I={vcf['vcf.gz']} \
            I={j.lifted_vcf['vcf.gz']} \
            O={j.output_vcf['vcf.gz']}
    """

    j.command(command(cmd))

    return j


def merge_mutect_stats(
    first_stats_file: hb.ResourceFile,
    second_stats_file: hb.ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """
    Merge stats files from two mutect runs

    Args:
        first_stats_file: Mutect stats ResourceFile
        second_stats_file: Mutect stats ResourceFile

    Outputs:
        combined_stats: Combined stats file

    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL573-L598C2
    """
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='gatk_MergeMutectStats')
    j = batch_instance.new_job('merge_stats', job_attrs)
    j.image(image_path('gatk'))

    resources.STANDARD.set_resources(j, ncpu=4)

    cmd = f"""
        gatk MergeMutectStats \
            --stats {first_stats_file} \
            --stats {second_stats_file} -O {j.combined_stats}
    """

    j.command(command(cmd))

    return j


def filter_variants(
    vcf: hb.ResourceGroup,
    merged_mutect_stats: hb.ResourceFile,
    contamination_estimate: hb.ResourceFile | None = None,
    min_allele_fraction: int | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Filter mutect variant calls

    Uses:
      gatk FilterMutectCalls to filter on variant properties
      gatk VariantFiltration to exclude a (broad supplied) black list of problem variants.

    Args:
        vcf: input vcf
        merged_mutect_stats: Mutect statistics file
        max_alt_allele_count: Maximum alt alleles per site (VariantFiltration).
        min_allele_fraction: Hard cutoff for minimum allele fraction. All sites with VAF
            less than this cutoff will be filtered. (VariantFiltration).
        contamination_estimate: Estimated sample contamination level (VariantFiltration).
            This is passed as a single float in a file as it is calculated at an earlier
            step in the pipeline.
        f_score_beta: F score beta, the relative weight of recall to precision,
            used if OPTIMAL_F_SCORE strategy is chosen (VariantFiltration). Default: 1.0
            is default provided by VariantFiltration.

    Outputs:
        job.output_vcf: filtered vcf file.

    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL486-L571C2
    Note:
        contamination_estimate pre-calculation has been moved out of this function.
    """
    max_alt_allele_count = config.config_retrieve(['mito_snv', 'max_alt_allele_count'])
    if min_allele_fraction is None:
        min_allele_fraction = config.config_retrieve(['mito_snv', 'vaf_filter_threshold'])
    f_score_beta = config.config_retrieve(['mito_snv', 'f_score_beta'])
    reference = get_mito_references()
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='gatk_FilterMutectCalls')
    j = batch_instance.new_job('filter_variants', job_attrs)
    j.image(image_path('gatk'))

    resources.STANDARD.set_resources(j, ncpu=4)

    blacklisted_sites = batch_instance.read_input_group(
        bed=reference_path('gnomad_mito/blacklist_sites'),
        idx=reference_path('gnomad_mito/blacklist_sites') + '.idx',
    )

    j.declare_resource_group(filtered_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    j.declare_resource_group(output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})

    if contamination_estimate:
        contamination_estimate_string = f'--contamination-estimate  $(cat {contamination_estimate})'
    else:
        contamination_estimate_string = ''

    cmd = f"""
      gatk --java-options "-Xmx2500m" FilterMutectCalls -V {vcf['vcf.gz']} \\
        -R {reference.base} \\
        -O {j.filtered_vcf['vcf.gz']} \\
        --stats {merged_mutect_stats} \\
        --max-alt-allele-count {max_alt_allele_count} \\
        --mitochondria-mode \\
        --min-allele-fraction  {min_allele_fraction} \\
        --f-score-beta  {f_score_beta} \\
        {contamination_estimate_string}
      gatk VariantFiltration -V {j.filtered_vcf['vcf.gz']} \\
        -O {j.output_vcf['vcf.gz']} \\
        --apply-allele-specific-filters \\
        --mask {blacklisted_sites.bed} \\
        --mask-name "blacklisted_site"
    """

    j.command(command(cmd, define_retry_function=True))

    return j


def split_multi_allelics(
    vcf: hb.ResourceGroup,
    remove_non_pass_sites: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Splits multi allelics and removes non pass sites
    Uses LeftAlignAndTrimVariants to split then optionally use SelectVariants to select
    only passing variants.
    Args:
        vcf: Input vcf file.
        reference: chrM reference fasta.
        remove_non_pass_sites:
    Output:
        output_vcf: Final vcf file.
    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L600
    """
    reference = get_mito_references()
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='gatk_SelectVariants')
    j = batch_instance.new_job('split_multi_allelics', job_attrs)
    j.image(image_path('gatk'))

    resources.STANDARD.set_resources(j, ncpu=4)

    j.declare_resource_group(split_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    # Downstream hail steps prefer explicit .bgz suffix
    j.declare_resource_group(output_vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    cmd = f"""
        gatk LeftAlignAndTrimVariants \
            -R {reference.base} \
            -V {vcf['vcf.gz']} \
            -O {j.split_vcf['vcf.gz']} \
            --split-multi-allelics \
            --dont-trim-alleles \
            --keep-original-ac
        """
    if remove_non_pass_sites:
        cmd += f"""
            gatk SelectVariants \
                -V {j.split_vcf['vcf.gz']} \
                -O {j.output_vcf['vcf.bgz']} \
                --exclude-filtered
        """
    else:
        cmd += f"""
            mv {j.split_vcf['vcf.gz']} {j.output_vcf['vcf.bgz']}
            mv {j.split_vcf['vcf.gz.tbi']} {j.output_vcf['vcf.bgz.tbi']}
        """

    j.command(command(cmd, define_retry_function=True))

    return j


def get_contamination(
    vcf: hb.ResourceGroup,
    haplocheck_output: Path | None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Uses new HaplocheckerCLI to estimate levels of contamination in mitochondria based
    on known mitochondrial haplotypes.

    Args:
        vcf: input vcf of passing variants with multi-allelics split.
        haplocheck_output: Path to write results file.

    Outputs:
        job.haplocheck_output: ResourceFile containing HaplocheckerCLI tsv report.

    Cmd from:
      https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L239
    """
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {} | dict(tool='haplocheckcli')
    j = batch_instance.new_job('get_contamination', job_attrs)
    j.image(image_path('haplocheckcli'))

    resources.STANDARD.set_resources(j, ncpu=2)

    cmd = f"""
        mv {vcf['vcf.bgz']} {vcf}.vcf.gz
        PARENT_DIR="$(dirname "{vcf['vcf.bgz']}")"
        java -jar /haplocheckCLI.jar "$PARENT_DIR"
        mv output {j.haplocheck_output}
        """

    j.command(command(cmd, define_retry_function=True))
    if haplocheck_output:
        batch_instance.write_output(j.haplocheck_output, str(haplocheck_output))

    return j


def parse_contamination_results(
    haplocheck_output: hb.ResourceFile,
    verifybamid_output: hb.ResourceFile | None = None,
    job_attrs: dict | None = None,
) -> tuple[Job, PythonResult]:
    """
    Post process halpocheck and (optionally) verifybamid reports to determine single value
    for estimated contamination that can be used for variant filtering.

    Inputs:
        haplocheck_report: native output from haplocheckCLI
        verifybamid_output: [optional] native output from verifyBamID

        Based on logic here:
        https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL523-L524

    Output:
        returns contamination level as a PythonResult
    """
    batch_instance = hail_batch.get_batch()
    job_attrs = job_attrs or {}
    j = batch_instance.new_python_job('parse_contamination_results', job_attrs)
    j.image(get_config()['workflow']['driver_image'])

    resources.STANDARD.set_resources(j, ncpu=4)

    def parse_contamination_worker(haplocheck_report: str, verifybamid_report: str | None) -> float:
        """
        Process haplocheckCLI and verifyBamID outputs to get contamination level as a
        single float.
        """
        cleaned_lines = []
        with open(haplocheck_report) as haplocheck:
            # Split line on tabs and strip double quotes
            for line in haplocheck:
                cleaned_lines.append([x.strip('"') for x in line.strip().split('\t')])
        # sanity check and reformat
        assert len(cleaned_lines) == 2, 'haplocheck report is unexpected format'
        assert len(cleaned_lines[0]) == 17, 'haplocheck report is unexpected format'
        report = dict(zip(cleaned_lines[0], cleaned_lines[1]))

        # Determine final contamination level
        if report['Contamination'] == 'YES':
            if float(report['MeanHetLevelMajor']) == 0:
                max_contamination = float(report['MeanHetLevelMinor'])
            else:
                max_contamination = 1.0 - float(report['MeanHetLevelMajor'])
        else:
            max_contamination = 0.0

        # If verifybamid_report is provided, chose the higher of the two
        if verifybamid_report:
            with open(verifybamid_report) as verifybamid:
                lines = [line.split('\t') for line in verifybamid.readlines()]
                assert len(lines) == 2
                report = dict(zip(lines[0], lines[1]))

            verifybamid_estimate = float(report['FREEMIX'])

            if verifybamid_estimate > max_contamination:
                max_contamination = verifybamid_estimate

        return max_contamination

    # Call parse_contamination_worker as pythonJob which returns contamination_level
    # as a hail PythonResult.
    contamination_level = j.call(parse_contamination_worker, haplocheck_output, verifybamid_output)

    return j, contamination_level


def mitoreport(
    sequencing_group: targets.SequencingGroup,
    vcf_path: Path,
    cram_path: Path,
    output_path: Path,
    job_attrs: dict,
) -> Job:
    """
    Run Mitoreport to generate html report of mito variants
    """
    mito_ref = get_mito_references()
    batch_instance = hail_batch.get_batch()
    j = batch_instance.new_job('mitoreport', job_attrs | {'tool': 'mitoreport'})
    j.image(image_path('mitoreport'))

    res = resources.STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    vcf = batch_instance.read_input_group(**{'vcf.gz': str(vcf_path)})
    cram = batch_instance.read_input_group(
        **{
            'cram': str(cram_path),
            'cram.crai': str(cram_path.with_suffix('.cram.crai')),
        },
    )

    cmd = f"""
        samtools view -T {mito_ref.base} -b -o {sequencing_group.id}.bam {cram['cram']}
        samtools index {sequencing_group.id}.bam

        java -jar mitoreport.jar mito-report \
            -sample {sequencing_group.id} \
            -mann resources/mito_map_annotations.json \
            -gnomad resources/gnomad.genomes.v3.1.sites.chrM.vcf.bgz \
            -vcf {vcf['vcf.gz']} \
            {sequencing_group.id}.bam ./resources/controls/*.bam

        gsutil -m cp -r 'mitoreport-{sequencing_group.id}/*' {output_path.parent}
        """

    j.command(
        command(
            cmd,
            setup_gcp=True,
        ),
    )

    return j
