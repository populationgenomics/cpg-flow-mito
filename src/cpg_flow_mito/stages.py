"""
Stage to call SNVs in the mitochondrial genome of a single sequencing_group.

Reimplemented version of;
https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
"""

from cpg_flow import stage, targets, workflow
from cpg_utils import Path, config, hail_batch

from cpg_flow_mito.jobs import bcftools, mito, picard, vep


@stage.stage(
    analysis_type='mito-cram',
    analysis_keys=['non_shifted_cram'],
)
class RealignMito(stage.SequencingGroupStage):
    """
    Re-align mitochondrial genome of a single sequencing_group.

    This is a re-implementation of the broad pipeline (as of 03/22) used for
    gnomAD v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
    A default config file here:
    https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json

    Mitochondrial variant calling can be subject to subtle artifacts resulting
    from mapping errors and other complexities such as NUMTs. In an attempt to minimize
    these issues, this pipeline attempts to re-implement as much of the logic, tools and
    configuration as possible.

    The main phases of analysis include:
        - Extraction of reads mapping to chrM from the main sequencing_group cram.
        - Realignment of chrM reads to chrM reference using bwa.
        - Realignment of chrM reads to a "shifted" chrM reference using bwa. The "shifted"
            reference is the same sequence but with a different linearisation point. This
            is used to overcome mapping artifacts at the boundary (aka the "control
            region" of the mt genome).
        - Generation of overall and  base level alignment statistics.

    Mitochondrial reference indexes and lift over chains were sourced from Broad
    (gs://gcp-public-data--broad-references/hg38/v0/chrM/).

    Requires:
        sequencing_group cram from Align stage.

    Outputs:
        non_shifted_cram: Sorted cram with duplicates marked aligned to the standard
            hg38 chrM.
        shifted_cram: Sorted cram with duplicates marked aligned to a copy of chrM
            with the linearisation location "shifted"
        base_level_coverage_metrics: per base coverage needed to differentiate between
            0/0 and ./. genotypes until mutect can generate a gVCF.
        coverage_metrics: CollectWgsMetrics output.
        theoretical_sensitivity_metrics: CollectWgsMetrics output.
    """

    def expected_outputs(self, sequencing_group: stage.SequencingGroup) -> dict[str, Path]:
        main = sequencing_group.dataset.prefix()
        analysis = sequencing_group.dataset.analysis_prefix()
        return {
            'non_shifted_cram': main / 'mito' / f'{sequencing_group.id}.mito.cram',
            'shifted_cram': main / 'mito' / f'{sequencing_group.id}.shifted_mito.cram',
            'base_level_coverage_metrics': main / 'mito' / f'{sequencing_group.id}.base_level_coverage.tsv',
            'coverage_metrics': analysis / 'mito' / f'{sequencing_group.id}.coverage_metrics.txt',
            'coverage_mean': analysis / 'mito' / f'{sequencing_group.id}.coverage_mean.txt',
            'coverage_median': analysis / 'mito' / f'{sequencing_group.id}.coverage_median.txt',
            'theoretical_sensitivity_metrics': analysis / 'mito' / f'{sequencing_group.id}.theoretical_sensitivity.txt',
        }

    def queue_jobs(
        self,
        sequencing_group: stage.SequencingGroup,
        inputs: stage.StageInput,  # noqa: ARG002
    ) -> stage.StageOutput:
        outputs = self.expected_outputs(sequencing_group)
        job_attrs = self.get_job_attrs(sequencing_group)
        jobs = []

        # Extract reads mapped to chrM
        subset_bam_j = mito.subset_cram_to_chrm(
            cram_path=sequencing_group.cram,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(subset_bam_j)

        # Align extracted reads to chrM using bwa
        realign_j = mito.mito_realign(
            sequencing_group_id=sequencing_group.id,
            input_bam=subset_bam_j.output_bam,
            job_attrs=job_attrs,
        )
        jobs.append(realign_j)

        # Mark duplicates (non-shifted)
        realign_mkdup_j = picard.markdup(
            sorted_bam=realign_j.output_cram,
            output_path=outputs['non_shifted_cram'],
            job_attrs=job_attrs,
        )
        jobs.append(realign_mkdup_j)

        # Align extracted reads to "shifted" chrM using bwa
        shifted_realign_j = mito.mito_realign(
            sequencing_group_id=sequencing_group.id,
            input_bam=subset_bam_j.output_bam,
            job_attrs=job_attrs,
            shifted=True,
        )
        jobs.append(shifted_realign_j)

        # Mark duplicates (shifted)
        shifted_mkdup_j = picard.markdup(
            sorted_bam=shifted_realign_j.output_cram,
            output_path=outputs['shifted_cram'],
            job_attrs=job_attrs,
            shifted=True,
        )
        jobs.append(shifted_mkdup_j)

        # Collect coverage metrics (only on non-shifted)
        coverage_metrics_job = mito.collect_coverage_metrics(
            cram=realign_mkdup_j.output_cram,
            metrics=outputs['coverage_metrics'],
            theoretical_sensitivity=outputs['theoretical_sensitivity_metrics'],
            job_attrs=job_attrs,
        )
        jobs.append(coverage_metrics_job)

        # Extract mean and median coverage (only on non-shifted)
        extract_coverage_mean_j = mito.extract_coverage_mean(
            metrics=coverage_metrics_job.metrics,
            mean_path=outputs['coverage_mean'],
            median_path=outputs['coverage_median'],
            job_attrs=job_attrs,
        )
        jobs.append(extract_coverage_mean_j)

        # CoverageAtEveryBase
        # Provides coverage at each base so low coverage sites can be
        # considered ./. rather than 0/0.
        non_control_region_coverage_j = mito.coverage_at_every_base(
            cram=realign_mkdup_j.output_cram,
            job_attrs=job_attrs,
        )
        jobs.append(non_control_region_coverage_j)

        shifted_control_region_coverage_j = mito.coverage_at_every_base(
            cram=shifted_mkdup_j.output_cram,
            job_attrs=job_attrs,
            shifted=True,
        )
        jobs.append(shifted_control_region_coverage_j)

        # Merge coverage stats
        merge_coverage_j = mito.merge_coverage(
            non_cr_coverage=non_control_region_coverage_j.per_base_coverage,
            shifted_cr_coverage=shifted_control_region_coverage_j.per_base_coverage,
            merged_coverage=outputs['base_level_coverage_metrics'],
            job_attrs=job_attrs,
        )
        jobs.append(merge_coverage_j)

        return self.make_outputs(sequencing_group, data=outputs, jobs=jobs)


@stage.stage(required_stages=[RealignMito])
class GenotypeMito(stage.SequencingGroupStage):
    """
    Call SNVs in the mitochondrial genome of a single sequencing_group.
    This is a re-implementation of the broad pipeline as of 03/22 that was used for
    gnomAD v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
    A default config file here:
    https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json

    The main phases of analysis include:
        - Calling of variants from non-shifted cram using mutect2
        - Calling of variants in control region in shifted cram using mutect2
        - Merging of the two call sets into a single vcf in normal chrM coordinate space.
        - Variant filtering part 1: Exclude black list of known problem sites and high
            alt allele counts.
        - Estimate contamination with haplocheckCLI (using known mito haplotypes)
        - Variant filtering part 2: exclude variants with VAF below contamination
            estimate.
        - Export final vcf

     Mitochondrial reference indexes and filtering black lists were copied from Broad
    (gs://gcp-public-data--broad-references/hg38/v0/chrM/).

    Requires:
        sequencing_group cram from Align stage.

    Outputs:
        out_vcf: the final filtered vcf for downstream use
        haplocheck_metrics: Metrics generated from the haplocheckCLI tool including an
            estimate of contamination and the predicted mitochondrial haplotype found.

    Configuration options:
        The following are surfaced as configurable parameters in the Broad WDL. Other
        parameters hardcoded in the WDL are also hardcoded in this pipeline.
        references.vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
        references.f_score_beta: "F-Score beta balances the filtering strategy between
            recall and precision. The relative weight of recall to precision."
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        main = sequencing_group.dataset.prefix()
        analysis = sequencing_group.dataset.analysis_prefix()
        return {
            'out_vcf': main / 'mito' / f'{sequencing_group.id}.mito.vcf.bgz',
            'haplocheck_metrics': analysis / 'mito' / f'{sequencing_group.id}.haplocheck.txt',
        }

    def queue_jobs(
        self,
        sequencing_group: targets.SequencingGroup,
        inputs: stage.StageInput,
    ) -> stage.StageOutput:
        outputs = self.expected_outputs(sequencing_group)
        jobs = []

        # Get input resources
        non_shifted_cram = inputs.as_str(sequencing_group, RealignMito, 'non_shifted_cram')
        non_shifted_cram_local = hail_batch.get_batch().read_input_group(
            cram=non_shifted_cram,
            crai=non_shifted_cram + '.crai',
        )

        shifted_cram = inputs.as_str(sequencing_group, RealignMito, 'shifted_cram')
        shifted_cram_local = hail_batch.get_batch().read_input_group(
            cram=shifted_cram,
            crai=shifted_cram + '.crai',
        )

        if config.config_retrieve(['workflow', 'use_verifybamid']):
            verify_bamid_path = (
                sequencing_group.dataset.prefix() / 'qc' / 'verify_bamid' / f'{sequencing_group.id}.verify-bamid.selfSM'
            )
            if not verify_bamid_path.exists():
                raise FileNotFoundError(f'Configured to use verifyBamID but file not found: {verify_bamid_path}')
            verifybamid_output = hail_batch.get_batch().read_input(verify_bamid_path)
        else:
            verifybamid_output = None

        # Call variants on WT genome
        call_j = mito.mito_mutect2(
            cram=non_shifted_cram_local,
            region='chrM:576-16024',  # Exclude the control region.
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(call_j)

        # Call variants in ONLY the control region using the shifted reference
        shifted_call_j = mito.mito_mutect2(
            cram=shifted_cram_local,
            region='chrM:8025-9144',  # Only call inside the control region.
            job_attrs=self.get_job_attrs(sequencing_group),
            shifted=True,
        )
        jobs.append(shifted_call_j)

        # Merge the wt and shifted VCFs
        merge_j = mito.liftover_and_combine_vcfs(
            vcf=call_j.output_vcf,
            shifted_vcf=shifted_call_j.output_vcf,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_j)

        # Merge the mutect stats output files (needed for filtering)
        merge_stats_j = mito.merge_mutect_stats(
            first_stats_file=call_j.output_vcf['vcf.gz.stats'],
            second_stats_file=shifted_call_j.output_vcf['vcf.gz.stats'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_stats_j)

        # Initial round of filtering to exclude blacklist and high alt alleles
        initial_filter_j = mito.filter_variants(
            vcf=merge_j.output_vcf,
            merged_mutect_stats=merge_stats_j.combined_stats,
            # alt_allele and vaf config hardcoded in this round of filtering as per
            # https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
            min_allele_fraction=0,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(initial_filter_j)

        # SplitMultiAllelics AND remove non-passing sites
        # Output is only used for input to haplocheck
        split_multiallelics_j = mito.split_multi_allelics(
            vcf=initial_filter_j.output_vcf,
            remove_non_pass_sites=True,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(split_multiallelics_j)

        # Estimate level of contamination from mito reads
        get_contamination_j = mito.get_contamination(
            vcf=split_multiallelics_j.output_vcf,
            haplocheck_output=outputs['haplocheck_metrics'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(get_contamination_j)

        # Parse contamination estimate reports
        parse_contamination_j, contamination_level = mito.parse_contamination_results(
            haplocheck_output=get_contamination_j.haplocheck_output,
            verifybamid_output=verifybamid_output,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(parse_contamination_j)

        # Filter round 2 - remove variants with VAF below estimated contamination
        second_filter_j = mito.filter_variants(
            vcf=initial_filter_j.output_vcf,
            merged_mutect_stats=merge_stats_j.combined_stats,
            contamination_estimate=contamination_level.as_str(),
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(second_filter_j)

        # Generate final output vcf
        split_multiallelics_j = mito.split_multi_allelics(
            vcf=second_filter_j.output_vcf,
            remove_non_pass_sites=False,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(split_multiallelics_j)

        # Write the final vcf to the bucket
        hail_batch.get_batch().write_output(
            split_multiallelics_j.output_vcf,
            str(outputs['out_vcf']).replace('.vcf.bgz', ''),
        )
        return self.make_outputs(sequencing_group, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[RealignMito, GenotypeMito],
    analysis_type='web',
    analysis_keys=['mitoreport'],
)
class MitoReport(stage.SequencingGroupStage):
    """
    Run the standalone MitoReport program on each SG

    This is not part of the Broad Mito pipeline, but generates an alternative non-seqr
    based interpretable html report of mito variants.

    Requires vep annotated individual vcf.

    https://github.com/bioinfomethods/mitoreport
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        main = sequencing_group.dataset.prefix()
        web = sequencing_group.dataset.web_prefix()
        return {
            'vep_vcf': main / 'mito' / f'{sequencing_group.id}.mito.vep.vcf.gz',
            'mitoreport': web / 'mito' / f'mitoreport-{sequencing_group.id}' / 'index.html',
        }

    def queue_jobs(
        self,
        sequencing_group: targets.SequencingGroup,
        inputs: stage.StageInput,
    ) -> stage.StageOutput:
        outputs = self.expected_outputs(sequencing_group)

        jobs = []

        vep_j = vep.vep_one(
            vcf=inputs.as_path(sequencing_group, GenotypeMito, 'out_vcf'),
            out_path=str(outputs['vep_vcf']),
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        if vep_j:
            jobs.append(vep_j)

        mitoreport_j = mito.mitoreport(
            sequencing_group=sequencing_group,
            vcf_path=outputs['vep_vcf'],
            cram_path=inputs.as_path(sequencing_group, RealignMito, 'non_shifted_cram'),
            output_path=outputs['mitoreport'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        if mitoreport_j:
            mitoreport_j.depends_on(*jobs)
            jobs.append(mitoreport_j)

        return self.make_outputs(sequencing_group, data=outputs, jobs=jobs)


@stage.stage(required_stages=[GenotypeMito], analysis_keys=['joint_vcf'], analysis_type='vcf')
class GenerateMitoJointCall(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        return {
            'joint_vcf': dataset.prefix()
            / workflow.get_workflow().name
            / dataset.get_alignment_inputs_hash()
            / self.name
            / 'merged_mito.vcf.bgz'
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        input_dict = inputs.as_dict_by_target(GenotypeMito)

        bcftools_job = bcftools.naive_merge_vcfs(
            input_list=[input_dict[sg.id]['out_vcf'] for sg in dataset.get_sequencing_groups()],
            job_attrs=self.get_job_attrs(dataset),
            output_file=str(outputs['joint_vcf']),
        )

        return self.make_outputs(dataset, data=outputs, jobs=bcftools_job)


@stage.stage(required_stages=[GenerateMitoJointCall], analysis_keys=['annotated_vcf'], analysis_type='vcf')
class AnnotateMitoJointCall(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        return {
            'annotated_vcf': dataset.prefix()
            / workflow.get_workflow().name
            / dataset.get_alignment_inputs_hash()
            / self.name
            / 'annotated_mito.vcf.bgz'
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        vep_j = vep.vep_one(
            vcf=inputs.as_path(dataset, GenerateMitoJointCall, 'joint_vcf'),
            out_path=str(outputs['annotated_vcf']),
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=vep_j)
