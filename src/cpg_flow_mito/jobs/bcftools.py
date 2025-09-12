"""
general tasks using bcftools in a pipeline context
"""

from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob


def naive_merge_vcfs(
    input_list: list[str | Path],
    output_file: str,
    cpu: int = 4,
    memory: str = '16Gi',
    storage: str = '50Gi',
    job_attrs: dict[str, str] | None = None,
) -> BashJob:
    """
    A generic method to merge multiple VCFs with non-overlapping sample sets to create one multi-sample VCF.
    Sample names should be unique across all files.

    Args:
        input_list (list[str]): all VCFs to merge, can be pre-localised (see vcfs_localised[bool)
        output_file (str): path to a vcf.bgz file to write to
        cpu (int): number of cores (threads when merging)
        memory (str): RAM requirement
        storage (str): storage requirement for the task
        job_attrs (dict[str, str]): attributes to pass to the job
    """

    batch_instance = hail_batch.get_batch()

    batch_vcfs = [batch_instance.read_input(each_vcf) for each_vcf in input_list]

    merge_job = batch_instance.new_job('Merge VCFs', attributes=job_attrs or {} | {'tool': 'bcftools'})
    merge_job.image(config.config_retrieve(['images', 'bcftools']))

    # guessing at resource requirements
    merge_job.cpu(cpu)
    merge_job.memory(memory)
    merge_job.storage(storage)
    merge_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # step 1, purge the VCFs of any post-splitting previously multiallelic sites
    reduced_vcfs = []
    for index, vcf in enumerate(batch_vcfs):
        # "bcftools view -H mito_CPG262550.mito.vcf.gz | awk -F'\t' '{split($10, a, ":"); if (gsub("/", "", a[1]) < 2) print}'"
        merge_job.command(f"""
        bcftools view -h {vcf} -Oz -o ${{BATCH_TMPDIR}}/{index}.vcf.bgz
        bcftools view -H {vcf} | awk -F'\t' '{{split($10, a, ":"); if (gsub("/", "", a[1]) < 2) print}}' | bgzip >> ${{BATCH_TMPDIR}}/{index}.vcf.bgz
        bcftools index -t ${{BATCH_TMPDIR}}/{index}.vcf.bgz
        """)
        reduced_vcfs.append(f'${{BATCH_TMPDIR}}/{index}.vcf.bgz')

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # -m: merge strategy
    # -0: missing-calls-to-ref, not important for inheritance checking, but useful for AC/AN/AF accuracy
    merge_job.command(
        f'bcftools merge {" ".join(reduced_vcfs)} -Oz -o {merge_job.output["vcf.bgz"]} --threads {cpu} -0 -W=tbi',
    )

    # write the result out
    batch_instance.write_output(merge_job.output, output_file.removesuffix('.vcf.bgz'))

    return merge_job
