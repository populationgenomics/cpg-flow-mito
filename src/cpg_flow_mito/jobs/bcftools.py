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
    missing_to_ref: bool = False,
    vcfs_localised: bool = False,
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
        missing_to_ref (bool): if specified, replace missing calls with unphased WT 0/0
        vcfs_localised (bool): if false, read into batch. If true, assume already in the batch
        job_attrs (dict[str, str]): attributes to pass to the job
    """

    batch_instance = hail_batch.get_batch()

    if vcfs_localised:
        batch_vcfs = input_list
    else:
        batch_vcfs = [
            batch_instance.read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
            for each_vcf in input_list
        ]

    merge_job = batch_instance.new_job('Merge VCFs', attributes=job_attrs or {} | {'tool': 'bcftools'})
    merge_job.image(config.config_retrieve(['images', 'bcftools']))

    # guessing at resource requirements
    merge_job.cpu(cpu)
    merge_job.memory(memory)
    merge_job.storage(storage)
    merge_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # -m: merge strategy
    # -0: missing-calls-to-ref (not used by default)
    merge_job.command(
        f'bcftools merge {" ".join(batch_vcfs)} -Oz -o '
        f'{merge_job.output["vcf.bgz"]} --threads {cpu} -m all {" -0" if missing_to_ref else ""} -W=tbi',
    )

    # write the result out
    batch_instance.write_output(merge_job.output, output_file.removesuffix('.vcf.bgz'))

    return merge_job
