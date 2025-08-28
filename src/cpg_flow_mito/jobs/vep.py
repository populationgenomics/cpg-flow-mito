"""
Creates a Hail Batch job to run the command line VEP tool.
"""

from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch, to_path

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def vep_one(
    vcf: Path,
    out_path: str,
    job_attrs: dict,
) -> 'BashJob':
    """Run a single VEP job."""

    local_vcf = hail_batch.get_batch().read_input(vcf)

    # check that the cache and image for this version exist
    vep_mount_path = to_path(config.config_retrieve(['references', 'vep_mount']))

    job = hail_batch.get_batch().new_bash_job('AnnotateFragmentedVcfWithVep', job_attrs | {'tool': 'VEP'})
    job.image(config.config_retrieve(['images', 'vep']))

    job.memory('16Gi').storage('15Gi').cpu(1)

    # gcsfuse works only with the root bucket, without prefix:
    data_mount = to_path(f'/{vep_mount_path.drive}')
    job.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

    job.command(
        f"""
    set -x

    vep \\
        --format vcf \\
        --vcf --compress_output bgzip \\
        -o {job.output} \\
        -i {local_vcf} \\
        --everything \\
        --mane_select \\
        --allele_number \\
        --minimal \\
        --species homo_sapiens \\
        --cache --offline --assembly GRCh38 \\
        --dir_cache {vep_dir}/vep/ \\
        --fasta /cpg-common-main/references/vep/110/mount/vep/homo_sapiens/110/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    """,
    )

    hail_batch.get_batch().write_output(job.output, out_path)

    return job
