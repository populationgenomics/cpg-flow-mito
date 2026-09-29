from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob


def download_latest_annotations(output_path: Path, job_attrs: dict[str, str]):
    """Download MitoMap annotations. Only called when no config default exists."""

    batch_instance = hail_batch.get_batch()
    job = batch_instance.new_bash_job('Monthly annotation update', job_attrs | {'tool': 'mitoreport'})
    job.image(config.config_retrieve(['mito_images', 'mitoreport']))

    job.command(f"""
    set +e
    n=0
    until [ "$n" -ge 5 ]
    do
       java -jar mitoreport.jar mito-map-download --output {job.output} && break
       n=$((n+1))
       sleep 20
    done
    if [ ! -s {job.output} ]; then
        echo "MitoMap download failed after 5 attempts." >&2
        echo "Set mito_references.mito_map_annotations in config to use a static reference." >&2
        exit 1
    fi
    """)

    batch_instance.write_output(job.output, str(output_path))
    return job


def reannotate_microproteins(annotations_input: Path, output_path: Path, job_attrs: dict[str, str]) -> list[BashJob]:
    """
    Replace microprotein annotations in the MitoMap data with consequences on their parent genes.

    Three jobs:
        1. write every allele annotated on a microprotein in the mapping to a sites-only VCF
        2. annotate that VCF with bcftools csq, using the mitochondrial codon table
        3. swap in the parent gene consequences, writing a new version of the annotations
    """

    batch_instance = hail_batch.get_batch()
    annotations = batch_instance.read_input(str(annotations_input))
    fasta = config.config_retrieve(['mito_references', 'fasta'])
    reference = batch_instance.read_input_group(fa=fasta, fai=f'{fasta}.fai')
    gff3 = batch_instance.read_input(config.config_retrieve(['mito_references', 'gff3']))

    vcf_job = batch_instance.new_bash_job('Microprotein alleles to VCF', job_attrs)
    vcf_job.image(config.config_retrieve(['workflow', 'driver_image']))
    vcf_job.command(
        f'python -m cpg_flow_mito.scripts.reannotate_microproteins to-vcf '
        f'-i {annotations} -r {reference.fa} -o {vcf_job.output}',
    )

    # --local-csq: each record is annotated independently, rather than as one haplotype
    # -C 2: vertebrate mitochondrial genetic code
    # --unify-chr-names: chrM in the VCF & fasta, M in the Ensembl GFF3
    csq_job = batch_instance.new_bash_job('Annotate microprotein alleles', job_attrs | {'tool': 'bcftools'})
    csq_job.image(config.config_retrieve(['mito_images', 'bcftools']))
    csq_job.command(
        f'bcftools csq --force --no-version -f {reference.fa} -g {gff3} --local-csq -C 2 '
        f"--unify-chr-names 'chr,-,chr' -Ov -o {csq_job.output} {vcf_job.output}",
    )

    apply_job = batch_instance.new_bash_job('Reannotate microproteins', job_attrs)
    apply_job.image(config.config_retrieve(['workflow', 'driver_image']))
    apply_job.command(
        f'python -m cpg_flow_mito.scripts.reannotate_microproteins apply '
        f'-i {annotations} -v {csq_job.output} -g {gff3} -o {apply_job.output}',
    )

    batch_instance.write_output(apply_job.output, str(output_path))
    return [vcf_job, csq_job, apply_job]
