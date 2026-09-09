from cpg_utils import Path, config, hail_batch


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


def remap_microproteins(annotations_input, reference_input, output_path: Path, job_attrs: dict[str, str]):
    """Remap microprotein loci to canonical genes in the downloaded MitoMap annotations."""

    batch_instance = hail_batch.get_batch()
    job = batch_instance.new_bash_job('Remap microprotein annotations', job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    job.command(
        f'python -m cpg_flow_mito.scripts.remap_microproteins '
        f'-i {annotations_input} -r {reference_input} -o {job.output}',
    )

    batch_instance.write_output(job.output, str(output_path))
    return job
