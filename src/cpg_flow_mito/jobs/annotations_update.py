from cpg_utils import Path, config, hail_batch


def download_latest_annotations(output_path: Path, job_attrs: dict[str, str]):
    """Download MitoMap annotations, fall back to config reference if blocked."""

    batch_instance = hail_batch.get_batch()
    job = batch_instance.new_bash_job('Monthly annotation update', job_attrs | {'tool': 'mitoreport'})
    job.image(config.config_retrieve(['mito_images', 'mitoreport']))

    fallback = batch_instance.read_input(config.config_retrieve(['mito_references', 'mito_map_annotations']))

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
        echo "MitoMap download failed, using fallback annotations from config"
        cp {fallback} {job.output}
    fi
    """)

    batch_instance.write_output(job.output, str(output_path))
    return job
