from cpg_utils import Path, config, hail_batch


def download_latest_annotations(output_path: Path, job_attrs: dict[str, str]):
    """Trigger the MitoMap download, save to GCP."""

    batch_instance = hail_batch.get_batch()
    job = batch_instance.new_bash_job('Monthly annotation update', job_attrs | {'tool': 'mitoreport'})
    job.image(config.config_retrieve(['mito_images', 'mitoreport']))

    # noted that this succeeds locally, but may be fragile. Perhaps a retry wrap makes sense.
    job.command(f"""
    n=0
    until [ "$n" -ge 5 ]
    do
       java -jar mitoreport.jar mito-map-download --output {job.output} && break
       n=$((n+1))
       sleep 20
    done
    """)
    batch_instance.write_output(job.output, output_path)
    return job
