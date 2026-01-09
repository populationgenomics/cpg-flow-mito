
from cpg_utils import Path, hail_batch, config

def download_latest_annotations(output_path: Path, job_attrs: dict[str, str]):
    """Trigger the MitoMap download, save to GCP."""

    batch_instance = hail_batch.get_batch()
    job = batch_instance.new_bash_job('Monthly annotation update', job_attrs | {'tool': 'mitoreport'})
    job.image(config.config_retrieve(['images', 'mitoreport']))

    # noted that this succeeds locally, but may be fragile. Perhaps a retry wrap makes sense.
    job.command(f'java -jar mitoreport.jar mito-map-download --output {job.output}')
    batch_instance.write_output(job.output, output_path)
    return job
