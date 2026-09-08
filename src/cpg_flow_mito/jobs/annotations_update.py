from cpg_utils import Path, config, hail_batch


def download_latest_annotations(output_path: Path, job_attrs: dict[str, str]):
    """
    Download MitoMap annotations. Only called when no config default exists.

    Since 09-2026 this is using a custom re-implementation, as the default MitoMap server refuses connections, and
    MitoReport's download class can't switch to a new data host URL.
    """

    batch_instance = hail_batch.get_batch()
    job = batch_instance.new_bash_job('Monthly annotation update', job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    mitomap = config.config_retrieve(['mito_references', 'mitomap_server'])
    job.command(f'python -m cpg_flow_mito.scripts.generate_annotation_data -o {job.output} --mito-map-host {mitomap}')
    batch_instance.write_output(job.output, str(output_path))
    return job
