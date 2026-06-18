from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob


def create_index_job(dataset: str, output: Path) -> BashJob:
    """Quick job to generate a MitoReport index generation job."""
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_bash_job(f'Generate index for {dataset}')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f'python -m cpg_flow_mito.scripts.build_index --dataset {dataset} --output {job.output}')
    batch_instance.write_output(job.output, output)
    return job
