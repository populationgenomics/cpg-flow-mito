import json

from cpg_utils import Path, config, hail_batch, to_path
from hailtop.batch.job import BashJob


def create_index_job(
    dataset: str,
    output: Path,
    reports: dict[str, dict[str, str]],
    tmp_prefix: Path,
) -> BashJob:
    """Create a Hail Batch job to generate the MitoReport index page."""
    batch_instance = hail_batch.get_batch()

    reports_gcs = tmp_prefix / f'{dataset}_reports.json'
    to_path(reports_gcs).write_text(json.dumps(reports))

    job = batch_instance.new_bash_job(f'Generate index for {dataset}')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    reports_local = batch_instance.read_input(str(reports_gcs))
    job.command(
        f'python -m cpg_flow_mito.scripts.build_index'
        f' --dataset {dataset}'
        f' --output {job.output}'
        f' --reports {reports_local}',
    )
    batch_instance.write_output(job.output, output)
    return job
