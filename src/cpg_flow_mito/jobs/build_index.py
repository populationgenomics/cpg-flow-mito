import json

from cpg_utils import Path, config, hail_batch, to_path
from hailtop.batch.job import BashJob


def create_index_job(
    dataset: str,
    output: Path,
    report_outputs: dict[str, dict[str, Path]],
    tmp_prefix: Path,
) -> BashJob:
    """Create a Hail Batch job to generate the MitoReport index page."""
    batch_instance = hail_batch.get_batch()

    web_url = config.config_retrieve(['storage', dataset, 'web_url']).rstrip('/')

    sg_ids = list(report_outputs.keys())
    reports: dict[str, dict[str, str]] = {}
    for sg_id, sg_outputs in report_outputs.items():
        reports[sg_id] = {
            'url': f'{web_url}/mito/mitoreport-{sg_id}/index.html',
            'path': str(sg_outputs['mitoreport']),
        }

    reports_gcs = tmp_prefix / f'{dataset}_reports.json'
    to_path(reports_gcs).write_text(json.dumps(reports))

    sg_ids_gcs = tmp_prefix / f'{dataset}_sg_ids.json'
    to_path(sg_ids_gcs).write_text(json.dumps(sg_ids))

    job = batch_instance.new_bash_job(f'Generate index for {dataset}')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    reports_local = batch_instance.read_input(str(reports_gcs))
    sg_ids_local = batch_instance.read_input(str(sg_ids_gcs))
    job.command(
        f'python -m cpg_flow_mito.scripts.build_index'
        f' --dataset {dataset}'
        f' --output {job.output}'
        f' --reports {reports_local}'
        f' --sg-ids {sg_ids_local}',
    )
    batch_instance.write_output(job.output, output)
    return job