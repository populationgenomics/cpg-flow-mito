import json

from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob


def create_index_job(dataset: str, output: Path, reports: dict[str, dict[str, str]]) -> BashJob:
    """Create a Hail Batch job to generate the MitoReport index page."""
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_bash_job(f'Generate index for {dataset}')
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    reports_json = json.dumps(reports)
    job.command(f"""cat <<'REPORTS_JSON_EOF' > /tmp/reports.json
{reports_json}
REPORTS_JSON_EOF
python -m cpg_flow_mito.scripts.build_index --dataset {dataset} --output {job.output} --reports /tmp/reports.json""")

    batch_instance.write_output(job.output, output)
    return job
