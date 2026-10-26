"""
Generate a MitoReport index page from pre-computed report metadata.

Report paths and participant IDs are passed at planning time via JSON.
File timestamps are read from GCS at execution time (st_mtime) to avoid
a race with Metamist registration jobs.
"""

import json
from argparse import ArgumentParser
from datetime import datetime, timezone
from pathlib import Path

import jinja2
from cpg_utils import to_path
from loguru import logger

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'


def main(dataset: str, output: str, reports: dict[str, dict[str, str]]) -> None:
    for sg_id, report in reports.items():
        try:
            mtime = to_path(report['path']).stat().st_mtime
            report['timestamp'] = datetime.fromtimestamp(mtime, tz=timezone.utc).strftime('%Y-%m-%d')
        except (OSError, TypeError):
            logger.warning(f'Could not stat report for {sg_id}, omitting timestamp')
            report['timestamp'] = ''

    template_context = {
        'title': f'MitoReport index for {dataset}',
        'reports': reports,
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('mito_index.html.jinja')
    content = template.render(**template_context)

    print(f'Writing {template_context["title"]} to {output}')
    to_path(output).write_text('\n'.join(line for line in content.split('\n') if line.strip()))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Dataset to generate index page for.')
    parser.add_argument('--output', required=True, help='Path to write new HTML file to.')
    parser.add_argument('--reports', required=True, help='Path to JSON file with report data.')
    args = parser.parse_args()

    with open(args.reports) as f:
        reports_data = json.load(f)

    main(args.dataset, args.output, reports_data)
