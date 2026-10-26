"""
Generate a MitoReport index page.

Report paths are passed at planning time via JSON. Participant IDs are
resolved from Metamist at execution time so the index covers all SGs
in the dataset, not just the current run.
"""

import json
from argparse import ArgumentParser
from datetime import datetime, timezone
from pathlib import Path

import jinja2
from cpg_utils import to_path
from loguru import logger
from metamist.graphql import gql, query

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

PARTICIPANT_QUERY = gql(
    """
    query GetParticipants($project: String!, $sgIds: [String!]!) {
        project(name: $project) {
            sequencingGroups(id: {in_: $sgIds}) {
                id
                sample {
                    participant {
                        externalId
                    }
                }
            }
        }
    }
    """,
)


def get_participant_ids(dataset: str, sg_ids: list[str]) -> dict[str, str]:
    result = query(PARTICIPANT_QUERY, variables={'project': dataset, 'sgIds': sg_ids})
    mapping: dict[str, str] = {}
    for sg in result.get('project', {}).get('sequencingGroups', []):
        try:
            mapping[sg['id']] = sg['sample']['participant']['externalId']
        except (KeyError, TypeError):
            pass
    return mapping


def main(dataset: str, output: str, reports: dict[str, dict[str, str]], sg_ids: list[str]) -> None:
    participant_ids = get_participant_ids(dataset, sg_ids)

    for sg_id, report in reports.items():
        report['participant'] = participant_ids.get(sg_id, sg_id)
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
    parser.add_argument('--sg-ids', required=True, help='Path to JSON file with SG ID list.')
    args = parser.parse_args()

    with open(args.reports) as f:
        reports_data = json.load(f)

    with open(args.sg_ids) as f:
        sg_ids_data = json.load(f)

    main(args.dataset, args.output, reports_data, sg_ids_data)