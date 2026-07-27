"""
This script is set to run at the end of the Mito workflow.

It queries Metamist for all analysis entries representing MitoReport results, reduces those down to 1/sequencing group

This will run scoped by dataset.
"""

from argparse import ArgumentParser
from pathlib import Path

import jinja2
from cpg_utils import config, to_path
from loguru import logger
from metamist import graphql

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

REPORT_QUERY = graphql.gql(
    """
    query ReportQuery($project: String!, $metaFilter: JSON) {
      project(name: $project) {
        sequencingGroups {
            id
            analyses(type: {eq: "web"}, meta: $metaFilter) {
              id
              meta
              outputs
              timestampCompleted
            }
        }
      }
    }"""
)


def query_for_reports(dataset: str) -> dict[str, dict[str, str]]:
    """
    Execute a GQL query for all relevant MitoReport HTMLs. Minimise to one-per-SGID.

    Args:
        dataset: str, the name of the dataset/project

    Returns:
        dict of SGID: Report
    """

    report_lookup: dict[str, dict[str, str]] = {}

    access_level = config.config_retrieve(['workflow', 'access_level'])

    results = graphql.query(
        REPORT_QUERY,
        variables={
            'project': f'{dataset}-test' if access_level == 'test' else dataset,
            'metaFilter': {
                'stage': 'MitoReport',
            },
        },
    )

    # iterate over the SGs and grab all with reports
    for sg in results['project']['sequencingGroups']:
        if len(sg['analyses']) > 1:
            logger.warning(f'Multiple MitoReport analyses found for {sg["id"]}, last one will be used in index.')
        for report in sg['analyses']:
            meta_dict = report['meta']
            outpath = to_path(report['outputs']['path']).blob # mito/mitoreport-CPGXXXX/index.html

            web_bucket = config.config_retrieve(['storage', dataset, 'web_url'])
            url = f'{web_bucket}/{outpath}'

            # capture a minrep of all the fields we want to present in the report
            report_lookup[sg['id']] = {
                'participant': meta_dict['participant_id'],
                'timestamp': report['timestampCompleted'].split('T')[0],
                'url': url,
            }

    return report_lookup


def main(dataset: str, output: str) -> None:

    template_context = {
        'title': f'MitoReport index for {dataset}',
        'reports': query_for_reports(dataset),
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('mito_index.html.jinja')
    content = template.render(**template_context)

    # write to common web bucket - either attached to a single dataset, or communal
    print(f'Writing {template_context["title"]} to {output}')
    to_path(output).write_text('\n'.join(line for line in content.split('\n') if line.strip()))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Dataset to generate index page for.')
    parser.add_argument('--output', required=True, help='Path to write new HTML file to.')
    args = parser.parse_args()
    main(args.dataset, args.output)
