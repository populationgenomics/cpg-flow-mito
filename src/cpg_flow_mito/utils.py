"""
suggested location for any utility methods or constants used across multiple stages
"""

from typing import TYPE_CHECKING
from functools import cache

from cpg_utils import hail_batch, config

from datetime import datetime

if TYPE_CHECKING:
    from hailtop.batch import ResourceGroup

DATE_STRING: str = datetime.now().strftime('%y-%m')  # noqa: DTZ005


@cache
def get_mito_references(shifted: bool = False) -> ResourceGroup:
    """
    get various mito config entries, reads them into the current batch
    single method switches between shifted and non-shifted references
    Args:
        shifted: whether to get the shifted reference files

    Returns:
        dict: mito config entries
    """
    shifted_str = 'shifted_' if shifted else ''
    mito_fa = config.config_retrieve(['references', f'{shifted_str}fasta'])
    return hail_batch.get_batch().read_input_group(
        dict=config.config_retrieve(['references', f'{shifted_str}dict']),
        base=mito_fa,
        amb=mito_fa + '.amb',
        ann=mito_fa + '.ann',
        bwt=mito_fa + '.bwt',
        fai=mito_fa + '.fai',
        pac=mito_fa + '.pac',
        sa=mito_fa + '.sa',
    )


@cache
def get_control_region_intervals() -> ResourceGroup:
    """
    get mito control region intervals
    Returns:
        dict: mito control region intervals
    """

    return hail_batch.get_batch().read_input_group(
        control_region_shifted=config.config_retrieve(['references', 'shifted_control_region_interval']),
        non_control_region=config.config_retrieve(['references', 'non_control_region_interval']),
    )
