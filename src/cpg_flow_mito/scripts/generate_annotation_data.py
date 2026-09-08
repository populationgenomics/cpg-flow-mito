"""
Download MitoMap annotation sources and integrate them into a single JSON file.

This is a Python port of the MitoReport ``MitoMapAnnotationsLoader`` Groovy class.
See https://github.com/bioinfomethods/mitoreport/blob/main/src/main/groovy/mitoreport/MitoMapAnnotationsLoader.groovy

It fetches:

- the Coding, Control and RNA-mutation variant tables (HTML pages embedding a DataTables JSON payload)
- the disease association TSV (``disease.cgi``)
- the MitoTIP scores TSV

and joins the disease and MitoTIP records onto each variant by its compact allele (e.g. ``A3243G``).
The result is written as a JSON list of annotation objects, one per unique variant.

The MitoMap host is configurable with ``--mito-map-host``; the page paths are fixed relative to that server.
"""

import json
import re
import shutil
import tempfile
from argparse import ArgumentParser
from dataclasses import dataclass, field
from datetime import datetime, timezone
from decimal import ROUND_HALF_EVEN, Decimal, InvalidOperation
from enum import Enum
from pathlib import Path
from time import sleep
from typing import Any

import requests
import urllib3
from loguru import logger

DEFAULT_MITO_MAP_HOST = 'https://fr.mitomap.org'
CODINGS_PAGE_PATH = '/foswiki/bin/view/MITOMAP/VariantsCoding'
CONTROLS_PAGE_PATH = '/foswiki/bin/view/MITOMAP/VariantsControl'
RNA_MUTATIONS_PAGE_PATH = '/foswiki/bin/view/MITOMAP/MutationsRNA'
DISEASES_PAGE_PATH = '/cgi-bin/disease.cgi'
MITO_TIPS_PAGE_PATH = '/downloads/mitotip_scores.txt'

USER_AGENT = 'https://www.mcri.edu.au'
DOWNLOAD_ATTEMPTS = 5
DOWNLOAD_RETRY_WAIT_SECONDS = 20
DOWNLOAD_TIMEOUT_SECONDS = 300

# Column titles (as they appear in the DataTables "columns" definition, after trimming) -> annotation attribute
TITLE_TO_PROPERTY_NAMES: dict[str, str] = {
    'Position': 'positionStr',
    'Locus': 'locusAnchor',
    'Nucleotide Change': 'alleleChange',
    'Allele': 'alleleStr',
    'Codon Number': 'codonNumber',
    'Codon Position': 'codonPosition',
    'Amino Acid Change': 'aminoAcidChange',
    "GB Freq<span class='mark'>&Dagger;</span>": 'gbFreqStr',
    "GB Freq<br><span style='white-space:nowrap;'>FL&nbsp;(CR)<span class='mark'>&ast;&Dagger;</span></span>": 'gbFreqStr',  # noqa: E501
    "GB&nbsp;Freq&nbsp;&nbsp;<br><span style='white-space:nowrap;'>FL&nbsp;(CR)<span class='mark'>&ast;&Dagger;</span></span>": 'gbFreqStr',  # noqa: E501
    'GB Seqs': 'gbSeqsAnchor',
    "GB Seqs<br><span style='white-space:nowrap;'>total&nbsp;(FL/CR)<span class='mark'>&ast;</span></span>": 'gbSeqsAnchor',  # noqa: E501
    "GB&nbsp;Seqs&nbsp;<br><span style='white-space:nowrap;'>FL&nbsp;(CR)<span class='mark'>&ast;</span></span>": 'gbSeqsAnchor',  # noqa: E501
    'Curated References': 'curatedRefsAnchor',
    'References': 'curatedRefsAnchor',
    'Disease': 'disease',
}

# Regexes used to pull the DataTables payload out of the variant HTML pages
DATA_PATTERN = re.compile(r'"data":(\[\s*?\[.*?\]\])', re.DOTALL)
COLUMNS_PATTERN = re.compile(r'"columns": (.*?}\])', re.DOTALL)

# Regexes used to derive attributes from raw cell values
HTML_ANCHOR_PATTERN = re.compile(r"""<a.*?href=['"](.*?)['"].*?>(\d+?\.?\d*?%*)</a>""")
HTML_TEXT_ANCHOR_PATTERN = re.compile(r"""<a.*?href=['"](.*?)["'].*?>(.*?)</a>""")
CONTROL_GB_FREQ_PATTERN = re.compile(r'\(.*?(\d+\.\d*?)%.*?\)')
ALLELE_CHANGE_PATTERN = re.compile(r'([ATCGU]+?)-([ATCGU|del]*)')
ALLELE_STR_PATTERN = re.compile(r'([ATCGU]+?)(\d+)([ATCGU|del]*)')

DELETION_MARKER = ':'
DELETION_STR = 'del'


class MitoTipQuartile(str, Enum):
    """MitoTIP score quartile, UNKNOWN when no MitoTIP record exists for an allele."""

    Q1 = 'Q1'
    Q2 = 'Q2'
    Q3 = 'Q3'
    Q4 = 'Q4'
    UNKNOWN = 'UNKNOWN'

    @classmethod
    def safe_value_of(cls, value: str | None) -> 'MitoTipQuartile':
        try:
            return cls(value)
        except ValueError:
            return cls.UNKNOWN


def to_decimal(value: str | None) -> Decimal | None:
    """Parse a string as a Decimal, returning None if it isn't a finite number."""
    if value is None:
        return None
    try:
        result = Decimal(value.strip())
    except InvalidOperation:
        return None
    return result if result.is_finite() else None


def to_int(value: str | None) -> int | None:
    """Parse a string as an int, returning None if it isn't one."""
    if value is None:
        return None
    try:
        return int(value.strip())
    except ValueError:
        return None


@dataclass
class MitoMapAnnotation:
    """
    A single MitoMap variant with derived attributes.

    Raw attributes mirror the columns scraped from the MitoMap variant pages; the properties derive
    the values that MitoReport consumes (compact allele, HGVS, GenBank frequency, curated refs, ...).
    Attribute names deliberately use the camelCase of the original Groovy model so the JSON output
    is compatible with MitoReport.
    """

    mitoMapHost: str  # noqa: N815
    regionType: str = 'CODING'  # noqa: N815
    positionStr: str | None = None  # noqa: N815
    locusAnchor: str | None = None  # noqa: N815
    alleleChange: str | None = None  # noqa: N815
    alleleStr: str | None = None  # noqa: N815
    codonNumber: str | None = None  # noqa: N815
    codonPosition: str | None = None  # noqa: N815
    aminoAcidChange: str | None = None  # noqa: N815
    gbFreqStr: str | None = None  # noqa: N815
    gbSeqsAnchor: str | None = None  # noqa: N815
    curatedRefsAnchor: str | None = None  # noqa: N815
    diseaseAaChange: str | None = None  # noqa: N815
    disease: str | None = None
    diseaseStatus: str | None = None  # noqa: N815
    mitoTipScore: Decimal | None = None  # noqa: N815
    mitoTipScorePercentile: Decimal | None = None  # noqa: N815
    mitoTipQuartile: MitoTipQuartile = MitoTipQuartile.UNKNOWN  # noqa: N815
    mitoTipCount: int | None = None  # noqa: N815
    mitoTipFreqPct: Decimal | None = None  # noqa: N815

    _parsed_allele: dict[str, str | None] | None = field(default=None, init=False, repr=False, compare=False)

    @property
    def identity(self) -> tuple[str | None, str | None, str | None]:
        """Key used to de-duplicate variants across the coding/control/RNA pages."""
        return (self.positionStr, self.alleleChange, self.alleleStr)

    @property
    def position(self) -> int:
        return to_int(self.positionStr) or 0

    @property
    def parsed_allele(self) -> dict[str, str | None]:
        """Ref/alt alleles parsed from 'Nucleotide Change' (e.g. 'G-A') or failing that 'Allele' (e.g. 'A3243G')."""
        if self._parsed_allele is not None:
            return self._parsed_allele

        parsed: dict[str, str | None] = {}
        if self.alleleChange and (match := ALLELE_CHANGE_PATTERN.search(self.alleleChange)):
            parsed = {'ref': match.group(1), 'alt': match.group(2) or None}
        elif self.alleleStr and (match := ALLELE_STR_PATTERN.search(self.alleleStr)):
            parsed = {'ref': match.group(1), 'alt': match.group(3) or None}

        self._parsed_allele = parsed
        return parsed

    @property
    def ref_allele(self) -> str | None:
        return self.parsed_allele.get('ref')

    @property
    def alt_allele(self) -> str | None:
        return self.parsed_allele.get('alt')

    @property
    def compact_allele(self) -> str | None:
        if not self.ref_allele or not self.position:
            return None
        return f'{self.ref_allele}{self.position}{self.alt_allele or ""}'

    @property
    def hgvs(self) -> str | None:
        if not self.position or not self.ref_allele or not self.alt_allele:
            return None
        if self.alt_allele == DELETION_STR:
            return f'm.{self.position}{self.ref_allele}del'
        return f'm.{self.position}{self.ref_allele}>{self.alt_allele}'

    @property
    def gb_freq_pct(self) -> Decimal:
        """GenBank frequency as a percentage. Control/RNA pages report 'FL%<br>(CR%)'; the CR value is used."""
        gb_freq_str = self.gbFreqStr or ''
        if self.regionType == 'CODING':
            anchor_match = HTML_ANCHOR_PATTERN.search(gb_freq_str)
            number_str = anchor_match.group(2) if anchor_match else gb_freq_str
            return to_decimal(number_str.replace('%', '')) or Decimal(0)
        if self.regionType in ('CONTROL', 'RNA_MUTATIONS'):
            control_match = CONTROL_GB_FREQ_PATTERN.search(gb_freq_str)
            if not control_match:
                return Decimal(0)
            return to_decimal(control_match.group(1).replace('%', '')) or Decimal(0)

        logger.warning(f'Unexpected regionType in gbFreqPct, regionType={self.regionType}, gbFreqStr={gb_freq_str}')
        return Decimal(0)

    @property
    def gb_freq(self) -> Decimal:
        return self.gb_freq_pct / Decimal(100)

    @property
    def curated_refs_count(self) -> int:
        anchor_match = HTML_ANCHOR_PATTERN.search(self.curatedRefsAnchor or '')
        if not anchor_match:
            return 0
        return to_int(anchor_match.group(2)) or 0

    @property
    def curated_refs_url(self) -> str:
        anchor_match = HTML_ANCHOR_PATTERN.search(self.curatedRefsAnchor or '')
        if not anchor_match:
            return self.mitoMapHost
        return f'{self.mitoMapHost}{anchor_match.group(1)}'

    @property
    def locus(self) -> str | None:
        anchor_match = HTML_TEXT_ANCHOR_PATTERN.search(self.locusAnchor or '')
        return anchor_match.group(2) if anchor_match else None

    @property
    def diseases(self) -> list[str]:
        if not self.disease:
            return []
        return [each.strip() for each in self.disease.split('+') if each.strip()]

    @property
    def disease_confirmed_pathogenic(self) -> bool:
        return (self.diseaseStatus or '').upper() == 'CFRM'

    @property
    def mito_tip_freq(self) -> Decimal | None:
        if self.mitoTipFreqPct is None:
            return None
        return self.mitoTipFreqPct / Decimal(100)

    def to_dict(self) -> dict[str, Any]:
        """Serialise raw and derived attributes, matching the Groovy JSON output."""
        return {
            'mitoMapHost': self.mitoMapHost,
            'regionType': self.regionType,
            'positionStr': self.positionStr,
            'position': self.position,
            'locusAnchor': self.locusAnchor,
            'locus': self.locus,
            'alleleChange': self.alleleChange,
            'alleleStr': self.alleleStr,
            'refAllele': self.ref_allele,
            'altAllele': self.alt_allele,
            'compactAllele': self.compact_allele,
            'allele': self.compact_allele,
            'hgvs': self.hgvs,
            'codonNumber': self.codonNumber,
            'codonPosition': self.codonPosition,
            'aminoAcidChange': self.aminoAcidChange,
            'gbFreqStr': self.gbFreqStr,
            'gbFreqPct': self.gb_freq_pct,
            'gbFreq': self.gb_freq,
            'gbSeqsAnchor': self.gbSeqsAnchor,
            'curatedRefsAnchor': self.curatedRefsAnchor,
            'curatedRefsCount': self.curated_refs_count,
            'curatedRefsUrl': self.curated_refs_url,
            'curatedRef': {'count': self.curated_refs_count, 'url': self.curated_refs_url},
            'diseaseAaChange': self.diseaseAaChange,
            'disease': self.disease,
            'diseases': self.diseases,
            'diseaseStatus': self.diseaseStatus,
            'diseaseConfirmedPathogenic': self.disease_confirmed_pathogenic,
            'mitoTipScore': self.mitoTipScore,
            'mitoTipScorePercentile': self.mitoTipScorePercentile,
            'mitoTipQuartile': self.mitoTipQuartile.value,
            'mitoTipCount': self.mitoTipCount,
            'mitoTipFreqPct': self.mitoTipFreqPct,
            'mitoTipFreq': self.mito_tip_freq,
        }


def download_page(page_url: str) -> str:
    """
    GET a page as text, retrying on failure.

    SSL verification is disabled to match the original loader, which ignored certificate issues on the MitoMap host.
    """
    logger.info(f'Downloading MitoMap page from {page_url}')
    last_error: Exception | None = None
    for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
        try:
            response = requests.get(
                page_url,
                headers={'User-Agent': USER_AGENT},
                verify=False,  # noqa: S501
                timeout=DOWNLOAD_TIMEOUT_SECONDS,
            )
            response.raise_for_status()
            return response.content.decode('utf-8', errors='replace')
        except requests.RequestException as err:
            last_error = err
            logger.warning(f'Attempt {attempt}/{DOWNLOAD_ATTEMPTS} failed for {page_url}: {err}')
            if attempt < DOWNLOAD_ATTEMPTS:
                sleep(DOWNLOAD_RETRY_WAIT_SECONDS)

    raise RuntimeError(f'Error downloading page {page_url}') from last_error


def parse_variants_html_page(html_text: str, region_type: str, mito_map_host: str) -> list[MitoMapAnnotation]:
    """
    Extract the DataTables 'columns' and 'data' JSON blobs embedded in a MitoMap variants page.

    Each data row is mapped to annotation attributes via the column titles in TITLE_TO_PROPERTY_NAMES.
    Columns with unrecognised titles are ignored.
    """
    data_match = DATA_PATTERN.search(html_text)
    columns_match = COLUMNS_PATTERN.search(html_text)
    if not data_match or not columns_match:
        raise ValueError(f'Could not locate DataTables "data"/"columns" payload in the {region_type} page')

    data = json.loads(data_match.group(1))
    # the columns block is written into a JS string, so single quotes arrive escaped (\') which is not valid JSON
    columns = json.loads(columns_match.group(1).replace("\\'", "'"))

    property_names = [TITLE_TO_PROPERTY_NAMES.get((column.get('title') or '').strip()) for column in columns]
    unmapped = [column.get('title') for column, name in zip(columns, property_names) if name is None]  # noqa: B905
    if unmapped:
        logger.debug(f'{region_type}: ignoring unmapped columns {unmapped}')

    annotations = []
    for row in data:
        attributes: dict[str, Any] = {}
        for property_name, value in zip(property_names, row):  # noqa: B905
            if property_name is None:
                continue
            attributes[property_name] = str(value).strip() if value is not None else None
        annotations.append(MitoMapAnnotation(mitoMapHost=mito_map_host, regionType=region_type, **attributes))

    logger.info(f'Parsed {len(annotations)} {region_type} variants')
    return annotations


def parse_diseases_tsv(diseases_tsv: str) -> dict[str, dict[str, str | None]]:
    """
    Parse disease.cgi output, keyed by compact allele.

    Columns: id pos ref alt aachange homoplasmy heteroplasmy disease status pubmed_ids gbcnt gbfreq
    Only aachange, disease and status are retained; the others are already present in other annotations.
    A ':' alt denotes a deletion and is normalised to 'del' to match the variant pages.
    """
    result: dict[str, dict[str, str | None]] = {}
    for line in diseases_tsv.splitlines()[1:]:
        if not line.strip():
            continue
        items = line.split('\t')
        pos = items[1] if len(items) > 1 else None
        ref = items[2] if len(items) > 2 else None
        alt = items[3] if len(items) > 3 else None
        aa_change = items[4] if len(items) > 4 else None
        disease = items[7] if len(items) > 7 else None
        mito_map_status = items[8] if len(items) > 8 else None
        if alt == DELETION_MARKER:
            alt = DELETION_STR
        compact_allele = f'{ref}{pos}{alt}'
        result[compact_allele] = {'diseaseAaChange': aa_change, 'disease': disease, 'mitoMapStatus': mito_map_status}

    logger.info(f'Parsed {len(result)} disease annotations')
    return result


def parse_mito_tips_tsv(mito_tips_tsv: str) -> dict[str, dict[str, Any]]:
    """
    Parse the MitoTIP scores TSV, keyed by compact allele, and assign each record a score percentile.

    Columns: Position rCRS Alt MitoTIP_Score Quartile Count Percentage Mitomap_Status
    Records are ranked by descending score (ties broken by allele) and the percentile is the share of
    records ranked strictly below, rounded half-even to 2 decimal places.
    """
    raw: dict[str, dict[str, Any]] = {}
    for line in mito_tips_tsv.splitlines()[1:]:
        if not line.strip():
            continue
        items = line.split('\t')
        pos = items[0] if len(items) > 0 else None
        ref = items[1] if len(items) > 1 else None
        alt = items[2] if len(items) > 2 else None
        score = (to_decimal(items[3]) if len(items) > 3 else None) or Decimal(0)
        quartile = MitoTipQuartile.safe_value_of(items[4]) if len(items) > 4 else MitoTipQuartile.UNKNOWN
        count = (to_decimal(items[5]) if len(items) > 5 else None) or Decimal(0)
        freq_pct = (to_decimal(items[6]) if len(items) > 6 else None) or Decimal(0)
        compact_allele = f'{ref}{pos}{DELETION_STR if alt == DELETION_MARKER else alt}'
        raw[compact_allele] = {
            'mitoTipScore': score,
            'mitoTipQuartile': quartile,
            'mitoTipCount': int(count),
            'mitoTipFreqPct': freq_pct,
        }

    total = len(raw)
    ranked = sorted(raw.items(), key=lambda item: (-item[1]['mitoTipScore'], item[0]))
    for index, (_allele, record) in enumerate(ranked):
        percentile = Decimal(total - index - 1) / Decimal(total) * 100
        record['mitoTipScorePercentile'] = percentile.quantize(Decimal('0.01'), rounding=ROUND_HALF_EVEN)

    logger.info(f'Parsed {total} MitoTIP annotations')
    return dict(ranked)


def build_annotations(mito_map_host: str) -> list[MitoMapAnnotation]:
    """Download all sources from the host and merge them into a de-duplicated list of annotations."""
    codings_html = download_page(f'{mito_map_host}{CODINGS_PAGE_PATH}')
    controls_html = download_page(f'{mito_map_host}{CONTROLS_PAGE_PATH}')
    rna_mutations_html = download_page(f'{mito_map_host}{RNA_MUTATIONS_PAGE_PATH}')

    codings = parse_variants_html_page(codings_html, 'CODING', mito_map_host)
    controls = parse_variants_html_page(controls_html, 'CONTROL', mito_map_host)
    rna_mutations = parse_variants_html_page(rna_mutations_html, 'RNA_MUTATIONS', mito_map_host)

    # de-duplicate on (position, nucleotide change, allele); first occurrence wins
    all_annotations: dict[tuple[str | None, str | None, str | None], MitoMapAnnotation] = {}
    for annotation in rna_mutations + codings + controls:
        all_annotations.setdefault(annotation.identity, annotation)
    logger.info(f'{len(all_annotations)} unique variants across all pages')

    diseases_annotations = parse_diseases_tsv(download_page(f'{mito_map_host}{DISEASES_PAGE_PATH}'))
    mito_tips_annotations = parse_mito_tips_tsv(download_page(f'{mito_map_host}{MITO_TIPS_PAGE_PATH}'))

    for annotation in all_annotations.values():
        disease_annotation = diseases_annotations.get(annotation.compact_allele or '', {})
        annotation.diseaseAaChange = disease_annotation.get('diseaseAaChange')
        annotation.disease = disease_annotation.get('disease')
        annotation.diseaseStatus = disease_annotation.get('mitoMapStatus')

        mito_tip_annotation = mito_tips_annotations.get(annotation.compact_allele or '', {})
        annotation.mitoTipScore = mito_tip_annotation.get('mitoTipScore')
        annotation.mitoTipScorePercentile = mito_tip_annotation.get('mitoTipScorePercentile')
        annotation.mitoTipQuartile = mito_tip_annotation.get('mitoTipQuartile', MitoTipQuartile.UNKNOWN)
        annotation.mitoTipCount = mito_tip_annotation.get('mitoTipCount')
        annotation.mitoTipFreqPct = mito_tip_annotation.get('mitoTipFreqPct')

    return list(all_annotations.values())


def json_default(value: Any) -> Any:
    """JSON encoder hook: Decimals are written as plain numbers."""
    if isinstance(value, Decimal):
        return int(value) if value == value.to_integral_value() and value.as_tuple().exponent >= 0 else float(value)  # type: ignore  # noqa: PGH003
    raise TypeError(f'Object of type {type(value).__name__} is not JSON serializable')


def write_to_file(output_path: Path, file_contents: str) -> None:
    """Write via a temp file then move into place, so a partial download never leaves a truncated output."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', prefix='mitoreport', dir=output_path.parent, delete=False) as handle:
        handle.write(file_contents)
        temp_path = Path(handle.name)
    shutil.move(str(temp_path), str(output_path))


def download_annotations(output_path: Path, mito_map_host: str) -> None:
    """Generate the annotations JSON at output_path, unless a non-empty file already exists there."""
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info(f'Skipping download, MitoMap Variants already exists at {output_path}')
        return

    annotations = build_annotations(mito_map_host)
    json_text = json.dumps([annotation.to_dict() for annotation in annotations], indent=2, default=json_default)
    write_to_file(output_path, json_text)
    logger.info(f'Wrote {len(annotations)} annotations to {output_path}')


def cli_main() -> None:
    parser = ArgumentParser(description='Download MitoMap annotations and integrate them into a single JSON file')
    parser.add_argument(
        '-o',
        '--output',
        type=Path,
        default=None,
        help='Path to save the annotations JSON to. Skips if a non-empty file already exists. '
        'Defaults to ./mito_map_annotations_<YYYYMMDD>.json',
    )
    parser.add_argument(
        '--mito-map-host',
        default=DEFAULT_MITO_MAP_HOST,
        help=f'Scheme and host serving the MitoMap pages (default: {DEFAULT_MITO_MAP_HOST})',
    )
    args = parser.parse_args()

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    output_path = args.output or Path.cwd() / f'mito_map_annotations_{datetime.now(tz=timezone.utc):%Y%m%d}.json'
    mito_map_host = args.mito_map_host.rstrip('/')
    logger.info(f'Downloading MitoMap annotations from {mito_map_host} to {output_path}')
    download_annotations(output_path, mito_map_host)


if __name__ == '__main__':
    cli_main()
