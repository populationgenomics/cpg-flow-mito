"""
Replace MitoMap annotations made against microprotein ORFs with consequences on their parent gene.

MitoMap assigns some variants to microprotein-derived open reading frames (e.g. MT-SHLP1, MT-GAU) which overlap
well-established genes (e.g. MT-RNR2, MT-CO1). Downstream these variants lose their canonical gene context.

Rather than re-implementing consequence prediction, this is done in three steps:

    1. `to-vcf`: find every annotation on a microprotein named in the mapping, and write each distinct allele to a
        sites-only VCF. The VCF ID column carries the MitoMap position:ref:alt key, used to join the results back
    2. `bcftools csq` annotates that VCF, using the vertebrate mitochondrial codon table (-C 2) - run externally
    3. `apply`: read the BCSQ consequences from the annotated VCF, and wherever there is a consequence on the
        parent gene, overwrite the locus and amino acid fields of the original annotation. Annotations without a
        consequence on the parent gene are left unchanged

The mapping of {microprotein: parent gene(s)} is read from config, at mito_references.microprotein_parents, e.g.
"MT-GAU" = "MT-CO1". A microprotein spanning several genes can have a list of parents, e.g.
"MT-SHMOOSE" = ["MT-ND5", "MT-TS2", "MT-TL2"]. Where an allele has
consequences on more than one of these, the earliest in the list is used.
"""

import gzip
import json
import re
from argparse import ArgumentParser
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from cpg_utils import config
from loguru import logger

# the fields in each BCSQ entry, see bcftools csq documentation
BCSQ_FIELDS = ['consequence', 'gene', 'transcript', 'biotype', 'strand', 'amino_acid_change', 'dna_change']

# matches simple csq protein changes, e.g. 2F>2L (missense), 29A (synonymous), 58V>58VV (inframe insertion)
SIMPLE_AA_CHANGE = re.compile(r'^(\d+)([A-Z*]+)(?:>\d+([A-Z*]+))?$')


@dataclass
class GeneModel:
    biotype: str
    cds_start: int | None = None
    cds_end: int | None = None
    strand: str | None = None


def open_maybe_gzipped(path: Path):
    """Open a text file, transparently decompressing if it is gzipped/bgzipped."""
    with open(path, 'rb') as handle:
        is_gzipped = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if is_gzipped else open(path)


def load_mapping() -> dict[str, list[str]]:
    """Read the {microprotein: parent gene(s)} mapping from config, with each value as a priority-ordered list."""
    mapping = config.config_retrieve(['mito_references', 'microprotein_parents'])
    return {protein: [parents] if isinstance(parents, str) else parents for protein, parents in mapping.items()}


def load_reference(fasta_path: Path) -> tuple[str, str]:
    """Read a single-contig FASTA, returning the contig name and sequence."""
    contig = None
    sequence = []
    with open_maybe_gzipped(fasta_path) as handle:
        for line in handle:
            if line.startswith('>'):
                if contig is not None:
                    raise ValueError(f'Expected a single contig in {fasta_path}')
                contig = line[1:].split()[0]
            else:
                sequence.append(line.strip().upper())
    if contig is None:
        raise ValueError(f'No contig found in {fasta_path}')
    return contig, ''.join(sequence)


def annotation_key(annotation: dict) -> str:
    return f'{annotation["position"]}:{annotation["refAllele"]}:{annotation["altAllele"]}'


def build_locus_anchor(gene_name: str) -> str:
    slug = gene_name.replace('-', '')
    return f'<a href="/MITOMAP/GenomeLoci#{slug}">{gene_name}</a>'


def to_vcf_alleles(annotation: dict, sequence: str) -> tuple[int, str, str] | None:
    """
    Convert a MitoMap allele to a VCF representation, validated against the reference.

    MitoMap deletions are given as the deleted bases with altAllele 'del', at the position of the first deleted base.
    In VCF these are left-anchored on the preceding reference base. All other alleles are already VCF-compatible.

    Returns:
        (POS, REF, ALT), or None if the allele can't be represented or doesn't match the reference
    """
    position = annotation.get('position')
    ref = (annotation.get('refAllele') or '').upper()
    alt = (annotation.get('altAllele') or '').upper()
    if not position or not ref or not alt:
        return None

    # 0-based index of the first base in ref
    start = position - 1
    if sequence[start : start + len(ref)] != ref:
        logger.warning(f'Reference mismatch for {annotation_key(annotation)}, skipping')
        return None

    if alt == 'DEL':
        if start == 0:
            return None
        anchor = sequence[start - 1]
        return position - 1, anchor + ref, anchor

    if not set(alt).issubset('ACGT'):
        logger.warning(f'Unsupported alt allele for {annotation_key(annotation)}, skipping')
        return None

    return position, ref, alt


def write_sites_vcf(annotations: list[dict], mapping: dict[str, list[str]], fasta: Path, output: Path) -> int:
    """
    Write every distinct allele annotated on a microprotein in the mapping to a sites-only VCF.

    Returns:
        the number of VCF records written
    """
    contig, sequence = load_reference(fasta)

    records: dict[str, tuple[int, str, str]] = {}
    for annotation in annotations:
        if annotation.get('locus') not in mapping:
            continue
        key = annotation_key(annotation)
        if key in records:
            continue
        if vcf_alleles := to_vcf_alleles(annotation, sequence):
            records[key] = vcf_alleles

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w') as handle:
        handle.write('##fileformat=VCFv4.2\n')
        handle.write(f'##contig=<ID={contig},length={len(sequence)}>\n')
        handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for key, (pos, ref, alt) in sorted(records.items(), key=lambda item: item[1]):
            handle.write(f'{contig}\t{pos}\t{key}\t{ref}\t{alt}\t.\t.\t.\n')

    return len(records)


def parse_gene_models(gff3: Path) -> dict[str, GeneModel]:
    """
    Get the biotype, and CDS span if protein-coding, of each mitochondrial gene from an Ensembl GFF3.

    Biotypes are taken from here rather than from BCSQ, as some bcftools versions (e.g. 1.22) misreport the biotype of
    non-coding transcripts. Mitochondrial genes are single-exon, so each CDS is taken as a single contiguous span.
    """
    genes: dict[str, GeneModel] = {}
    gene_names: dict[str, str] = {}
    transcript_genes: dict[str, str] = {}

    with open_maybe_gzipped(gff3) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            chrom, _source, feature, start, end, _score, strand, _phase, attributes = line.rstrip('\n').split('\t')
            if chrom not in {'M', 'MT', 'chrM'}:
                continue
            attrs = dict(attr.split('=', 1) for attr in attributes.split(';') if '=' in attr)
            if attrs.get('ID', '').startswith('gene:') and 'Name' in attrs:
                gene_names[attrs['ID']] = attrs['Name']
                genes[attrs['Name']] = GeneModel(biotype=attrs.get('biotype', ''))
            elif attrs.get('ID', '').startswith('transcript:'):
                transcript_genes[attrs['ID']] = attrs['Parent']
            elif feature == 'CDS':
                gene = genes[gene_names[transcript_genes[attrs['Parent']]]]
                cds_start, cds_end = gene.cds_start or int(start), gene.cds_end or int(end)
                gene.cds_start, gene.cds_end = min(cds_start, int(start)), max(cds_end, int(end))
                gene.strand = strand

    return genes


def parse_bcsq(vcf: Path) -> dict[str, list[dict[str, str]]]:
    """
    Read the BCSQ annotations from a bcftools csq annotated VCF.

    Returns:
        {MitoMap key (from the ID column): [each consequence, as a dict of BCSQ_FIELDS]}
    """
    consequences: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open_maybe_gzipped(vcf) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            key, info = fields[2], fields[7]
            for entry in info.split(';'):
                if not entry.startswith('BCSQ='):
                    continue
                for csq in entry.removeprefix('BCSQ=').split(','):
                    # @POS entries are references to a compound consequence reported at another position
                    if csq.startswith('@'):
                        continue
                    consequences[key].append(dict(zip(BCSQ_FIELDS, csq.split('|'), strict=False)))
    return consequences


def format_amino_acid_change(csq: dict[str, str]) -> str:
    """Convert a csq protein change to the MitoMap style, e.g. 2F>2L -> F-L, 29A -> A-A, 5W>5* -> W-Term."""
    if 'frameshift' in csq['consequence']:
        return 'frameshift'
    aa_change = csq.get('amino_acid_change', '')
    if match := SIMPLE_AA_CHANGE.match(aa_change):
        ref_aa, alt_aa = match.group(2), match.group(3) or match.group(2)
        return f'{ref_aa}-{alt_aa}'.replace('*', 'Term')
    return aa_change or csq['consequence']


def replacement_fields(csq: dict[str, str], position: int, gene_model: GeneModel) -> dict[str, str]:
    """Generate the replacement annotation fields for a MitoMap annotation, from a csq consequence on the parent."""
    gene = csq['gene']
    fields = {'locus': gene, 'locusAnchor': build_locus_anchor(gene)}

    if gene_model.cds_start is None or gene_model.cds_end is None:
        biotype = gene_model.biotype.lower()
        fields['aminoAcidChange'] = 'rRNA' if 'rrna' in biotype else 'tRNA' if 'trna' in biotype else gene_model.biotype
        fields['codonNumber'] = '-'
        fields['codonPosition'] = '-'
        return fields

    # codon number & position follow the MitoMap convention: the first affected base in the CDS
    offset = position - gene_model.cds_start if gene_model.strand == '+' else gene_model.cds_end - position
    fields['aminoAcidChange'] = format_amino_acid_change(csq)
    fields['codonNumber'] = str(offset // 3 + 1)
    fields['codonPosition'] = str(offset % 3 + 1)
    return fields


def apply_replacements(
    annotations: list[dict],
    mapping: dict[str, list[str]],
    consequences: dict[str, list[dict[str, str]]],
    gene_models: dict[str, GeneModel],
) -> tuple[Counter, Counter]:
    """
    Update each microprotein annotation in place, where a consequence on one of its parent genes is available.

    Parents are checked in the order given in the mapping, and the first with a consequence is used.

    Returns:
        Counters of replaced annotations keyed on (microprotein, parent), and unresolved keyed on the microprotein
    """
    replaced: Counter = Counter()
    unresolved: Counter = Counter()

    for annotation in annotations:
        microprotein = annotation.get('locus')
        if microprotein not in mapping:
            continue

        csq_by_gene: dict[str, dict[str, str]] = {}
        for csq in consequences.get(annotation_key(annotation), []):
            csq_by_gene.setdefault(csq.get('gene', ''), csq)

        parent = next((gene for gene in mapping[microprotein] if gene in csq_by_gene and gene in gene_models), None)
        if parent is None:
            unresolved[microprotein] += 1
            continue

        annotation.update(replacement_fields(csq_by_gene[parent], annotation['position'], gene_models[parent]))
        replaced[(microprotein, parent)] += 1

    return replaced, unresolved


def to_vcf_main(args) -> None:
    with open(args.input) as handle:
        annotations = json.load(handle)
    mapping = load_mapping()
    count = write_sites_vcf(annotations, mapping, args.reference, args.output)
    logger.info(f'Wrote {count} candidate alleles for {len(mapping)} microproteins to {args.output}')


def apply_main(args) -> None:
    with open(args.input) as handle:
        annotations = json.load(handle)
    mapping = load_mapping()
    consequences = parse_bcsq(args.vcf)
    gene_models = parse_gene_models(args.gff3)

    replaced, unresolved = apply_replacements(annotations, mapping, consequences, gene_models)

    for microprotein, parents in mapping.items():
        counts = ', '.join(f'{replaced[(microprotein, parent)]} -> {parent}' for parent in parents)
        logger.info(f'{microprotein}: {counts}, {unresolved[microprotein]} unresolved')

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as handle:
        json.dump(annotations, handle, indent=2)


def cli_main() -> None:
    parser = ArgumentParser(description='Replace MitoMap microprotein annotations with parent gene consequences')
    subparsers = parser.add_subparsers(required=True)

    to_vcf = subparsers.add_parser('to-vcf', help='Write candidate alleles to a sites-only VCF')
    to_vcf.add_argument('-i', '--input', type=Path, required=True, help='Input MitoMap annotations JSON')
    to_vcf.add_argument('-r', '--reference', type=Path, required=True, help='chrM FASTA')
    to_vcf.add_argument('-o', '--output', type=Path, required=True, help='Output sites-only VCF')
    to_vcf.set_defaults(func=to_vcf_main)

    apply = subparsers.add_parser('apply', help='Apply csq consequences on parent genes to the annotations')
    apply.add_argument('-i', '--input', type=Path, required=True, help='Input MitoMap annotations JSON')
    apply.add_argument('-v', '--vcf', type=Path, required=True, help='bcftools csq annotated VCF, from to-vcf output')
    apply.add_argument('-g', '--gff3', type=Path, required=True, help='GFF3 used for csq, for biotypes & CDS')
    apply.add_argument('-o', '--output', type=Path, required=True, help='Output annotations JSON')
    apply.set_defaults(func=apply_main)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    cli_main()
