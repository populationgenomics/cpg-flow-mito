"""
Remap microprotein loci to their canonical gene counterparts in MitoMap annotation JSON.

MitoMap assigns some variants to microprotein-derived open reading frames (e.g. MT-SHLP1,
MT-SHMOOSE, MT-TER) that overlap well-established genes (e.g. MT-RNR2, MT-TL1). This
script remaps those loci so that downstream tools see variants attributed to the canonical
gene, making them easier to interpret clinically.

The remapping rules are defined in a CSV with columns:
    microprotein, micro_start, micro_end, canonical_replacement, canonical_start, canonical_end, ...

A variant is remapped when its locusAnchor matches the microprotein name AND its position
falls within the canonical gene's range. Variants outside the canonical range (e.g. MT-OLR
entries outside MT-TC) are left unchanged.

For protein-coding canonical genes, codonNumber and codonPosition are recomputed from the
canonical gene's CDS boundaries using the mitochondrial genetic code (NCBI table 2).
aminoAcidChange is preserved for SNVs (verified to be frame-identical) and cleared for indels.
For tRNA/rRNA canonical genes, aminoAcidChange is set to 'tRNA'.
"""

import csv
import json
import re
import sys
from argparse import ArgumentParser
from importlib import resources
from pathlib import Path

MITO_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': '*', 'AGG': '*',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

PROTEIN_CODING_CDS = {
    'MT-ND1': (3307, 4262), 'MT-ND2': (4470, 5511),
    'MT-CO1': (5904, 7445), 'MT-CO2': (7586, 8269),
    'MT-ATP8': (8366, 8572), 'MT-ATP6': (8527, 9207),
    'MT-CO3': (9207, 9990), 'MT-ND3': (10059, 10404),
    'MT-ND4L': (10470, 10766), 'MT-ND4': (10760, 12137),
    'MT-ND5': (12337, 14148), 'MT-CYB': (14747, 15887),
}

RRNA_GENES = {'MT-RNR1', 'MT-RNR2'}

TRNA_GENES = {
    'MT-TF', 'MT-TV', 'MT-TL1', 'MT-TI', 'MT-TQ', 'MT-TM',
    'MT-TW', 'MT-TA', 'MT-TN', 'MT-TC', 'MT-TY', 'MT-TS1',
    'MT-TD', 'MT-TK', 'MT-TG', 'MT-TR', 'MT-TH', 'MT-TS2',
    'MT-TL2', 'MT-TE', 'MT-TT', 'MT-TP',
}

NON_CODING_GENES = RRNA_GENES | TRNA_GENES


def load_remapping_rules(csv_path: Path | None = None) -> list[dict]:
    if csv_path is not None:
        text = csv_path.read_text()
    else:
        ref = resources.files('cpg_flow_mito.data').joinpath('microprotein_remapping.csv')
        text = ref.read_text(encoding='utf-8')
    return list(csv.DictReader(text.splitlines()))


def load_reference(fasta_path: Path) -> str:
    with open(fasta_path) as f:
        return ''.join(line.strip() for line in f if not line.startswith('>'))


def build_locus_anchor(gene_name: str) -> str:
    slug = gene_name.replace('-', '')
    return f'<a href="/MITOMAP/GenomeLoci#{slug}">{gene_name}</a>'


def _is_snv(ref_allele: str | None, alt_allele: str | None) -> bool:
    return bool(ref_allele and alt_allele and len(ref_allele) == 1 and len(alt_allele) == 1 and alt_allele != 'del')


def _normalize_aa_change(aa_change: str | None) -> str | None:
    """Normalize formats like 'W163R' -> 'W-R' and 'W-Term' -> 'W-*'."""
    if not aa_change:
        return aa_change
    aa_change = aa_change.replace('Term', '*')
    m = re.match(r'^([A-Z*])(\d+)([A-Z*])$', aa_change)
    if m:
        return f'{m.group(1)}-{m.group(3)}'
    return aa_change


def compute_codon_info(position: int, alt_allele: str, gene: str, ref_seq: str) -> dict | None:
    """Compute codonNumber, codonPosition, and aminoAcidChange for a SNV in a protein-coding gene."""
    if gene not in PROTEIN_CODING_CDS:
        return None
    cds_start, cds_end = PROTEIN_CODING_CDS[gene]
    if not (cds_start <= position <= cds_end):
        return None

    offset = position - cds_start
    codon_number = (offset // 3) + 1
    codon_position = (offset % 3) + 1
    codon_start_0 = cds_start - 1 + (codon_number - 1) * 3
    ref_codon = ref_seq[codon_start_0:codon_start_0 + 3].upper()

    alt_codon = list(ref_codon)
    alt_codon[codon_position - 1] = alt_allele.upper()
    alt_codon = ''.join(alt_codon)

    ref_aa = MITO_CODON_TABLE.get(ref_codon, '?')
    alt_aa = MITO_CODON_TABLE.get(alt_codon, '?')

    return {
        'aminoAcidChange': f'{ref_aa}-{alt_aa}',
        'codonNumber': str(codon_number),
        'codonPosition': str(codon_position),
    }


def remap_annotations(
    annotations: list[dict],
    rules: list[dict],
    ref_seq: str | None = None,
) -> tuple[list[dict], int]:
    lookup: dict[str, list[tuple[int, int, str]]] = {}
    for rule in rules:
        mp = rule['microprotein']
        lookup.setdefault(mp, []).append((
            int(rule['canonical_start']),
            int(rule['canonical_end']),
            rule['canonical_replacement'],
        ))

    remapped_count = 0
    for annotation in annotations:
        locus = annotation.get('locus')
        if not locus or locus not in lookup:
            continue

        pos = annotation.get('position', 0)
        for c_start, c_end, canonical in lookup[locus]:
            if c_start <= pos <= c_end:
                annotation['locusAnchor'] = build_locus_anchor(canonical)
                annotation['locus'] = canonical

                # MitoMap catalogued these under CODING for the microprotein ORF; correct to
                # RNA_MUTATIONS when the canonical gene is rRNA/tRNA so the label stays accurate
                # even though mitoreport v1.1.0 does not currently read regionType.
                if canonical in NON_CODING_GENES and annotation.get('regionType') == 'CODING':
                    annotation['regionType'] = 'RNA_MUTATIONS'

                ref_allele = annotation.get('refAllele')
                alt_allele = annotation.get('altAllele')

                if canonical in NON_CODING_GENES:
                    annotation['aminoAcidChange'] = 'rRNA' if canonical in RRNA_GENES else 'tRNA'
                    annotation['codonNumber'] = '-'
                    annotation['codonPosition'] = '-'
                elif _is_snv(ref_allele, alt_allele) and ref_seq:
                    codon_info = compute_codon_info(pos, alt_allele, canonical, ref_seq)
                    if codon_info:
                        annotation['aminoAcidChange'] = codon_info['aminoAcidChange']
                        annotation['codonNumber'] = codon_info['codonNumber']
                        annotation['codonPosition'] = codon_info['codonPosition']
                    else:
                        annotation['aminoAcidChange'] = _normalize_aa_change(annotation.get('aminoAcidChange'))
                else:
                    annotation['aminoAcidChange'] = None
                    annotation['codonNumber'] = '-'
                    annotation['codonPosition'] = '-'

                remapped_count += 1
                break

    return annotations, remapped_count


def cli_main() -> None:
    parser = ArgumentParser(description='Remap microprotein loci to canonical genes in MitoMap annotations')
    parser.add_argument('-i', '--input', type=Path, required=True, help='Input annotations JSON')
    parser.add_argument('-o', '--output', type=Path, required=True, help='Output annotations JSON')
    parser.add_argument('-c', '--csv', type=Path, default=None, help='Remapping CSV (default: bundled)')
    parser.add_argument('-r', '--reference', type=Path, default=None, help='chrM FASTA for codon recomputation')
    args = parser.parse_args()

    with open(args.input) as f:
        annotations = json.load(f)

    rules = load_remapping_rules(args.csv)
    ref_seq = load_reference(args.reference) if args.reference else None

    annotations, count = remap_annotations(annotations, rules, ref_seq)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(annotations, f, indent=2)

    print(f'Remapped {count} annotations across {len(rules)} rules', file=sys.stderr)


if __name__ == '__main__':
    cli_main()
