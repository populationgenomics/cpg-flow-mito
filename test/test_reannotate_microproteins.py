from cpg_flow_mito.scripts import reannotate_microproteins
from cpg_flow_mito.scripts.reannotate_microproteins import (
    GeneModel,
    apply_replacements,
    format_amino_acid_change,
    load_mapping,
    parse_bcsq,
    parse_gene_models,
    to_vcf_alleles,
    write_sites_vcf,
)

# positions 1-12 of a toy chrM
SEQUENCE = 'GATCACAGGTCT'

GFF3 = '\n'.join(
    [
        '##gff-version 3',
        'M\tinsdc\tgene\t2\t10\t.\t+\t.\tID=gene:G1;Name=MT-CO1;biotype=protein_coding',
        'M\tinsdc\tmRNA\t2\t10\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding',
        'M\tinsdc\tCDS\t2\t10\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1',
        'M\tinsdc\tncRNA_gene\t11\t12\t.\t+\t.\tID=gene:G2;Name=MT-RNR2;biotype=Mt_rRNA',
        'M\tinsdc\tMt_rRNA\t11\t12\t.\t+\t.\tID=transcript:T2;Parent=gene:G2;biotype=Mt_rRNA',
        '1\tinsdc\tgene\t2\t10\t.\t+\t.\tID=gene:G3;Name=OTHER;biotype=protein_coding',
    ],
)


def annotation(position: int, ref: str, alt: str, locus: str = 'MT-GAU') -> dict:
    return {
        'position': position,
        'refAllele': ref,
        'altAllele': alt,
        'locus': locus,
        'locusAnchor': 'original',
        'aminoAcidChange': 'original',
        'codonNumber': '-',
        'codonPosition': '-',
    }


def test_to_vcf_alleles_snv_and_insertion():
    assert to_vcf_alleles(annotation(4, 'C', 'T'), SEQUENCE) == (4, 'C', 'T')
    assert to_vcf_alleles(annotation(4, 'C', 'CC'), SEQUENCE) == (4, 'C', 'CC')


def test_to_vcf_alleles_deletion_is_left_anchored():
    # MitoMap m.4CAdel -> VCF 3 TCA>T
    assert to_vcf_alleles(annotation(4, 'CA', 'del'), SEQUENCE) == (3, 'TCA', 'T')


def test_to_vcf_alleles_rejects_reference_mismatch_and_missing():
    assert to_vcf_alleles(annotation(4, 'G', 'T'), SEQUENCE) is None
    assert to_vcf_alleles(annotation(4, None, None), SEQUENCE) is None


def test_write_sites_vcf_only_candidates_deduplicated(tmp_path):
    fasta = tmp_path / 'chrM.fa'
    fasta.write_text(f'>chrM\n{SEQUENCE}\n')
    output = tmp_path / 'out.vcf'
    annotations = [
        annotation(5, 'A', 'G'),
        annotation(5, 'A', 'G', locus='MT-SHLP1'),
        annotation(4, 'C', 'T', locus='MT-CO1'),
    ]
    count = write_sites_vcf(annotations, {'MT-GAU': 'MT-CO1', 'MT-SHLP1': 'MT-RNR2'}, fasta, output)
    assert count == 1
    records = [line.split('\t') for line in output.read_text().splitlines() if not line.startswith('#')]
    assert records == [['chrM', '5', '5:A:G', 'A', 'G', '.', '.', '.']]
    assert '##contig=<ID=chrM,length=12>' in output.read_text()


def test_parse_gene_models(tmp_path):
    gff3 = tmp_path / 'test.gff3'
    gff3.write_text(GFF3)
    models = parse_gene_models(gff3)
    assert models == {
        'MT-CO1': GeneModel(biotype='protein_coding', cds_start=2, cds_end=10, strand='+'),
        'MT-RNR2': GeneModel(biotype='Mt_rRNA'),
    }


def test_format_amino_acid_change():
    assert format_amino_acid_change({'consequence': 'missense', 'amino_acid_change': '2F>2L'}) == 'F-L'
    assert format_amino_acid_change({'consequence': 'synonymous', 'amino_acid_change': '29A'}) == 'A-A'
    assert format_amino_acid_change({'consequence': 'stop_gained', 'amino_acid_change': '31W>31*'}) == 'W-Term'
    assert format_amino_acid_change({'consequence': 'inframe_insertion', 'amino_acid_change': '58V>58VV'}) == 'V-VV'
    assert format_amino_acid_change({'consequence': 'stop_gained&frameshift', 'amino_acid_change': '1PLG>1P*'}) == (
        'frameshift'
    )


def test_parse_bcsq_and_apply(tmp_path):
    vcf = tmp_path / 'annotated.vcf'
    vcf.write_text(
        '\n'.join(
            [
                '##fileformat=VCFv4.2',
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
                'chrM\t5\t5:A:G\tA\tG\t.\t.\tBCSQ=missense|MT-CO1|T1|protein_coding|+|2T>2A|5A>G',
                # csq 1.22 misreports this biotype, which should be ignored in favour of the GFF3
                'chrM\t11\t11:C:T\tC\tT\t.\t.\tBCSQ=non_coding|MT-RNR2||MT_tRNA',
                'chrM\t12\t12:T:C\tT\tC\t.\t.\t.',
            ],
        ),
    )
    gff3 = tmp_path / 'test.gff3'
    gff3.write_text(GFF3)

    annotations = [
        annotation(5, 'A', 'G'),
        annotation(11, 'C', 'T', locus='MT-SHLP1'),
        annotation(12, 'T', 'C', locus='MT-SHLP1'),
        annotation(5, 'A', 'G', locus='MT-CO1'),
    ]
    mapping = {'MT-GAU': ['MT-CO1'], 'MT-SHLP1': ['MT-RNR2']}
    replaced, unresolved = apply_replacements(annotations, mapping, parse_bcsq(vcf), parse_gene_models(gff3))

    assert replaced == {('MT-GAU', 'MT-CO1'): 1, ('MT-SHLP1', 'MT-RNR2'): 1}
    assert unresolved == {'MT-SHLP1': 1}

    coding, rrna, no_csq, untouched = annotations
    assert coding['locus'] == 'MT-CO1'
    assert coding['locusAnchor'] == '<a href="/MITOMAP/GenomeLoci#MTCO1">MT-CO1</a>'
    assert (coding['aminoAcidChange'], coding['codonNumber'], coding['codonPosition']) == ('T-A', '2', '1')
    assert rrna['locus'] == 'MT-RNR2'
    assert (rrna['aminoAcidChange'], rrna['codonNumber'], rrna['codonPosition']) == ('rRNA', '-', '-')
    assert no_csq == annotation(12, 'T', 'C', locus='MT-SHLP1')
    assert untouched == annotation(5, 'A', 'G', locus='MT-CO1')


def test_load_mapping_normalises_to_lists(monkeypatch):
    mapping = {'MT-GAU': 'MT-CO1', 'MT-SHMOOSE': ['MT-ND5', 'MT-TS2']}
    monkeypatch.setattr(reannotate_microproteins.config, 'config_retrieve', lambda _key: mapping)
    assert load_mapping() == {'MT-GAU': ['MT-CO1'], 'MT-SHMOOSE': ['MT-ND5', 'MT-TS2']}


def test_apply_multiple_parents_uses_first_listed(tmp_path):
    vcf = tmp_path / 'annotated.vcf'
    vcf.write_text(
        '\n'.join(
            [
                '##fileformat=VCFv4.2',
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
                # consequences on both parents, the first listed parent (MT-CO1) should be used
                'chrM\t5\t5:A:G\tA\tG\t.\t.\tBCSQ=non_coding|MT-RNR2||Mt_rRNA,missense|MT-CO1|T1|protein_coding|+|2T>2A|5A>G',
                # only a consequence on the second listed parent
                'chrM\t11\t11:C:T\tC\tT\t.\t.\tBCSQ=non_coding|MT-RNR2||Mt_rRNA',
            ],
        ),
    )
    gff3 = tmp_path / 'test.gff3'
    gff3.write_text(GFF3)

    annotations = [annotation(5, 'A', 'G', locus='MT-SHMOOSE'), annotation(11, 'C', 'T', locus='MT-SHMOOSE')]
    mapping = {'MT-SHMOOSE': ['MT-CO1', 'MT-RNR2']}
    replaced, unresolved = apply_replacements(annotations, mapping, parse_bcsq(vcf), parse_gene_models(gff3))

    assert replaced == {('MT-SHMOOSE', 'MT-CO1'): 1, ('MT-SHMOOSE', 'MT-RNR2'): 1}
    assert not unresolved
    assert annotations[0]['locus'] == 'MT-CO1'
    assert annotations[0]['aminoAcidChange'] == 'T-A'
    assert annotations[1]['locus'] == 'MT-RNR2'
    assert annotations[1]['aminoAcidChange'] == 'rRNA'
