# cpg-flow-mito

A re-implementation of the CPG's RD Mitochondrial Analysis workflow, using cpg-flow. The original CPG-workflows implementation is [here](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/stages/mito.py)

The production-pipelines version of this workflow was itself a re-implementation of the original Broad workflow,
originally created in WDL [here](https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl)

```text
src
├── __init__.py
└── cpg_flow_mito
    ├── __init__.py
    ├── config_template.toml
    ├── jobs
    │   ├── __init__.py
    │   ├── mito.py
    │   ├── picard.py
    │   └── vep.py
    ├── run_workflow.py
    ├── scripts
    │   └── __init__.py
    ├── stages.py
    └── utils.py
```

## Error handling

In the event of variant calling on a Sequencing Group generating no calls, or exclusively filtered calls, the contamination
sub-workflow fill fail, derailing the whole run. As a solution the config setting `workflow.skip_contamination` can be
used to omit this extra analysis.

Recommendation
-
- on the first pass, always include the contamination workflow
- in the event of recurrent failure at the parse_contamination_results step, run with the skip_contamination flag

All samples with results already will be unaffected, and remaining samples will be run without contamination metrics.

Example invocation:

```bash
analysis-runner \
    --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-mito:0.3.1-1\
    --config your laptop path/cpg-flow-mito/src/cpg_flow_mito/config_template.toml\
    --dataset seqr \
    --description 'mitoindex' \
    --access-level full \
    --output-dir mitoindex\
  python3 src/cpg_flow_mito/run_workflow.py
```
The config file should be modifed to have these at minimum:

```toml
[workflow]
input_cohorts = ["<cohort_id>"]
sequencing_type = "<exome|genome>"
```

## MitoMap microprotein annotation remapping

MitoMap has begun annotating variants under recently proposed microprotein and alternative ORF gene names
(e.g. MT-GAU, MT-SHMOOSE, MT-CYTB-187AA) that overlap with established canonical mitochondrial genes.
This causes variants in well-characterised coding regions (such as MT-CO1, MT-CYB, MT-ND4) to lose their
codon context and appear as non-coding in downstream reports.

The `RemapAnnotations` stage replaces these with consequences on the parent gene, using `bcftools csq` rather than
re-implementing consequence prediction. The mapping of microprotein to parent gene is set in config, as
`mito_references.microprotein_parents`. Microproteins spanning several genes can list multiple parents, e.g.
`"MT-SHMOOSE" = ["MT-ND5", "MT-TS2", "MT-TL2"]`; where an allele has consequences on more than one, the first listed
is used. The script `scripts/reannotate_microproteins.py` runs in three steps:

1. `to-vcf`: every allele annotated on a microprotein in the mapping is written to a sites-only VCF
2. `bcftools csq` annotates that VCF against `mito_references.gff3`, using the mitochondrial codon table (`-C 2`)
3. `apply`: wherever csq reports a consequence on the parent gene, the `locus`, `locusAnchor`, `aminoAcidChange`,
   `codonNumber`, and `codonPosition` fields are replaced. Annotations without a consequence on any of their parent
   genes are left unchanged
