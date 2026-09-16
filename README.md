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
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-mito:0.2.4-1\
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

The script `remap_microprotein_annotations.py` remaps these annotations back to the canonical gene, but
only when the variant position falls within the canonical gene's rCRS (NC_012920.1) coordinates.

### Remapping rules

| Microprotein | Canonical gene | Source |
|---|---|---|
| MT-GAU | MT-CO1 | Faure et al. 2011, Biology Direct 6:56 |
| MT-CYTB-187AA | MT-CYB | |
| MT-ALTND4 | MT-ND4 | |
| MT-SHMOOSE | MT-TS2 / MT-TL2 / MT-ND5 | Position-dependent (spans 3 genes) |
| MT-HN, MT-Hum | MT-RNR2 | Humanin |
| MT-MOTSc | MT-RNR1 | |
| MT-SHLP1–6 | MT-RNR2 | |
| MT-TER | MT-TL1 | |
| MT-RNR3 | MT-RNR2 | |
| MT-OLR | MT-TC | Only positions within MT-TC (5761–5826); non-coding OL positions are left unchanged |
