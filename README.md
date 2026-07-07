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
