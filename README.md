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
