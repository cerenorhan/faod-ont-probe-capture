# FAOD ONT Probe-Capture Workflow

Reproducibility repository accompanying the manuscript evaluating an established amplicon-based long-read PCR/ONT workflow and a custom probe-capture ONT workflow for molecular analysis of fatty acid oxidation disorders.

## Repository contents

- `pipelines/`: final V2 LR-PCR/ONT and V3 probe-capture/ONT processing pipelines
- `scripts/`: custom scripts used for analyses reported in the manuscript
- `data/`: publication-derived tab-delimited datasets underlying reported analyses
- `qc_reports/`: workflow-level sequencing quality-control reports for the V2 and V3 datasets
- `figures/`: computationally generated figures reported in the manuscript or Supplementary Information

Only publication-relevant derived datasets, workflow-level quality-control reports, analytical outputs, and final processing pipelines are included.

Individual-level raw sequencing data, BAM files, VCF files, and other potentially identifiable human genomic data are not publicly distributed.

## Data scope

The repository includes derived data supporting:

- cohort and paired-sample selection
- sequencing and workflow quality control
- shared-target coverage and interval-level coverage rescue
- clinically relevant variant comparison
- allele-frequency agreement and genotype-consistency analyses
- candidate-call quality profiling
- technical reproducibility analyses

## Pipelines

The repository includes the final clinical-status v3.2 processing pipelines used for the V2 LR-PCR/ONT and V3 probe-capture/ONT workflows, including the BAM-level safety-net interrogation framework.

## Citation

A versioned release of this repository will be archived in Zenodo and assigned a DOI.
