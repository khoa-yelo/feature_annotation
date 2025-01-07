# Feature Annotation

## Overview

## Method

### BLAST

### Lookup table 

### Inputs

### Outputs

## Usage

1. Load Nextflow

```
ml biology load nextflow/23.04.3
```

2. Set input parameter at `configs/params.yaml`
```
# output directory
outdir: "/oak/stanford/groups/horence/khoa/scratch/repos/feature_annotation/test_data/outs"
# bin contain scripts that run processes
bin: "/oak/stanford/groups/horence/khoa/scratch/repos/feature_annotation/workflows/bin"

blast: True # whether to use bast for annotation
lookup_table: True # whether to use lookup table for annotation
# path to lookup tables used for query (path1:path2:path3...)
lookup_table_paths: '/oak/stanford/groups/horence/khoa/shares/lookup_tables/artifacts.slt:/oak/stanford/groups/horence/khoa/shares/lookup_tables/common_microbe_w_trnx.slt'

# input_type either fasta or tsv
input_type: "tsv"
# which column in tsv is used for annotation
annotate_column: "anchor"
# input: either tsv or fasta
input: /oak/stanford/groups/horence/khoa/scratch/repos/feature_annotation/test_data/test_inputs/seqs.tsv

```
3. Run pipeline
```
nextflow run workflows/feature_annotate.nf -c configs/run.config -params-file configs/params.yml -resume
```

## Example

