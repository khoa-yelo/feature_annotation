# Feature Annotation

## Overview
Tool to annotate DNA sequence with BLAST and Lookup table.
User can optionally run either BLAST, Lookup table or both.
Lookup table is ran on a list of specified input tables.

Repo structure

- configs: folder containing input params and Nextflow run config
- containers: stores Docker image used in pipeline
- modules: contains Processes in pipeline
- test_data: test input data and test outs
- workflows: pipeline definition and scripts for each Process in `bin` directory

## Method

### Inputs

- Either `fasta` or `tsv` file with specified column used for annotation (e.g. anchor)
- A list of built lookup tables to query from

Input params should be specified in `configs/params.yml`

### Outputs

```
.
├── blast_features_outs
│   # contain blast feature output of each splitted chunks
├── blast_outs
│   # contain blast output of each splitted chunks
├── lookup_out
│   # lookup table output for all specified tables
├── preprocessed_inputs
│   ├── input_chunks
│   # splitted fasta files
│   └── input.fasta # in case of tsv input, also output fasta of specified column
└── summary
    └── annotation.tsv # summarized annotation file
```

Details of outputs in each folder below


### Process: Preprocess

This process reads in input file and split the input into chunks, each chunk contains 50 sequences, for parallel processing. Output recored in `preprocessed_inputs` folder, containing splitted fasta and a copy of the input fasta file or the fasta converted version of the specified column in case of tsv input. See `workflows/bin/preprocess.py`

### Process: Lookup_table

This process runs lookup tables directly on the input with the following cmd

```
lookup_table query --truncate_paths --stats_fmt with_stats {input_file} {lookup_table} {output_file}
```
Output recored in `lookup_out` folder. 3 output will be generated for each table
- [table_name].lookup_out.tsv: direct output of lookup table
- [table_name].lookup_out_cleaned.tsv: cleaned version of `[table_name].lookup_out.tsv`

    This file contain 3 columns: `sequence`, `stat`, and `matches`
    - sequence: input sequence
    - stat: count of each category
    - matches: the most abundance annotation within each category


- [table_name].lookup_out_all.tsv: merged of `[table_name].lookup_out.tsv` and `[table_name].lookup_out_cleaned.tsv`

See `workflows/bin/lookup_table.py`


### Process: Blast

BLAST is ran with the following command on each splitted chunks and outputs are recorded at `blast_outs`

```
fmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseqid sgi sacc slen staxids stitle"

cmd = f"blastn -outfmt '{fmt}' -query {input_fasta} -remote -db nt -out {blast_out} -evalue 0.1 -task blastn -dust no -word_size 24 -reward 1 -penalty -3 -max_target_seqs 4"

```

BLAST output columns includes:
```
["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",\
 "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "sgi", \
 "sacc", "slen", "staxids", "stitle"]
```

See `workflows/bin/blast.py`


### Process: Blast extract feature

This process extracts annotations of subject sequence overlapped or within 5000 nucleotide of query sequence using `Bio.Entrex`. These are features such as CDS, genes, amino acid sequences. See `workflows/bin/blast_extract_features.py`

Outputs are recorded at `blast_features_outs`
Output file contains 4 columns:
- query: input sequence
- identity: identity to subject sequence	
- features: annotation from Entrez of the overlapped region
- features_5000_window: annotaions within 5000 nt of the overlapped region


### Process: Summarize

This process merge all outputs of lookup tables (use *.lookup_out_all.tsv files), blast, and blast features into a single `annotation.tsv` file recorded at `summary`. See `workflows/bin/summarize.py`

In case of multiple hits, subject sequence with highest identity is kept.


## Usage

1. Load Nextflow

```
ml biology load nextflow/23.04.3
```

2. Set input parameter at `configs/params.yaml`

3. Run pipeline
```
nextflow run workflows/feature_annotate.nf -c configs/run.config -params-file configs/params.yml -resume
```

## Example

In this example, input params is set to run on 200 anchors in `seqs.tsv` in `test_data`. 200 sequences will be spliited into 4 chunks for parallel processing.

2 look up tables is specified:
- `/oak/stanford/groups/horence/khoa/shares/lookup_tables/artifacts.slt` was built from `/oak/stanford/groups/horence/khoa/shares/lookup_tables/artifacts.txt`
- `/oak/stanford/groups/horence/khoa/shares/lookup_tables/common_microbe_w_trnx.slt` was built from 
    - `/oak/stanford/groups/horence/khoa/shares/lookup_tables/microbial_list.txt` 
    - `/oak/stanford/groups/horence/khoa/shares/lookup_tables/microbes_transcriptomes.txt`(transcriptomes)




`params.yml` file:

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
Run log

```
(base) [khoang99@sh03-16n22 /oak/stanford/groups/horence/khoa/scratch/repos/feature_annotation] (job 57825140) $ nextflow run workflows/feature_annotate.nf -c configs/run.config -params-file configs/params.yml -resume
N E X T F L O W  ~  version 23.04.3
Launching `workflows/feature_annotate.nf` [zen_sammet] DSL2 - revision: ed3520c699
[d1/de68b3] process > Preprocess (Preprocess null)                [100%] 1 of 1, cached: 1 ✔
[81/0f26ef] process > Lookup_table (Lookup table input.fasta)     [100%] 1 of 1, cached: 1 ✔
[c1/e5ba69] process > Blast (BLAST null)                          [100%] 4 of 4, cached: 4 ✔
[68/128557] process > Blast_extract_features (BLAST feature null) [100%] 4 of 4, cached: 4 ✔
[4a/bff089] process > Summarize (Summarize)                       [100%] 1 of 1, cached: 1 ✔
```

Output 

```
.
├── blast_features_outs
│   ├── split_0.blast_out.features.tsv
│   ├── split_1.blast_out.features.tsv
│   ├── split_2.blast_out.features.tsv
│   └── split_3.blast_out.features.tsv
├── blast_outs
│   ├── split_0.blast_out.tsv
│   ├── split_1.blast_out.tsv
│   ├── split_2.blast_out.tsv
│   └── split_3.blast_out.tsv
├── lookup_out
│   ├── artifacts.lookup_out_all.tsv
│   ├── artifacts.lookup_out_cleaned.tsv
│   ├── artifacts.lookup_out.tsv
│   ├── common_microbe_w_trnx.lookup_out_all.tsv
│   ├── common_microbe_w_trnx.lookup_out_cleaned.tsv
│   └── common_microbe_w_trnx.lookup_out.tsv
├── preprocessed_inputs
│   ├── input_chunks
│   │   ├── split_0.fasta
│   │   ├── split_1.fasta
│   │   ├── split_2.fasta
│   │   └── split_3.fasta
│   └── input.fasta
└── summary
    └── annotation.tsv
```