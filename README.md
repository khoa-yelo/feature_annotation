# Feature Annotation

## Overview
This pipeline annotates DNA sequences using two complementary methods: BLAST and lookup table queries. Users can choose to run either BLAST, lookup table, or both. In addition, a lookup table process supports multiple input lookup tables. 

Repository Structure

- `configs`: Contains input parameters and Nextflow run configuration files.
- `containers`: Stores the Docker image used in the pipeline.
- `modules`: Contains Nextflow processes that execute individual tasks.
- `test_data`: Provides test inputs and example outputs for verification.
- `workflows`: Contains the main pipeline definition and scripts located in the bin directory.


## Method

### Inputs

- A `fasta` or a `tsv` file with specified column used for annotation (e.g. anchor)
- A list of built lookup tables to query against.

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
### What it does:

- Reads the input file (FASTA or TSV) and splits it into chunks of 50 sequences for parallel processing.
- Produces a folder (`preprocessed_inputs`) containing:
  - A directory (`input_chunks`) with individual split FASTA files.
  - A master FASTA file (`input.fasta`), which is either the original FASTA or a converted version of the specified column from the TSV input.
### Script:
```
workflows/bin/preprocess.py
```

### Process: Lookup_table

### What it does:
For each provided lookup table, the process runs a query using:
```
lookup_table query --truncate_paths --stats_fmt with_stats {input_file} {lookup_table} {output_file}
```
It produces three outputs per lookup table in the `lookup_out` folder:
- [table_name].lookup_out.tsv: direct output of lookup table

Further explanation on the meaning of each category can be found [here](https://github.com/refresh-bio/SPLASH/wiki/Lookup_table)

- [table_name].lookup_out_cleaned.tsv: cleaned version of `[table_name].lookup_out.tsv`

    This file contain 3 columns: `sequence`, `stat`, and `matches`
    - sequence: input sequence
    - stat: frequency of each category
    - matches: annotation with the highest count in each category.

For example

Input sequence
1) `sequence`: TTTTCTTTCACATTATAATGAAATAAG

5 kmer category 1, 2 kmer category 2, ect...

2)  `stat`:	{'1': 5, '2': 2, '3': 1, '4': 1, 'U': 1}

Out of 5 kmer in category 1, 3 kmer (higest count) belongs to GCA_004000535.1_ASM400053v1_genomic.fna,CP034509.1 Eukaryotic synthetic construct chromosome X.

3) `matches`:	{1: [(' GCA_004000535.1_ASM400053v1_genomic.fna,CP034509.1 Eukaryotic synthetic construct chromosome X', 3)], 2: [(' GCA_004000535.1_ASM400053v1_genomic.fna,3', 1)], 3: [(' final_purged_primary.fasta ctg.000095F,botznik-chr.fa chr5,GCA_004000535.1_ASM400053v1_genomic.fna CP034498.1 Eukaryotic synthetic construct chromosome 2,Carp_GCA_019924925.1_HZGC01_genomic.fna CM034330.1 Ctenopharyngodon idella isolate HZGC_01 chromosome 6,whole genome shotgun sequence', 1)]}


- [table_name].lookup_out_all.tsv: merged of `[table_name].lookup_out.tsv` and `[table_name].lookup_out_cleaned.tsv`

See `workflows/bin/lookup_table.py`


### Process: Blast
### What it does:
Executes BLAST on each split chunk using the following command:

```
blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseqid sgi sacc slen staxids stitle' -query {input_fasta} -remote -db nt -out {blast_out} -evalue 0.1 -task blastn -dust no -word_size 24 -reward 1 -penalty -3 -max_target_seqs 4

```
Definitions of output format can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK279684/#_appendices_Option_5_)

The outputs are stored in `blast_outs`. BLAST output columns includes:
```
["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",\
 "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "sgi", \
 "sacc", "slen", "staxids", "stitle"]
```

### Script:
```
workflows/bin/blast.py
```

## Process: BLAST Extract Features

### What it does:
- Extracts annotations from NCBI using the `Bio.Entrez` library.
- **Annotation Extraction:**
  - **Overlapped Region:** Region on the subject sequence that directly overlaps the BLAST alignment.
  - **5000 nt Window:** Annotations from features located within 5000 nucleotides upstream or downstream of the BLAST hit.
- **Conventions and Details:**
  - If multiple features (e.g., CDS, gene, regulatory elements) are found in the overlapped region, they are reported in a single field, separated by a delimiter (e.g., comma).
  - Features include information on whether they are in the sense or antisense orientation relative to the query.
  - The output is stored in `blast_features_outs` with four columns:
    - `query`: The input sequence.
    - `identity`: Percent identity to the subject sequence.
    - `features`: Annotations of the directly overlapped region.
    - `features_5000_window`: Annotations within a 5000 nucleotide window around the BLAST hit.


### Process: Summarize

### What it does:
- Merges the outputs from lookup tables (`*.lookup_out_all.tsv` files), BLAST results, and BLAST feature extractions.
- Retains the subject with the highest percent identity in cases of multiple BLAST hits.
- Outputs consolidated results as `annotation.tsv` in the `summary` folder.


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

In this example, input params is set to run on 200 anchors in `seqs.tsv` (this is an output from SPLASH) in `test_data`. 200 sequences will be spliited into 4 chunks for parallel processing.

2 look up tables is specified:
- `/oak/stanford/groups/horence/khoa/shares/lookup_tables/artifacts.slt` was built from `/oak/stanford/groups/horence/khoa/shares/lookup_tables/artifacts.txt`

### Lookup Table Build Process:
Lookup tables are built using the provided `build_lookup_table.py` script.

#### Example: Artifacts Lookup Table
```
build_lookup_table.py --poly_ACGT_len 6 --kmer_len 18 --bin_path /path/to/splash/bin --outname /path/to/lookup_tables/artifacts.slt /path/to/lookup_tables/artifacts.txt
```

#### Example: Common Microbe with Transcriptomes Lookup Table
```
build_lookup_table.py --poly_ACGT_len 6 --kmer_len 18 --bin_path /path/to/splash/bin --outname /path/to/lookup_tables/common_microbe_w_trnx.slt --transcriptomes /path/to/lookup_tables/microbes_transcriptomes.txt /path/to/lookup_tables/microbial_list.txt
```


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