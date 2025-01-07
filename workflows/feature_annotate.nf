#!/usr/bin/env nextflow

/* Feature Annotation workflow */

include { Preprocess } from '../modules/preprocess/preprocess.nf'
include { Blast } from '../modules/blast/blast.nf'
include { Lookup_table } from '../modules/lookup_table/lookup_table.nf'
include { Blast_extract_features } from '../modules/blast_extract_features/blast_extract_features.nf'
include { Summarize } from '../modules/summarize/summarize.nf'

workflow {

    ch_processed_input = Preprocess(params.input, params.input_type, params.annotate_column)

    if (params.lookup_table) {
        ch_lookup_table_out = Lookup_table(ch_processed_input.input_fasta)
    }
    if (params.blast) {
        ch_blast_out = Blast(ch_processed_input.input_chunks.flatten())
        ch_blast_feat_out = Blast_extract_features(ch_blast_out.blast_outs.flatten())
    }
    Summarize(params.input, params.input_type, params.annotate_column, ch_processed_input.input_fasta, ch_lookup_table_out.lookup_out_all.collect(), ch_blast_out.blast_outs.collect(), ch_blast_feat_out.blast_outs_features.collect())
}
