/*  Run Lookup table on input folder containing fasta files
*/

process Lookup_table {
    tag "Lookup table $input_fasta"
    publishDir "${params.outdir}/lookup_out", mode: 'copy'

    input:
    path(input_fasta)
    
    output:
    path("*.lookup_out.tsv"), emit: lookup_out
    path("*.lookup_out_cleaned.tsv"), emit: lookup_out_cleaned

    script:
    
    """
    echo $params.lookup_table_paths
    python $params.bin/lookup_table.py --input $input_fasta --lookup_table_paths $params.lookup_table_paths --output "${input_fasta}.lookup_out.tsv"
    """
}

workflow {
    Lookup_table(input_fasta)
}

