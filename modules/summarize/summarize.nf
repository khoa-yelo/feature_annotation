/* Summarize annotation outputs*/

process Summarize {
    tag "Summarize"
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path(input_file)
    val(input_type)
    val(annotate_column)
    path(input_fasta)
    path(lookup_outs)
    path(blast_outs)
    path(blast_outs_features)

    
    output:
    path("annotation.tsv"), emit: annotation

    script:
    """
    echo "cmd: python $params.bin/preprocess.py --input $params.input --input_type $params.input_type --annotate_column $params.annotate_column --output input_chunks"
    python $params.bin/preprocess.py --input $params.input --input_type $params.input_type --annotate_column $params.annotate_column --output input_chunks
    """
}