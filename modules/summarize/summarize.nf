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
    echo "CMD: python $params.bin/summarize.py --input_file $input_file --input_type $input_type --annotate_column $annotate_column --input_fasta $input_fasta --lookup_outs ${lookup_outs.collect { "'${it}'" }.join(' ')} --blast_outs ${blast_outs.collect { "'${it}'" }.join(' ')} --blast_outs_features ${blast_outs_features.collect { "'${it}'" }.join(' ')}"
    python $params.bin/summarize.py --input_file $input_file --input_type $input_type --annotate_column $annotate_column --input_fasta $input_fasta --lookup_outs ${lookup_outs.collect { "'${it}'" }.join(' ')} --blast_outs ${blast_outs.collect { "'${it}'" }.join(' ')} --blast_outs_features ${blast_outs_features.collect { "'${it}'" }.join(' ')}
    """
}