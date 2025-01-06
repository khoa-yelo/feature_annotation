/* Preprocess inputs for parallel processing in later steps*/

process Preprocess {
    tag "Preprocess $input"
    publishDir "${params.outdir}/preprocessed_inputs", mode: 'copy'

    input:
    path(input_file)
    val(input_type)
    val(annotate_column)
    
    output:
    path("input_chunks/*.fasta"), emit: input_chunks
    path("input.fasta"), emit: input_fasta

    script:
    """
    echo "cmd: python $params.bin/preprocess.py --input $params.input --input_type $params.input_type --annotate_column $params.annotate_column --output input_chunks"
    python $params.bin/preprocess.py --input $params.input --input_type $params.input_type --annotate_column $params.annotate_column --output input_chunks
    """
}

workflow {
    println params.input
    // check if annotate_column is provided when input_type is "tsv"
    if (params.input_type == "tsv" && params.annotate_column == null) {
        error "annotate_column must be provided when input_type is tsv"
    }
    Preprocess(params.input, params.input_type, params.annotate_column)
}

