/*  Run BLAST on input folder containing fasta files
*/

process Blast {
    tag "BLAST $input_chunk"
    publishDir "${params.outdir}/blast_outs", mode: 'copy'

    input:
    path(input_fasta)
    
    output:
    path("*.blast_out.tsv"), emit: blast_outs

    script:
    """
    echo "cmd: python $params.bin/blast.py --input $input_fasta --output "${input_fasta}.blast_out.tsv""
    python $params.bin/blast.py --input $input_fasta
    """
}

workflow {
    Blast(input_fasta)
}

