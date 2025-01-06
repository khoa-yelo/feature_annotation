/*  Extract features from BLAST subject sequence
*/

process Blast_extract_features {
    tag "BLAST feature $blast_outs"
    publishDir "${params.outdir}/blast_features_outs", mode: 'copy'

    input:
    path(blast_out)
    
    output:
    path("*.features.tsv"), emit: blast_outs_features

    script:
    """
    echo "cmd: python $params.bin/blast_extract_features.py --input $blast_out"
    python $params.bin/blast_extract_features.py --blast_out $blast_out
    """
}

workflow {
    Blast_extract_features(input_fasta)
}

