import os
from os.path import join, basename
import subprocess
import shutil
import sys
import argparse

# SLURM_ARRAY_TASK_COUNT = int(os.environ['SLURM_ARRAY_TASK_COUNT'])
# SLURM_ARRAY_TASK_ID = int(os.environ['SLURM_ARRAY_TASK_ID'])

def run_blast(splitted_fasta, blast_out):

    fmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseqid sgi sacc slen staxids stitle"

    if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
        print(f"Skipping {splitted_fasta} as blast output already exists")
        sys.exit(0)
    cmd = f"blastn -outfmt '{fmt}' -query {splitted_fasta} -remote -db nt -out {blast_out} -evalue 0.1 -task blastn -dust no -word_size 24 -reward 1 -penalty -3 -max_target_seqs 4"
    subprocess.run(cmd, shell = True, check=True)
    print(f"Blast complete for {splitted_fasta}")

def main():
    parser = argparse.ArgumentParser(description='Run blast on the input data')
    parser.add_argument('--input', type=str, required=True, help='Input directory')
    args = parser.parse_args()

    print(args)
    print(f"Running blast on input data: {args.input}")
    output = args.input.rsplit(".", 1)[0] + ".blast_out.tsv"
    run_blast(args.input, output)

if __name__ == "__main__":
    main()