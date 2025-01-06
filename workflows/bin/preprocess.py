"""
Module to preprocess input for parallel processing in later steps
Khoa Hoang
01/06/25
"""
import os
from os.path import join
import argparse
import pandas as pd
from Bio import SeqIO
import shutil
import pandas as pd

SPLIT_THRESH = 100
SPLIT_EACH = 50 

def read_fasta(fasta_file, output_type="dict"):
    """
    Read a fasta file and return a dictionary with the sequence id as key and the sequence as value.
    """
    if output_type == "list":
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
        return sequences
    elif output_type == "dict":
        sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = record.seq
        return sequences
    elif output_type == "pandas":
        sequences = []
        description = []
        ids = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
            description.append(record.description)
            ids.append(record.id)
        return pd.DataFrame({"ID": ids, "Description": description, "Sequence": sequences})


def split_fasta(fasta_file, output_dir, num_seq=1):
    """
    Split a fasta file into multiple files.
    """
    os.makedirs(output_dir, exist_ok=True)
    if num_seq == 1:
        for record in SeqIO.parse(fasta_file, "fasta"):
            output_file = os.path.join(output_dir, record.id + ".fasta")
            with open(output_file, "w") as f:
                f.write(">" + record.description + "\n")
                f.write(str(record.seq) + "\n")
    else:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        num_files = len(records) // num_seq
        for i in range(num_files):
            output_file = os.path.join(output_dir, f"split_{i}.fasta")
            with open(output_file, "w") as f:
                for record in records[i*num_seq:(i+1)*num_seq]:
                    f.write(">" + record.description + "\n")
                    f.write(str(record.seq) + "\n")

def write_fasta(df, output_file, id_col, description_col, sequence_col):
    """
    Write a pandas dataframe to a fasta file.
    """
    with open(output_file, "w") as f:
        for index, row in df.iterrows():
            f.write(">" + str(row[id_col]) + " " + str(row[description_col]) + "\n")
            f.write(str(row[sequence_col]) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Preprocess the input data')
    parser.add_argument('--input', type=str, required=True, help='Input file')
    parser.add_argument('--input_type', type=str, required=True, help='Input type')
    parser.add_argument('--annotate_column', type=str, required=True, help='column to annotate')
    parser.add_argument('--output', type=str, required=True, help='Output directory')
    args = parser.parse_args()

    print(args)
    os.makedirs(args.output, exist_ok=True)
    print(f"Preprocessing input data: {args.input}")

    fasta_file = args.input
    if args.input_type == 'tsv':
        fasta_file = "temp.fasta"
        df = pd.read_csv(args.input, sep='\t')
        seqs = df[args.annotate_column].tolist()
        write_fasta(df, fasta_file, args.annotate_column, args.annotate_column, args.annotate_column)
    elif args.input_type == 'fasta':
        seqs = read_fasta(fasta_file)
    shutil.copy(fasta_file, "input.fasta")
    print(f"Number of sequences: {len(seqs)}")
    print("Fasta files: ", fasta_file)
    # Split the sequences
    if len(seqs) > SPLIT_THRESH:
        split_fasta(fasta_file, args.output, SPLIT_EACH)
    else:
        shutil.copy(fasta_file, args.output)
        sys.exit(0)

if __name__ == '__main__':
    main()