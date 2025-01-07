"""
Module to summarize the results of the feature annotation workflow.
"""

import os
import sys
import pandas as pd
import argparse

    # df.columns = ["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",\
    #                  "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "sgi", \
    #                     "sacc", "slen", "staxids", "stitle"]

def parse_args():

    parser = argparse.ArgumentParser(description="Summarize the results of the feature annotation workflow.")
    parser.add_argument("--input_file", type=str, required = True, help="The input fasta or tsv file.")
    parser.add_argument("--input_type", type=str, required = True, help="The type of input file. Either 'fasta' or 'tsv'.")
    parser.add_argument("--annotate_column", type=str, required = True, help="The column in the input file that contains sequences to annotate.")
    parser.add_argument("--input_fasta", type=str, required = True, help="The input file converted to fasta format.")
    parser.add_argument("--lookup_outs", nargs="+", required = True, help="lookup table outputs.")
    parser.add_argument("--blast_outs", nargs="+", required = True, help="blast outputs")
    parser.add_argument("--blast_outs_features", nargs="+",required = True, help="blast feature outputs")
    return parser.parse_args()

def summarize_lookup(df_lookups, df_lookups_names, seq_col):
    for i, df in enumerate(df_lookups):
        #rename first column to seq_col
        df.rename(columns={df.columns[0]: seq_col}, inplace=True)
        #append name other columns with the name of the lookup table
        new_cols = []
        for col in df.columns:
            if col != seq_col:
                new_cols.append(col + "." + df_lookups_names[i])
                df[col + "." + df_lookups_names[i]] = df[col]
        #drop the original columns
        df.drop(columns=[col for col in df.columns if col not in [seq_col] + new_cols], inplace=True)
        # drop duplicated columns
        df = df.loc[:,~df.columns.duplicated()]
    df_lookup_summarized = pd.concat(df_lookups, axis=0)
    # drop duplicates
    df_lookup_summarized = df_lookup_summarized.drop_duplicates(subset=seq_col)
    return df_lookup_summarized

def summarize_blast(df_blasts, df_blast_features, seq_col):
    # concat all the blast outputs
    df_blast = pd.concat(df_blasts, axis=0)
    df_blast.columns = ["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",\
                            "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "sgi", \
                                "sacc", "slen", "staxids", "stitle"]
    df_blast_features = pd.concat(df_blast_features)
    # drop identity col in df_blast_features
    df_blast_features.drop(columns=["identity"], inplace=True)
    df_blast_summarized = pd.merge(df_blast, df_blast_features, on="query", how="left")
    # groupby query select row with max identity
    df_blast_summarized = df_blast_summarized.loc[df_blast_summarized.groupby("query")["identity"].idxmax()]
    # rename query to seq_col
    df_blast_summarized.rename(columns={"query": seq_col}, inplace=True)
    # group by seq_col and select row with max identity
    df_blast_summarized = df_blast_summarized.loc[df_blast_summarized.groupby(seq_col)["identity"].idxmax()]
    
    return df_blast_summarized

def summarize_blast_features(df_blast_features, df_input, seq_col):
    pass


def main():
    args = parse_args()
    print(args)

    df_summary = pd.DataFrame()
    df_blasts = [pd.read_csv(blast_out, sep="\t", header=None) for blast_out in args.blast_outs]
    df_features = [pd.read_csv(blast_out_features, sep="\t") for blast_out_features in args.blast_outs_features]
    df_lookups = [pd.read_csv(lookup_out, sep="\t") for lookup_out in args.lookup_outs]
    df_lookups_names = [os.path.basename(lookup_out).split(".")[0] for lookup_out in args.lookup_outs]

    seq_col = args.annotate_column if args.input_type == "tsv" else "sequence"
    if args.input_type == "tsv":
        df_input = pd.read_csv(args.input_file, sep="\t")
    elif args.input_type == "fasta":
        with open(args.input_file, "r") as f:
            lines = f.readlines()
        headers = [line for line in lines if line.startswith(">")]
        sequences = [line for line in lines if not line.startswith(">")]
        sequences = [seq.replace("\n", "") for seq in sequences]
        df_input = pd.DataFrame({"header": headers, "sequence": sequences})
        df_input["header"] = df_input["header"].str.replace(">", "")
    
    df_lookup_summarized = summarize_lookup(df_lookups, df_lookups_names, seq_col)
    df_blast_summarized = summarize_blast(df_blasts, df_features, seq_col)
    df_summary = pd.merge(df_lookup_summarized, df_blast_summarized, on=seq_col, how="outer")
    df_summary.to_csv("annotation.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()