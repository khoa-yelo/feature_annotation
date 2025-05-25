from Bio import Entrez, SeqIO
import pandas as pd
import time
import sys
import os
from os.path import join, basename
import argparse

Entrez.email = "khoang99@stanford.edu"
MAX_RETRIES = 100
BLAST_WINDOW = 5000

def fetch_sequence(seq_id):
    print(f"Fetching sequence {seq_id}")
    # retry until success, sleep for 5 seconds between each retry, max retries = MAX_RETRIES
    for i in range(MAX_RETRIES):
        try:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            break
        except:
            time.sleep(5)
    return record
 
def find_overlapping_features(record, window_start, window_end):
    overlapping_features = []
    for feature in record.features:
        if feature.type == "source":
            continue
        feature_start = feature.location.start
        feature_end = feature.location.end

        # Check if the feature overlaps with the specified window
        overlap = (feature_start <= window_end) and (feature_end >= window_start)

        if overlap:
            overlapping_features.append({
                "type": feature.type,
                "start": str(feature_start),
                "end": str(feature_end),
                "gene": feature.qualifiers.get("gene"),
                "product": feature.qualifiers.get("product"),
                "protein_seq": feature.qualifiers.get("translation")
            })

    return overlapping_features

def featurize_blast_out(blast_out, window=5000):
    # Check for empty file
    if os.path.getsize(blast_out) == 0:
        return pd.DataFrame(columns=["query", "identity", "features", f"features_{window}_window"])

    df = pd.read_csv(blast_out, sep="\t", header=None)
    df.columns = [
        "query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",
        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "sgi",
        "sacc", "slen", "staxids", "stitle"
    ]

    df["features"] = None
    df[f"features_{window}_window"] = None
    sacc_records = {}

    for sacc in df["sacc"].unique():
        sacc_records[sacc] = fetch_sequence(sacc)

    for index, row in df.iterrows():
        record = sacc_records[row["sacc"]]
        features = find_overlapping_features(record, row["s_start"], row["s_end"])
        df.at[index, "features"] = features

        # get the window
        window_start = max(row["s_start"] - window, 0)
        window_end = row["s_end"] + window
        features = find_overlapping_features(record, window_start, window_end)
        df.at[index, f"features_{window}_window"] = features

    return df[["query", "identity", "features", f"features_{window}_window"]]

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Featurize blast output')
    parser.add_argument('--blast_out', type=str, required=True, help='Blast out file')
    args = parser.parse_args()
    blast_feat_out = args.blast_out.rsplit(".", 1)[0] + ".features.tsv"   
    df_features = featurize_blast_out(args.blast_out, window=BLAST_WINDOW)
    df_features.to_csv(blast_feat_out, index = None, sep = "\t")
    print(f"Featurize blast output complete for {args.blast_out}. Output file: {blast_feat_out}")
