#!/usr/bin/env python3
"""
check_origins.py
Annotate a merged TSV with all possible database origins of each Peak Data sequence.

Input:
  - merged TSV with a 'Peak Data' column
  - FASTA database used by FragmentFinder

Output:
  - <merged_basename>_checked.tsv with:
      Alt_Origins_Count
      Possible Origin 1 ... Possible Origin N  (each as "ID:start-end")

Coordinates are 0-based, half-open [start:end).
"""

import os
import sys
import argparse
import pandas as pd

def iter_fasta(path):
    """Yield (record_id, sequence) from a FASTA file (pure Python, no Biopython)."""
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        yield header, "".join(seq_chunks)

def find_all_occurrences(haystack: str, needle: str):
    """Return list of (start, end) 0-based half-open coords for all (including overlapping) matches."""
    if not needle:
        return []
    i = 0
    n = len(needle)
    hits = []
    while True:
        j = haystack.find(needle, i)
        if j < 0:
            break
        hits.append((j, j + n))
        i = j + 1
    return hits

def build_peak_to_origins_map(db_fasta_path: str, peak_strings):
    result = {p: [] for p in peak_strings}
    if not peak_strings:
        return result
    for rec_id, seq in iter_fasta(db_fasta_path):
        for p in peak_strings:
            hits = find_all_occurrences(seq, p)
            if hits:
                result[p].extend((rec_id, s, e) for (s, e) in hits)
    return result

def main(merged_tsv: str, db_fasta_path: str, out_path: str | None = None):
    if not os.path.isfile(merged_tsv):
        sys.exit(f"ERROR: merged TSV not found: {merged_tsv}")
    if not os.path.isfile(db_fasta_path):
        sys.exit(f"ERROR: database FASTA not found: {db_fasta_path}")

    df = pd.read_csv(merged_tsv, sep="\t")
    if "Peak Data" not in df.columns:
        sys.exit("ERROR: input TSV missing required column 'Peak Data'.")

    unique_peaks = df["Peak Data"].dropna().astype(str).unique().tolist()
    peak2origins = build_peak_to_origins_map(db_fasta_path, unique_peaks)

    origins_counts = []
    origins_lists = []
    for _, row in df.iterrows():
        pk = "" if pd.isna(row.get("Peak Data")) else str(row.get("Peak Data"))
        origins = peak2origins.get(pk, [])
        origins_counts.append(len(origins))
        origins_lists.append(origins)

    df["Alt_Origins_Count"] = origins_counts

    max_k = max(origins_counts) if origins_counts else 0
    for k in range(1, max_k + 1):
        col = f"Possible Origin {k}"
        out_col = []
        for origins in origins_lists:
            if len(origins) >= k:
                rid, s, e = origins[k - 1]
                out_col.append(f"{rid}:{s}-{e}")
            else:
                out_col.append("")
        df[col] = out_col

    if out_path is None:
        base, ext = os.path.splitext(merged_tsv)
        out_path = f"{base}_checked.tsv"

    df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--merged", required=True)
    ap.add_argument("-d", "--db", required=True)
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()
    main(a.merged, a.db, a.out)
