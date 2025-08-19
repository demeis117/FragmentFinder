import numpy as np
from scipy.signal import find_peaks
import csv
import os
from typing import Dict, List, Union
from random import choices

def _val_check(vec: np.ndarray, idx: int, threshold: float):
    left = idx
    right = idx
    for i in range(idx - 1, -1, -1):
        if vec[i] > threshold:
            left = i
        else:
            break
    for i in range(idx + 1, len(vec)):
        if vec[i] > threshold:
            right = i
        else:
            break
    return left, right

def get_peaks(
    results_dict: Dict[str, np.ndarray],
    original_record_dict: Dict[str, str],
    fasta_file_name: str,
    db_file_name: str,
    total_reads: Union[int, float],
    fasta_path: str = "",
):
    rows = []
    info_name = f"{fasta_file_name}({db_file_name})"

    if total_reads == 0:
        raise ValueError("Total reads is zero. Cannot compute RPM.")

    global_background = []

    # First pass: identify peaks and build isolated per-peak background
    peak_metadata = []
    for uid, vec in results_dict.items():
        vec = np.asarray(vec, dtype=np.float64)
        sd = np.std(vec)
        peaks, _ = find_peaks(vec, prominence=sd * 0.85, distance=18)
        seq = original_record_dict[uid]

        for idx in peaks:
            left, right = _val_check(vec, idx, vec[idx] - sd * 0.1)
            peak_start = max(0, left)
            peak_end = min(right + 18, len(vec))
            rpm = round(float(vec[idx] / total_reads * 1_000_000), 4)
            if rpm < 5.0:
                continue

            # Create an isolated exclusion window around current peak
            exclude_mask = np.zeros(len(vec), dtype=bool)
            start = max(0, left - 5)
            end = min(len(vec), right + 5 + 18)
            exclude_mask[start:end] = True
            local_background = vec[~exclude_mask]

            if len(local_background) >= 2:
                global_background.extend(local_background)

            peak_seq = seq[peak_start:peak_end]
            peak_end = peak_start + len(peak_seq)
            peak_metadata.append((uid, idx, vec[idx], rpm, peak_start, peak_end, peak_seq))

    global_background = np.array(global_background)
    num_peaks = len(peak_metadata)

    global_mean = np.mean(global_background)
    global_std = np.std(global_background)

    # Dynamic permutation count
    if num_peaks < 10000:
        n_perm = 1000000
    elif num_peaks > 10000:
        n_perm = 200000

    for uid, idx, peak_val, rpm, peak_start, peak_end, peak_seq in peak_metadata:
        # Permutation test
        permutations = choices(global_background, k=n_perm)
        count_extreme = sum([1 for val in permutations if val >= peak_val])
        p_val = (1 + count_extreme) / (1 + n_perm)

        rows.append([
            uid,
            rpm,
            peak_start,
            peak_end,
            peak_seq,
            f"{p_val:.50e}",
        round((peak_val - global_mean) / global_std, 4),
            int(total_reads),
            info_name
        ])

    if not rows:
        print("[Util] No qualifying peaks found. No CSV written.")
        return []

    outfile = f"{fasta_file_name}({db_file_name}).tsv"
    with open(outfile, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow([
        "miRNA ID",
        "Reads per Million",
        "Peak Start",
        "Peak End",
        "Peak Data",
        "Permutation P-Value",
        "Z-score Rank",
            "Total Reads in NGS File",
            "File Name"
        ])
        writer.writerows(rows)

    print(f"[Util] CSV written → {os.path.abspath(outfile)}")
    return rows