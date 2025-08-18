import pandas as pd
from Bio import SeqIO
import glob
import os
import ast
import re
import sys

# Import the origin checker
from check_origins import main as check_origins_main

EXCLUDE_BASENAMES = {"t1.tsv", "grouped_by_peak.tsv", "files_and_reads.tsv"}

def sanitize_dir(p: str) -> str:
    p = (p or "").strip().strip('"').strip("'")
    return os.path.normpath(p)

def should_exclude(path: str) -> bool:
    base = os.path.basename(path)
    # Exclude intermediates and outputs
    if base in EXCLUDE_BASENAMES:
        return True
    if base.endswith("_checked.tsv"):
        return True
    # Anything under Joined_Results should not be re-merged
    if "Joined_Results" in path.replace("\\", "/").split("/"):
        return True
    return False

def get_files(dir_name):
    dir_name = sanitize_dir(dir_name)
    patterns = [os.path.join(dir_name, "*.tsv"), os.path.join(dir_name, "*.TSV")]
    joined_list = []
    for pat in patterns:
        joined_list.extend(glob.glob(pat))
    # Filter
    joined_list = [p for p in joined_list if not should_exclude(p)]
    return sorted(set(joined_list))

def show_files(lst):
    for file in lst:
        print(file)

def get_list_count(lst):
    return len(lst)

def get_new_out_file_name(lst):
    if not lst:
        raise ValueError("No input TSV files found; cannot derive output file name.")
    s = lst[0]
    base = os.path.basename(s)
    paren = re.findall(r'\(.*?\)', base)
    if paren:
        return paren[0]
    return os.path.splitext(base)[0]

def read_many_tsv(files):
    frames = []
    for f in files:
        try:
            if os.path.getsize(f) == 0:
                print(f"[warn] Skipping zero-byte file: {f}")
                continue
            df = pd.read_csv(f, sep='\t')
            # basic sanity: must include core columns
            required = {"miRNA ID", "Peak Start", "Peak End", "Reads per Million", "File Name"}
            if not required.issubset(df.columns):
                print(f"[warn] Skipping {f} (missing required columns: {required - set(df.columns)})")
                continue
            frames.append(df)
        except pd.errors.EmptyDataError:
            print(f"[warn] Skipping empty/invalid TSV: {f}")
        except Exception as e:
            print(f"[warn] Skipping {f}: {e}")
    if not frames:
        raise RuntimeError("After filtering, no valid TSVs remained to merge.")
    return pd.concat(frames, ignore_index=True)

def merge_files(dir_name):
    dir_name = sanitize_dir(dir_name)
    joined_list = get_files(dir_name)
    if not joined_list:
        raise FileNotFoundError(f"No eligible .tsv files found in: {dir_name}")
    df = read_many_tsv(joined_list)
    df.sort_values(['miRNA ID', 'Peak Start'], ascending=True, inplace=True)
    df.to_csv(f'{dir_name}/t1.tsv', sep='\t', index=False)

# Define a function to check if two values are within a given range
def within_range(val1, val2, range_value):
    if val1 is None or val2 is None:
        return False
    return abs(val1 - val2) <= range_value

def init_group(dir_name):
    dir_name = sanitize_dir(dir_name)
    tsv_file = f'{dir_name}/t1.tsv'
    if not os.path.isfile(tsv_file):
        raise FileNotFoundError(f"Expected intermediate file not found: {tsv_file}")
    df = pd.read_csv(tsv_file, sep='\t')

    # Sort the DataFrame by 'Peak Start'
    df.sort_values(by=['miRNA ID', 'Peak Start'], inplace=True)

    # Initialize variables
    grouped_data = []
    current_group = []
    previous_start = None
    previous_ID = None
    current_group_start = None

    # Iterate through the DataFrame rows
    for _, row in df.iterrows():
        peak_start = row['Peak Start']
        peak_end = row['Peak End']
        id_val = row['miRNA ID']

        if (previous_ID is None) or (id_val != previous_ID) or (current_group_start is None) or (not within_range(peak_start, previous_start, 10)) or (peak_end - current_group_start > 40):
            if current_group:
                grouped_data.append(current_group)
            current_group = [row]
            current_group_start = peak_start
        else:
            current_group.append(row)
            
        previous_start = peak_start
        previous_ID = id_val

    # Append the last group
    if current_group:
        grouped_data.append(current_group)

    output_tsv = f'{dir_name}/grouped_by_peak.tsv'
    with open(output_tsv, 'w', encoding='utf-8') as f:
        f.write('Group\tmiRNA ID\tPeak Start\tPeak End\tReads per Million\tFile Name\n')
        for idx, group in enumerate(grouped_data, start=1):
            for item in group:
                f.write(f"Group {idx}\t{item['miRNA ID']}\t{item['Peak Start']}\t{item['Peak End']}\t{item['Reads per Million']}\t{item['File Name']}\n")

def get_orig_strings(db):
    db = sanitize_dir(db)
    if not os.path.isfile(db):
        raise FileNotFoundError(f"Database FASTA not found: {db}")
    records = {}
    with open(db, 'r', encoding='utf-8') as f:
        for record in SeqIO.parse(f, "fasta"):
            s = str(record.seq)
            if s:
                records[record.id] = s
    if not records:
        raise ValueError("Parsed 0 records from the database FASTA.")
    return records

def final_group(group, arg2):
    min_start = group['Peak Start'].min()
    max_end = group['Peak End'].max()

    seq_id = group['miRNA ID'].iloc[0]
    if seq_id not in arg2:
        raise KeyError(f"Sequence ID '{seq_id}' not found in database FASTA.")
    seq = arg2[seq_id]
    if max_end > len(seq):
        max_end = len(seq)
    peak_data = seq[int(min_start):int(max_end)]
    
    reads_and_files = []
    for _, row in group.iterrows():
        reads_and_files.append({row['File Name']: row['Reads per Million']})
    
    return pd.DataFrame({
            'miRNA ID': [seq_id],
            'Peak Start': [min_start],
            'Peak End': [max_end],
            'Files : Reads': [reads_and_files],
            'Peak Data': [peak_data]
        })

def get_final_group(dir_name, hmap):
    dir_name = sanitize_dir(dir_name)
    file = f'{dir_name}/grouped_by_peak.tsv'
    if not os.path.isfile(file):
        raise FileNotFoundError(f"Expected intermediate file not found: {file}")
    df2 = pd.read_csv(file, sep='\t')
    grouped = df2.groupby(['Group', 'miRNA ID']).apply(final_group, arg2=hmap)
    grouped.rename_axis(index={'miRNA ID': 'ID'}, inplace=True)
    grouped = grouped.reset_index(drop=True)
    grouped = grouped.sort_values(by=['miRNA ID', 'Peak Start'])
    grouped.to_csv(f'{dir_name}/files_and_reads.tsv', sep='\t', index=False)

def export_final_tsv(dir_name, count, fname):
    dir_name = sanitize_dir(dir_name)
    tsv_file = f'{dir_name}/files_and_reads.tsv'
    if not os.path.isfile(tsv_file):
        raise FileNotFoundError(f"Expected intermediate file not found: {tsv_file}")
    original = pd.read_csv(tsv_file, sep='\t')

    original['Files : Reads'] = original['Files : Reads'].apply(ast.literal_eval)

    all_keys = set()
    for _, row in original.iterrows():
        reads_list = row['Files : Reads']
        for dictionary in reads_list:
            all_keys.update(dictionary.keys())

    data_dict = {key: [] for key in all_keys}

    for _, row in original.iterrows():
        reads_list = row['Files : Reads']
        row_data = {key: 0 for key in all_keys}
        for dictionary in reads_list:
            for key, value in dictionary.items():
                if float(row_data[key]) < float(value):
                    row_data[key] = float(value)
        for key, value in row_data.items():
            data_dict[key].append(value)

    read_df = pd.DataFrame(data_dict)
    extracted_col = original['Peak Data']
    original = original.drop(columns=['Peak Data'])
    final_df = pd.concat([original, read_df], axis=1)
    final_df = final_df.drop(columns=['Files : Reads'])
    final_df.insert(count+3, 'Peak Data', extracted_col)
    final_df = final_df.sort_values(by=['miRNA ID', 'Peak Start'])
    os.makedirs(f'{dir_name}/Joined_Results', exist_ok=True)
    final_df.to_csv(f'{dir_name}/Joined_Results/{fname}.tsv', sep='\t', index=False)

def main():
    directory_name = input("Enter Full Path to Directory Containing Files to Merge: ")
    directory_name = sanitize_dir(directory_name)
    db_file_path = input("Enter Full Path to Database File Used for Sequence Alignment: ")
    db_file_path = sanitize_dir(db_file_path)

    file_list = get_files(directory_name)
    if not file_list:
        print(f"No eligible .tsv files found in: {directory_name}")
        sys.exit(1)

    show_files(file_list)

    print()
    print('Merging and Grouping Files...')
    print()

    out_file = get_new_out_file_name(file_list)
    file_count = get_list_count(file_list)

    merge_files(directory_name)
    init_group(directory_name)
    hmap = get_orig_strings(db_file_path)
    get_final_group(directory_name, hmap)
    export_final_tsv(directory_name, file_count, out_file)

    # --- Cross-check Peak Data origins against the database ---
    merged_path = f"{directory_name}/Joined_Results/{out_file}.tsv"
    db_path = db_file_path  # same FASTA used above
    checked_out = os.path.splitext(merged_path)[0] + "_checked.tsv"

    check_origins_main(merged_path, db_path, checked_out)
    print(f"[merge] Origin cross-check complete -> {checked_out}")

    print(f"Files Merged and Output Exported to 'Joined_Results/{out_file}.tsv'")
    try:
        os.remove(f'{directory_name}/files_and_reads.tsv')
        os.remove(f'{directory_name}/grouped_by_peak.tsv')
        os.remove(f'{directory_name}/t1.tsv')
    except OSError:
        pass

if __name__ == "__main__":
    main()
