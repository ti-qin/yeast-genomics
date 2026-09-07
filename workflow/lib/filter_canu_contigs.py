import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process folders and filter FASTA based on coverage using MAD method.")
    parser.add_argument("data_dir", type=str, help="Directory containing the subfolders to process.")
    return parser.parse_args()

# ✅ 最小改动：不直接过滤；返回 outlier mask + 阈值，便于在外面组合条件
def mad_outlier_mask(df, column, threshold=5):
    median = np.median(df[column])
    mad = np.median(np.abs(df[column] - median))

    upper = threshold * mad + median
    lower = -threshold * mad + median

    print(f"Median coverage: {median}")
    print(f"Upper threshold: {upper}")
    print(f"Lower threshold: {lower}")

    # MAD=0 时，不要把所有点都当离群；此时没有离群值
    if mad == 0:
        print("[WARN] MAD is 0, no coverage outliers will be removed.")
        return pd.Series(False, index=df.index), lower, upper

    is_outlier = (df[column] < lower) | (df[column] > upper)
    return is_outlier, lower, upper

def extract_fasta(fasta_file, formatted_tig_ids):
    output_fasta = []

    # 小优化（不改变逻辑）：用 set 加速 in 判断
    keep_set = set(formatted_tig_ids)

    for record in SeqIO.parse(fasta_file, "fasta"):
        tig_id = record.id
        if tig_id in keep_set:
            output_fasta.append(f">{tig_id}")
            output_fasta.append(str(record.seq))

    return output_fasta

def process_folder(folder, data_dir):
    sample = os.path.basename(folder.rstrip("/"))
    tiginfo_file = os.path.join(folder, f"{sample}.contigs.layout.tigInfo")
    fasta_file = os.path.join(folder, f"{sample}.contigs.fasta")

    if not os.path.exists(tiginfo_file) or not os.path.exists(fasta_file):
        print(f"Missing files in folder: {folder}")
        return

    tiginfo_df = pd.read_csv(tiginfo_file, sep='\t', comment='#', header=None)
    tiginfo_df = tiginfo_df.iloc[:, :8].copy()
    tiginfo_df.columns = ["tigID", "tigLen", "coverage", "tigClass", "sugRept", "sugBubb", "sugCirc", "numChildren"]

    print(f"Processing {folder}...")
    print(f"First 2 rows of tigInfo file:\n{tiginfo_df.head(2)}")

    # 基础筛选（不含 tigLen）
    base = tiginfo_df[
        (tiginfo_df['tigClass'] == 'contig') &
        (tiginfo_df['sugRept'] == 'no') &
        (tiginfo_df['sugBubb'] == 'no') &
        (tiginfo_df['sugCirc'] == 'no') &
        (tiginfo_df['coverage'] >= 18) &
        (tiginfo_df['tigLen'] >= 6000)
    ].copy()

    if base.empty:
        print("[WARN] No contigs after structural+coverage filter, skip.")
        return

    is_outlier, lower, upper = mad_outlier_mask(base, 'coverage', threshold=5)

    to_remove = is_outlier & (base['tigLen'] < 200000)

    print(f"Total contigs after base filter: {len(base)}")
    print(f"Coverage outliers: {int(is_outlier.sum())}")
    print(f"Removed (outlier & tigLen<200000): {int(to_remove.sum())}")
    print(f"Kept contigs: {len(base) - int(to_remove.sum())}")

    kept = base.loc[~to_remove].copy()

    # 生成保留的 tigID 并格式化成 tig00000001
    valid_tig_ids = kept['tigID'].astype(int).astype(str).tolist()
    formatted_tig_ids = [f"tig{int(tig_id):08d}" for tig_id in valid_tig_ids]

    filtered_fasta = extract_fasta(fasta_file, formatted_tig_ids)
    output_fasta_file = os.path.join(data_dir, f"{sample}.canu.fasta")
    with open(output_fasta_file, 'w') as f_out:
        for line in filtered_fasta:
            f_out.write(line + '\n')

    print(f"Filtered fasta saved for {folder} to {output_fasta_file}")

def main():
    args = parse_arguments()
    data_dir = args.data_dir

    for folder_name in os.listdir(data_dir):
        folder_path = os.path.join(data_dir, folder_name)
        if os.path.isdir(folder_path):
            process_folder(folder_path, data_dir)

if __name__ == "__main__":
    main()
