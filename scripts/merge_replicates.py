#!/usr/bin/env python3
"""
merge_replicates.py — Merge technical-replicate SRR runs into biological samples.

For the 18 multi-lane samples (66 SRR runs total), this script:
  1. Sums featureCounts read counts (gene_counts_clean.txt)
  2. Sums featureCounts assignment summary (gene_counts.txt.summary)
  3. Sums Kallisto estimated counts (count_matrix_kallisto.csv)
  4. Recalculates Kallisto TPM from merged counts (tpm_matrix_kallisto.csv)
  5. Consolidates variant_summary.tsv (sums variant counts per group)

Usage:
    python3 scripts/merge_replicates.py
"""

import os
import sys
import pandas as pd
import numpy as np

# ── Project paths ────────────────────────────────────────────────────────────
PROJECT_DIR = os.environ.get(
    "PROJECT_DIR",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)

# Input paths (original SRR-level data)
FC_COUNTS   = os.path.join(PROJECT_DIR, "counts", "gene_counts_clean.txt")
FC_SUMMARY  = os.path.join(PROJECT_DIR, "counts", "gene_counts.txt.summary")
KAL_COUNTS  = os.path.join(PROJECT_DIR, "kallisto_output", "count_matrix_kallisto.csv")
KAL_TPM     = os.path.join(PROJECT_DIR, "kallisto_output", "tpm_matrix_kallisto.csv")
VAR_SUMMARY = os.path.join(PROJECT_DIR, "variants", "variant_summary.tsv")

# Output paths (git-tracked merged data)
OUT_FC_COUNTS   = os.path.join(PROJECT_DIR, "counts_git", "gene_counts_clean.txt")
OUT_FC_SUMMARY  = os.path.join(PROJECT_DIR, "counts_git", "gene_counts.txt.summary")
OUT_KAL_COUNTS  = os.path.join(PROJECT_DIR, "kallisto_git", "count_matrix_kallisto.csv")
OUT_KAL_TPM     = os.path.join(PROJECT_DIR, "kallisto_git", "tpm_matrix_kallisto.csv")
OUT_VAR_SUMMARY = os.path.join(PROJECT_DIR, "variants_git", "variant_summary.tsv")

# ── Merge groups ─────────────────────────────────────────────────────────────
# Format: representative_SRR -> [all constituent SRRs]
# The representative keeps the name; others are summed in and dropped.
MERGE_GROUPS = {
    "SRR38105633": ["SRR38105633", "SRR38105634", "SRR38105635", "SRR38105636"],
    "SRR38105637": ["SRR38105637", "SRR38105638", "SRR38105640", "SRR38105641"],
    "SRR38105642": ["SRR38105642", "SRR38105643", "SRR38105644", "SRR38105645"],
    "SRR38105646": ["SRR38105646", "SRR38105647", "SRR38105648", "SRR38105649"],
    "SRR38105651": ["SRR38105651", "SRR38105652", "SRR38105695", "SRR38105696"],
    "SRR38105697": ["SRR38105697", "SRR38105698", "SRR38105699", "SRR38105700"],
    "SRR38105701": ["SRR38105701", "SRR38105702", "SRR38105704", "SRR38105705"],
    "SRR38105706": ["SRR38105706", "SRR38105707", "SRR38105708", "SRR38105709"],
    "SRR38105710": ["SRR38105710", "SRR38105711", "SRR38105712", "SRR38105713"],
    "SRR38105715": ["SRR38105715", "SRR38105716", "SRR38105818", "SRR38105819"],
    "SRR38105824": ["SRR38105824", "SRR38105825"],
    "SRR38105835": ["SRR38105835", "SRR38105836"],
    "SRR38105839": ["SRR38105839", "SRR38105840"],
    "SRR38105820": ["SRR38105820", "SRR38105821", "SRR38105822", "SRR38105823"],
    "SRR38105827": ["SRR38105827", "SRR38105828", "SRR38105829", "SRR38105830"],
    "SRR38105831": ["SRR38105831", "SRR38105832", "SRR38105833", "SRR38105834"],
    "SRR38105841": ["SRR38105841", "SRR38105891", "SRR38105892", "SRR38105893"],
    "SRR38105894": ["SRR38105894", "SRR38105895", "SRR38105896", "SRR38105897"],
}

# Build reverse lookup: SRR -> representative
SRR_TO_REP = {}
ALL_REPLICATES = set()
for rep, members in MERGE_GROUPS.items():
    for srr in members:
        SRR_TO_REP[srr] = rep
        if srr != rep:
            ALL_REPLICATES.add(srr)

print(f"Merge groups: {len(MERGE_GROUPS)}")
print(f"Total SRRs in groups: {sum(len(v) for v in MERGE_GROUPS.values())}")
print(f"SRRs to be merged away: {len(ALL_REPLICATES)}")


def merge_tsv_matrix(inpath, outpath, sep="\t"):
    """Merge a gene × sample matrix (TSV or CSV) by summing replicate columns."""
    if not os.path.exists(inpath):
        print(f"  ⚠ Skipping (not found): {inpath}")
        return None

    df = pd.read_csv(inpath, sep=sep, index_col=0)
    before_cols = df.shape[1]

    # Find which columns are in merge groups
    cols_present = set(df.columns)

    for rep, members in MERGE_GROUPS.items():
        present = [m for m in members if m in cols_present]
        if len(present) <= 1:
            continue  # No merging needed

        # Sum all member columns into the representative
        # Handle column names that may have path prefixes
        df[rep] = df[present].sum(axis=1)
        # Drop the non-representative columns
        to_drop = [m for m in present if m != rep]
        df.drop(columns=to_drop, inplace=True)

    after_cols = df.shape[1]
    merged = before_cols - after_cols

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    df.to_csv(outpath, sep=sep)

    print(f"  ✓ {os.path.basename(outpath)}: {before_cols} → {after_cols} samples (merged {merged})")
    return df


def merge_fc_matrix(inpath, outpath):
    """Merge featureCounts matrix — columns may have full path names."""
    if not os.path.exists(inpath):
        print(f"  ⚠ Skipping (not found): {inpath}")
        return None

    df = pd.read_csv(inpath, sep="\t", index_col=0)
    before_cols = df.shape[1]

    # featureCounts columns often have full path names — extract SRR IDs
    col_map = {}
    for col in df.columns:
        # Extract SRR ID from paths like /path/to/SRR38105552/SRR38105552
        for part in col.split("/"):
            if part.startswith("SRR"):
                col_map[col] = part
                break
        if col not in col_map:
            col_map[col] = col

    # Reverse map: SRR -> original column name
    srr_to_col = {v: k for k, v in col_map.items()}

    for rep, members in MERGE_GROUPS.items():
        present_cols = []
        for m in members:
            if m in srr_to_col:
                present_cols.append(srr_to_col[m])

        if len(present_cols) <= 1:
            continue

        rep_col = srr_to_col.get(rep)
        if rep_col is None:
            continue

        # Sum all member columns into the representative
        df[rep_col] = df[present_cols].sum(axis=1)
        to_drop = [c for c in present_cols if c != rep_col]
        df.drop(columns=to_drop, inplace=True)

    after_cols = df.shape[1]
    merged = before_cols - after_cols

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    df.to_csv(outpath, sep="\t")

    print(f"  ✓ {os.path.basename(outpath)}: {before_cols} → {after_cols} samples (merged {merged})")
    return df


def merge_variant_summary(inpath, outpath):
    """Merge variant summary TSV by summing counts for replicate groups."""
    if not os.path.exists(inpath):
        print(f"  ⚠ Skipping (not found): {inpath}")
        return None

    df = pd.read_csv(inpath, sep="\t")
    before_rows = len(df)

    rows_to_keep = []
    for _, row in df.iterrows():
        srr = row["Sample"]
        if srr in ALL_REPLICATES:
            continue  # Skip — will be summed into representative
        rows_to_keep.append(row)

    # Sum replicate group stats
    for rep, members in MERGE_GROUPS.items():
        member_rows = df[df["Sample"].isin(members)]
        if member_rows.empty:
            continue

        merged_row = {
            "Sample": rep,
            "Raw_Variants": member_rows["Raw_Variants"].sum(),
            "Filtered_Variants": member_rows["Filtered_Variants"].sum(),
            "SNPs": member_rows["SNPs"].sum(),
            "Indels": member_rows["Indels"].sum(),
        }
        rows_to_keep.append(pd.Series(merged_row))

    result = pd.DataFrame(rows_to_keep)
    # Remove duplicate representative rows (keep the merged one)
    result = result.drop_duplicates(subset="Sample", keep="last")
    result = result.sort_values("Sample").reset_index(drop=True)

    after_rows = len(result)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    result.to_csv(outpath, sep="\t", index=False)

    print(f"  ✓ {os.path.basename(outpath)}: {before_rows} → {after_rows} samples (merged {before_rows - after_rows})")
    return result


def recalculate_tpm(count_df, outpath):
    """Recalculate TPM from merged counts. For a single-gene-set (all same organism),
    TPM = (count / sum_of_counts) * 1e6 per sample."""
    if count_df is None:
        print(f"  ⚠ Skipping TPM recalculation (no count data)")
        return

    tpm_df = count_df.copy()
    for col in tpm_df.columns:
        total = tpm_df[col].sum()
        if total > 0:
            tpm_df[col] = (tpm_df[col] / total) * 1e6
        else:
            tpm_df[col] = 0.0

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    tpm_df.to_csv(outpath)

    print(f"  ✓ {os.path.basename(outpath)}: TPM recalculated for {tpm_df.shape[1]} samples")


def print_summary(fc_df, kal_df, var_df):
    """Print post-merge summary statistics."""
    print("\n" + "═" * 60)
    print("  Post-Merge Summary")
    print("═" * 60)

    if fc_df is not None:
        total_reads = fc_df.values.sum()
        n_samples = fc_df.shape[1]
        n_detected = (fc_df.sum(axis=0) > 0).sum()
        print(f"  featureCounts: {fc_df.shape[0]} genes × {n_samples} samples")
        print(f"    Total assigned reads: {total_reads:,.0f}")
        print(f"    Samples with detected virus: {n_detected}/{n_samples} ({100*n_detected/n_samples:.1f}%)")

    if kal_df is not None:
        n_samples = kal_df.shape[1]
        print(f"  Kallisto: {kal_df.shape[0]} genes × {n_samples} samples")

    if var_df is not None:
        n_with_var = (var_df["Filtered_Variants"] > 0).sum()
        total_snps = var_df["SNPs"].sum()
        total_indels = var_df["Indels"].sum()
        print(f"  Variants: {len(var_df)} samples")
        print(f"    With filtered variants: {n_with_var}/{len(var_df)} ({100*n_with_var/len(var_df):.1f}%)")
        print(f"    Total SNPs: {total_snps:,.0f}")
        print(f"    Total Indels: {total_indels:,.0f}")

    print("═" * 60)


def main():
    print("═" * 60)
    print("  Merging Technical Replicates")
    print("═" * 60)

    print("\n[1/5] featureCounts count matrix...")
    fc_df = merge_fc_matrix(FC_COUNTS, OUT_FC_COUNTS)

    print("\n[2/5] featureCounts summary...")
    merge_fc_matrix(FC_SUMMARY, OUT_FC_SUMMARY)

    print("\n[3/5] Kallisto count matrix...")
    kal_df = merge_tsv_matrix(KAL_COUNTS, OUT_KAL_COUNTS, sep=",")

    print("\n[4/5] Kallisto TPM matrix (recalculating)...")
    recalculate_tpm(kal_df, OUT_KAL_TPM)

    print("\n[5/5] Variant summary...")
    var_df = merge_variant_summary(VAR_SUMMARY, OUT_VAR_SUMMARY)

    print_summary(fc_df, kal_df, var_df)

    print("\nDone! Merged files written to *_git/ directories.")
    print("Next steps:")
    print("  1. python3 scripts/generate_plots.py")
    print("  2. python3 scripts/generate_plots_v2.py")
    print("  3. git add -f counts_git/ kallisto_git/ variants_git/ report/figures/")
    print("  4. git commit -m 'data: merge 18 technical replicate groups'")
    print("  5. git push origin mufakir/fix-pipeline-356-hybrid")


if __name__ == "__main__":
    main()
