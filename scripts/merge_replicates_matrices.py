#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

PROJECT_DIR = "/home/ansari/Desktop/Assignments/ebola_class_project"
COUNTS_DIR = os.path.join(PROJECT_DIR, "counts_git")
KALLISTO_DIR = os.path.join(PROJECT_DIR, "kallisto_git")
VARIANTS_DIR = os.path.join(PROJECT_DIR, "variants_git")

merge_groups = {
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
    "SRR38105894": ["SRR38105894", "SRR38105895", "SRR38105896", "SRR38105897"]
}

def merge_dataframe(df, is_tpm=False):
    # Some columns might not be SRRs (like Geneid or Status).
    # We will process SRR columns only.
    cols_to_keep = []
    
    # Create a mapping from target to its existing sources in the df
    srr_cols = [c for c in df.columns if "SRR" in c]
    other_cols = [c for c in df.columns if "SRR" not in c]
    
    # We want to rename things properly and sum/mean
    new_data = {}
    for col in other_cols:
        new_data[col] = df[col]
        
    processed_srrs = set()
    
    for target, sources in merge_groups.items():
        # Find which sources actually exist in the dataframe
        valid_sources = [s for s in sources if s in srr_cols or any(s in c for c in srr_cols)]
        if not valid_sources:
            continue
            
        # Get the exact column names in df that correspond to valid_sources
        exact_cols = []
        for vs in valid_sources:
            for c in srr_cols:
                if vs in c:
                    exact_cols.append(c)
                    processed_srrs.add(c)
        
        if exact_cols:
            if is_tpm:
                # Average TPM
                new_data[target] = df[exact_cols].mean(axis=1)
            else:
                # Sum counts
                new_data[target] = df[exact_cols].sum(axis=1)
                
    # Add back all the SRRs that were not part of any merge group
    for c in srr_cols:
        if c not in processed_srrs:
            # We must rename columns that have long paths in featureCounts
            # e.g. /users/PWSU0516/mufakiransari/ebola_class_project/aligned/SRR38105552/SRR38105552.dedup.bam -> SRR38105552
            clean_name = c
            if ".bam" in c:
                clean_name = os.path.basename(c).split(".")[0]
            new_data[clean_name] = df[c]
            
    # Reconstruct dataframe
    new_df = pd.DataFrame(new_data)
    return new_df

def process_featurecounts():
    fc_file = os.path.join(COUNTS_DIR, "gene_counts_clean.txt")
    if os.path.exists(fc_file):
        df = pd.read_csv(fc_file, sep="\t")
        merged_df = merge_dataframe(df, is_tpm=False)
        merged_df.to_csv(fc_file, sep="\t", index=False)
        print(f"featureCounts: {df.shape[1]-1} -> {merged_df.shape[1]-1} samples")

    summ_file = os.path.join(COUNTS_DIR, "gene_counts.txt.summary")
    if os.path.exists(summ_file):
        df = pd.read_csv(summ_file, sep="\t")
        merged_df = merge_dataframe(df, is_tpm=False)
        merged_df.to_csv(summ_file, sep="\t", index=False)

def process_kallisto():
    counts_file = os.path.join(KALLISTO_DIR, "count_matrix_kallisto.csv")
    if os.path.exists(counts_file):
        df = pd.read_csv(counts_file)
        # Rename Unnamed: 0 to gene
        if df.columns[0] == "Unnamed: 0":
            df = df.rename(columns={"Unnamed: 0": "gene"})
        merged_df = merge_dataframe(df, is_tpm=False)
        merged_df.to_csv(counts_file, index=False)
        print(f"Kallisto counts: {df.shape[1]-1} -> {merged_df.shape[1]-1} samples")

    tpm_file = os.path.join(KALLISTO_DIR, "tpm_matrix_kallisto.csv")
    if os.path.exists(tpm_file):
        df = pd.read_csv(tpm_file)
        if df.columns[0] == "Unnamed: 0":
            df = df.rename(columns={"Unnamed: 0": "gene"})
        merged_df = merge_dataframe(df, is_tpm=True)
        merged_df.to_csv(tpm_file, index=False)

def process_variants():
    var_file = os.path.join(VARIANTS_DIR, "variant_summary.tsv")
    if os.path.exists(var_file):
        df = pd.read_csv(var_file, sep="\t")
        
        # We will keep the row with the max Filtered_Variants for each group
        to_drop = set()
        for target, sources in merge_groups.items():
            valid_sources = [s for s in sources if s in df["Sample"].values]
            if not valid_sources:
                continue
            
            # Find the row with max Filtered_Variants
            subset = df[df["Sample"].isin(valid_sources)]
            max_idx = subset["Filtered_Variants"].idxmax()
            
            # Keep max_idx, drop others
            for idx in subset.index:
                if idx != max_idx:
                    to_drop.add(idx)
            
            # Rename the kept one to the target name
            df.at[max_idx, "Sample"] = target
            
        merged_df = df.drop(index=list(to_drop))
        merged_df.to_csv(var_file, sep="\t", index=False)
        print(f"Variants: {df.shape[0]} -> {merged_df.shape[0]} samples")

if __name__ == "__main__":
    process_featurecounts()
    process_kallisto()
    process_variants()
    print("Done merging matrices.")
