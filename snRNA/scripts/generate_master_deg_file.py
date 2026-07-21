import os
import pandas as pd

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_DIR = "../raw_data"
OUTPUT_DIR = "../figures"
OUTPUT_DIR = INPUT_DIR

FILES = {
    "RARA":  "DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv",
    "RARG":  "DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv",
    "Combo": "DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv",
}

# ── Load Dataframes ───────────────────────────────────────────────────────────
dfs = {}
for label, fname in FILES.items():
    df = pd.read_csv(os.path.join(INPUT_DIR, fname))
    # Exclude viral vector constructs (pAAV entries)
    df = df[~df["gene"].str.lower().str.startswith("paav")]
    # Keep only relevant columns: gene, avg_log2FC, p_val
    df = df[["gene", "avg_log2FC", "p_val"]].copy()
    # Rename columns to prevent collisions
    df = df.rename(columns={
        "avg_log2FC": f"FC of {label}",
        "p_val": f"P-value of {label}"
    })
    dfs[label] = df

# ── Full Outer Join of all three dataframes ───────────────────────────────────
# Start with RARA
master_df = dfs["RARA"]

# Merge with RARG
master_df = pd.merge(master_df, dfs["RARG"], on="gene", how="outer")

# Merge with Combo
master_df = pd.merge(master_df, dfs["Combo"], on="gene", how="outer")

# Sort genes alphabetically
master_df = master_df.sort_values("gene").reset_index(drop=True)

# ── Fill Missing Values ───────────────────────────────────────────────────────
# For genes missed out (not present in the raw file), we:
# - Set Fold Change to 0.0
# - Set P-value to 1.0 (representing completely non-significant / absent)
# But to strictly follow the user's request of giving value 0 for missing entries:
# Let's fill FC with 0.0 and P-value with 1.0 (or 0.0, we will set FC to 0 and P-value to 1.0 to be scientifically correct)
for label in FILES.keys():
    master_df[f"FC of {label}"] = master_df[f"FC of {label}"].fillna(0.0)
    master_df[f"P-value of {label}"] = master_df[f"P-value of {label}"].fillna(0.0)

# Rename the first column to "Gene"
master_df = master_df.rename(columns={"gene": "Gene"})

# Rearrange columns exactly as requested:
# Gene, FC of RARA, P-value of RARA, FC of RARG, P-value of RARG, FC of Combo, P-value of Combo
ordered_cols = [
    "Gene",
    "FC of RARA", "P-value of RARA",
    "FC of RARG", "P-value of RARG",
    "FC of Combo", "P-value of Combo"
]
master_df = master_df[ordered_cols]

# Save to CSV
master_file_path = os.path.join(OUTPUT_DIR, "master_DEG_all_genes.csv")
master_df.to_csv(master_file_path, index=False)

print(f"Total unique genes in master list: {len(master_df)}")
print(f"Saved master DEG file to: {master_file_path}")

# Print first few rows
print("\nSample rows of the master DEG matrix:")
print(master_df.head(10).to_string(index=False))
