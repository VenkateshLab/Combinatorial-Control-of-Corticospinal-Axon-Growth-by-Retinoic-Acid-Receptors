import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_DIR = "../raw_data"
OUTPUT_DIR = "../figures"
DAVID_DIR = INPUT_DIR

FILES = {
    "Up":   "DAVIDChartReport_All_3_Up_2026-07-08.csv",
    "Down": "DAVIDChartReport_All_3_Down_2026-07-08.csv",
}

# The specific 5 terms selected per direction (relevant to SCI and regeneration)
SELECTED_TERMS = {
    "Up": [
        "cellular response to L-leucine",
        "insulin-like growth factor receptor signaling pathway",
        "regulation of neuronal synaptic plasticity",
        "G1/S transition of mitotic cell cycle",
        "positive regulation of ERK1 and ERK2 cascade"
    ],
    "Down": [
        "axon extension",
        "actin cytoskeleton organization",
        "cell differentiation",
        "positive regulation of cell size",
        "regulation of multicellular organism growth"
    ]
}

CMAPS = {
    "Up":   "Purples",
    "Down": "Blues",
}

def extract_go_id(full_term):
    if isinstance(full_term, str) and "~" in full_term:
        return full_term.split("~")[0]
    return full_term

def capitalize_first(term_str):
    if not term_str:
        return ""
    # Capitalize ONLY the first character to preserve internal uppercase letters (like ERK1, L-leucine)
    return term_str[0].upper() + term_str[1:]

# ── Process Data ──────────────────────────────────────────────────────────────
dfs = {}
color_lists = {}

for direction, filename in FILES.items():
    filepath = os.path.join(INPUT_DIR, filename)
    if not os.path.exists(filepath):
        print(f"Warning: File {filepath} not found.")
        continue
    df = pd.read_csv(filepath)
    if df.empty:
        print(f"Warning: File {filepath} is empty.")
        continue
        
    df["-log10p"] = -np.log10(df["P-Value"].clip(lower=1e-300))
    selected = SELECTED_TERMS[direction]
    df_filtered = df[df["Term"].isin(selected)].copy()
    df_sorted = df_filtered.sort_values("Fold Enrichment", ascending=True)
    df_sorted["GO_ID"] = df_sorted["Full Term"].apply(extract_go_id)
    # Capitalize first letter and wrap to max 2 lines (width=30 ensures max 2 lines for all)
    df_sorted["Wrapped_Term"] = df_sorted["Term"].apply(lambda x: "\n".join(textwrap.wrap(capitalize_first(x), width=30)))
    
    # Color mapping
    cmap = plt.get_cmap(CMAPS[direction])
    p_vals = df_sorted["-log10p"].values
    vmin = max(p_vals.min() - 0.1, 0)
    vmax = p_vals.max() + 0.1
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    colors = [cmap(norm(v)) for v in p_vals]
    
    dfs[direction] = df_sorted
    color_lists[direction] = colors

# ── Generate Plots for both GOID and Names ────────────────────────────────────
for label_type in ["GOID", "Names"]:
    fig = plt.figure(figsize=(2.7, 3.5), dpi=300)
    ax = fig.add_axes([0.05, 0.10, 0.90, 0.82])
    
    y_down = np.arange(5)
    y_up = np.arange(5, 10)
    
    df_up = dfs["Up"]
    colors_up = color_lists["Up"]
    df_down = dfs["Down"]
    colors_down = color_lists["Down"]
    
    # Draw UP bars
    bars_up = ax.barh(
        y=y_up,
        width=df_up["Fold Enrichment"],
        color=colors_up,
        height=0.55,
        edgecolor="#333333",
        linewidth=0.3,
        zorder=2
    )
    
    # Draw DOWN bars
    bars_down = ax.barh(
        y=y_down,
        width=-df_down["Fold Enrichment"],
        color=colors_down,
        height=0.55,
        edgecolor="#333333",
        linewidth=0.3,
        zorder=2
    )
    
    # Draw central y-axis line at x=0
    ax.axvline(0, color='#000000', linewidth=0.6, zorder=3)
    ax.axhline(4.5, color='#cccccc', linestyle=':', linewidth=0.5, zorder=1)
    
    # Text Placement
    for yi, go_id, term in zip(y_up, df_up["GO_ID"], df_up["Wrapped_Term"]):
        lbl = go_id if label_type == "GOID" else term
        offset = -3.5 if label_type == "Names" else -2.0
        ax.text(
            offset, yi, lbl,
            ha='right', va='center',
            fontsize=3.4 if label_type == "Names" else 4.2, 
            fontweight='bold', color='#1a1a1a',
            linespacing=0.9
        )
        
    for yi, go_id, term in zip(y_down, df_down["GO_ID"], df_down["Wrapped_Term"]):
        lbl = go_id if label_type == "GOID" else term
        offset = 2.5 if label_type == "Names" else 2.0
        ax.text(
            offset, yi, lbl,
            ha='left', va='center',
            fontsize=3.4 if label_type == "Names" else 4.2, 
            fontweight='bold', color='#1a1a1a',
            linespacing=0.9
        )
        
    # Headers
    ax.text(0.25, 1.02, "Downregulated", transform=ax.transAxes, ha='center', va='bottom', 
            fontsize=5.0, fontweight='bold', color='#1d3557', clip_on=False)
    ax.text(0.75, 1.02, "Upregulated", transform=ax.transAxes, ha='center', va='bottom', 
            fontsize=5.0, fontweight='bold', color='#5a189a', clip_on=False)
            
    # Expanded x-limit to 140 as requested
    ax.set_xlim(-140, 140)
    ax.set_ylim(-0.6, 9.6)
    
    # Custom ticks for 140 range
    x_ticks = [-120, -60, 0, 60, 120]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(abs(t)) for t in x_ticks], fontsize=4.0)
    ax.set_yticks([])
    
    ax.set_xlabel("Fold Enrichment", fontsize=4.5, fontweight="bold", labelpad=2)
    ax.tick_params(axis='x', pad=1)
    
    for spine in ["top", "bottom", "left", "right"]:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_color('#000000')
        ax.spines[spine].set_linewidth(0.4)
        
    ax.xaxis.grid(True, linestyle=":", alpha=0.15, linewidth=0.25)
    ax.set_axisbelow(True)
    
    # Save Plot
    if label_type == "GOID":
        out_name = "DAVID_GO_All_3_stacked_GOID.png"
    else:
        out_name = "DAVID_GO_All_3_stacked_Names.png"
        
    if label_type == "Names":
        shutil_path = os.path.join(DAVID_DIR, "DAVID_GO_All_3_stacked.png")
        
    butterfly_path = os.path.join(OUTPUT_DIR, "Figure_X_F_Common_GO_butterfly.png")
    plt.savefig(butterfly_path, dpi=600, bbox_inches="tight", pad_inches=0.01)
    plt.close()
    
    if label_type == "Names":
        import shutil
        shutil.copy(butterfly_path, shutil_path)
        
    print(f"Saved combined butterfly barplot ({label_type}) to {butterfly_path}")
