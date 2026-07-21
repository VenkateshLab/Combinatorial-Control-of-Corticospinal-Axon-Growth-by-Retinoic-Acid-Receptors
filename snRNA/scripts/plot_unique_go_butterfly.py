import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_DIR = "../raw_data"
OUTPUT_DIR = "../figures"
DAVID_DIR = INPUT_DIR

FILES = {
    "RARA_Up":   "DAVIDChartReport_UP_RARA_2026-07-08.csv",
    "RARG_Up":   "DAVIDChartReport_UP_RARG_2026-07-08.csv",
    "Combo_Up":  "DAVIDChartReport_Combo_up_2026-07-08.csv",
    "RARA_Down": "DAVIDChartReport_Down_RARA_2026-07-08.csv",
    "RARG_Down": "DAVIDChartReport_DOwn_RARG_2026-07-08.csv",
    "Combo_Down":"DAVIDChartReport_Combo_Down_2026-07-08.csv"
}

# Graded Purple Palette for Upregulated unique pathways
UP_COLORS = {
    "RARA":  "#cca8fa",  # Light Purple
    "RARG":  "#9e68e4",  # Medium Purple
    "Combo": "#7030c0"   # Dark Purple
}

# Graded Grey/Slate Palette for Downregulated unique pathways
DOWN_COLORS = {
    "RARA":  "#cbd5e1",  # Light Grey
    "RARG":  "#94a3b8",  # Medium Grey
    "Combo": "#475569"   # Dark Grey
}

def extract_go_id(full_term):
    if isinstance(full_term, str) and "~" in full_term:
        return full_term.split("~")[0]
    return full_term

def clean_term(full_term):
    if isinstance(full_term, str) and "~" in full_term:
        return full_term.split("~")[1]
    return full_term

def capitalize_first(term_str):
    if not term_str:
        return ""
    return term_str[0].upper() + term_str[1:]

def load_top_terms(filename, condition_name, max_terms=3):
    path = os.path.join(INPUT_DIR, filename)
    if not os.path.exists(path):
        print(f"Warning: {path} not found.")
        return pd.DataFrame()
    df = pd.read_csv(path)
    if df.empty:
        return pd.DataFrame()
        
    df["GO_ID"] = df["Full Term"].apply(extract_go_id)
    df["Term_Clean"] = df["Full Term"].apply(clean_term)
    
    # Sort by P-Value ascending (most significant first)
    df_sorted = df.sort_values("P-Value", ascending=True)
    df_top = df_sorted.head(max_terms).copy()
    
    # We want to reverse the order for plotting from bottom-to-top in horizontal bar plot
    df_top = df_top.iloc[::-1]
    df_top["Condition"] = condition_name
    return df_top

# ── Load data ─────────────────────────────────────────────────────────────────
# UP Unique
df_rara_up  = load_top_terms(FILES["RARA_Up"], "RARA", max_terms=2) # Only 2 significant terms available
df_rarg_up  = load_top_terms(FILES["RARG_Up"], "RARG", max_terms=3)
df_combo_up = load_top_terms(FILES["Combo_Up"], "Combo", max_terms=3)

# DOWN Unique
df_rara_down  = load_top_terms(FILES["RARA_Down"], "RARA", max_terms=3)
df_rarg_down  = load_top_terms(FILES["RARG_Down"], "RARG", max_terms=3)
df_combo_down = load_top_terms(FILES["Combo_Down"], "Combo", max_terms=3)

# ── Align rows by condition ───────────────────────────────────────────────────
# Combo unique: Row 0, 1, 2 (at the bottom)
# RARG unique: Row 3, 4, 5 (in the middle)
# RARA unique: Row 6, 7, 8 (at the top, with Row 8 UP being empty)

build_condition_list = lambda df: df.to_dict("records") + [None]*(3-len(df)) if not df.empty else [None]*3
combo_up_list = build_condition_list(df_combo_up)
rarg_up_list = build_condition_list(df_rarg_up)
rara_up_list = build_condition_list(df_rara_up)
up_rows = combo_up_list + rarg_up_list + rara_up_list

combo_down_list = build_condition_list(df_combo_down)
rarg_down_list = build_condition_list(df_rarg_down)
rara_down_list = build_condition_list(df_rara_down)
down_rows = combo_down_list + rarg_down_list + rara_down_list

# ── Generate Plots for both GOID and Names ────────────────────────────────────
for label_type in ["GOID", "Names"]:
    # 2.9 x 3.5 inches
    fig = plt.figure(figsize=(2.9, 3.5), dpi=300)
    
    # Slightly larger left margin to accommodate the facet strip
    ax = fig.add_axes([0.08, 0.10, 0.88, 0.82])
    
    y_indices = np.arange(9)
    
    # Draw Bars and Text labels
    for y in y_indices:
        up_data = up_rows[y]
        down_data = down_rows[y]
        
        # 1. Draw UP bar (Positive side, extending right)
        if up_data is not None:
            w_up = up_data["Fold Enrichment"]
            cond_up = up_data["Condition"]
            go_id_up = up_data["GO_ID"]
            
            # Capitalize and wrap to maximum 2 lines (width=30 ensures maximum 2 lines)
            wrapped_term = "\n".join(textwrap.wrap(capitalize_first(up_data["Term_Clean"]), width=30))
            lbl = go_id_up if label_type == "GOID" else wrapped_term
            
            ax.barh(
                y=y,
                width=w_up,
                color=UP_COLORS[cond_up],
                height=0.55,
                edgecolor="#333333",
                linewidth=0.3,
                zorder=2
            )
            
            # Place label to the right of the bar tip
            ax.text(
                w_up + 2.5, y, lbl,
                ha='left', va='center',
                fontsize=3.2 if label_type == "Names" else 4.2, 
                fontweight='bold', color='#1a1a1a',
                linespacing=0.9
            )
            
        # 2. Draw DOWN bar (Negative side, extending left)
        if down_data is not None:
            w_down = down_data["Fold Enrichment"]
            cond_down = down_data["Condition"]
            go_id_down = down_data["GO_ID"]
            
            # Capitalize and wrap to maximum 2 lines (width=30 ensures maximum 2 lines)
            wrapped_term = "\n".join(textwrap.wrap(capitalize_first(down_data["Term_Clean"]), width=30))
            lbl = go_id_down if label_type == "GOID" else wrapped_term
            
            ax.barh(
                y=y,
                width=-w_down,
                color=DOWN_COLORS[cond_down],
                height=0.55,
                edgecolor="#333333",
                linewidth=0.3,
                zorder=2
            )
            
            # Place label to the left of the bar tip
            ax.text(
                -w_down - 2.5, y, lbl,
                ha='right', va='center',
                fontsize=3.2 if label_type == "Names" else 4.2, 
                fontweight='bold', color='#1a1a1a',
                linespacing=0.9
            )

    # Draw central y-axis line at x=0
    ax.axvline(0, color='#000000', linewidth=0.6, zorder=3)
    
    # Draw light grey dotted lines separating categories (only in the plotting area x > -150)
    ax.hlines(2.5, -150, 140, colors='#cccccc', linestyles=':', linewidths=0.5, zorder=1)
    ax.hlines(5.5, -150, 140, colors='#cccccc', linestyles=':', linewidths=0.5, zorder=1)

    # Draw Left Vertical facet-like category strips (occupying x = [-165, -150])
    # Combo unique: Row 0, 1, 2 (y range: [-0.5, 2.5])
    # RARG unique: Row 3, 4, 5 (y range: [2.5, 5.5])
    # RARA unique: Row 6, 7, 8 (y range: [5.5, 8.5])
    
    # Draw grey backgrounds for the strips
    rect_combo = mpatches.Rectangle((-165, -0.5), 15, 3.0, facecolor='#f1f5f9', edgecolor='#000000', linewidth=0.4, zorder=3)
    rect_rarg  = mpatches.Rectangle((-165, 2.5), 15, 3.0, facecolor='#f1f5f9', edgecolor='#000000', linewidth=0.4, zorder=3)
    rect_rara  = mpatches.Rectangle((-165, 5.5), 15, 3.0, facecolor='#f1f5f9', edgecolor='#000000', linewidth=0.4, zorder=3)
    
    ax.add_patch(rect_combo)
    ax.add_patch(rect_rarg)
    ax.add_patch(rect_rara)
    
    # Write Category Names inside the strips (rotated 90 degrees)
    ax.text(-157.5, 1.0, "Combo Unique", ha='center', va='center', rotation=90, fontsize=4.0, fontweight='bold', color='#1a1a1a', zorder=4)
    ax.text(-157.5, 4.0, "RARG Unique", ha='center', va='center', rotation=90, fontsize=4.0, fontweight='bold', color='#1a1a1a', zorder=4)
    ax.text(-157.5, 7.0, "RARA Unique", ha='center', va='center', rotation=90, fontsize=4.0, fontweight='bold', color='#1a1a1a', zorder=4)

    # Section Headers placed OUTSIDE the plot box (using transAxes coordinates, adjusted for shifted plot center)
    ax.text(0.28, 1.02, "Downregulated Unique", transform=ax.transAxes, ha='center', va='bottom', 
            fontsize=5.0, fontweight='bold', color='#475569', clip_on=False)
    ax.text(0.76, 1.02, "Upregulated Unique", transform=ax.transAxes, ha='center', va='bottom', 
            fontsize=5.0, fontweight='bold', color='#7030c0', clip_on=False)

    # Expanded x-limit to 140 on positive and -165 on negative to fit the left facet strip
    ax.set_xlim(-165, 140)
    ax.set_ylim(-0.6, 8.6)

    # Custom ticks (positive-only) inside the plotting area (x >= -150)
    x_ticks = [-120, -60, 0, 60, 120]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(abs(t)) for t in x_ticks], fontsize=4.0)
    ax.set_yticks([])  # Hide default y-ticks

    ax.set_xlabel("Fold Enrichment", fontsize=4.5, fontweight="bold", labelpad=2)
    ax.tick_params(axis='x', pad=1)

    # Boxed layout: show all 4 spines surrounding the plot area (including the facet strip)
    for spine in ["top", "bottom", "left", "right"]:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_color('#000000')
        ax.spines[spine].set_linewidth(0.4)

    # Subtle grid lines
    ax.xaxis.grid(True, linestyle=":", alpha=0.15, linewidth=0.25)
    ax.set_axisbelow(True)

    # Save Plot
    if label_type == "GOID":
        out_name = "DAVID_GO_Unique_stacked_GOID.png"
    else:
        out_name = "DAVID_GO_Unique_stacked_Names.png"
        
    if label_type == "Names":
        shutil_path = os.path.join(DAVID_DIR, "DAVID_GO_Unique_stacked.png")
        
    out_plot_path = os.path.join(OUTPUT_DIR, "Figure_X_H_Unique_GO_butterfly.png")
    plt.savefig(out_plot_path, dpi=600, bbox_inches="tight", pad_inches=0.01)
    plt.close()
    
    if label_type == "Names":
        import shutil
        # shutil.copy(out_plot_path, shutil_path)
        
    print(f"Saved combined unique GO butterfly barplot ({label_type}) to {out_plot_path}")
