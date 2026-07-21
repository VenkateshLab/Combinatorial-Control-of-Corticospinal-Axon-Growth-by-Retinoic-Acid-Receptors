import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# File paths in workspace
INPUT_DIR = "../raw_data"
OUTPUT_DIR = "../figures"
files = {
    "RARA vs GFP": os.path.join(working_dir, "DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv"),
    "RARG vs GFP": os.path.join(working_dir, "DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv"),
    "Combo vs GFP": os.path.join(working_dir, "DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv")
}

# Thresholds
P_VAL_CUTOFF = 0.05
LOG2FC_CUTOFF = 0.5

# Calculate counts
results = []
for name, filepath in files.items():
    if not os.path.exists(filepath):
        print(f"Warning: File {filepath} not found.")
        continue
    df = pd.read_csv(filepath)
    up_genes = df[(df['p_val'] < P_VAL_CUTOFF) & (df['avg_log2FC'] >= LOG2FC_CUTOFF)]
    down_genes = df[(df['p_val'] < P_VAL_CUTOFF) & (df['avg_log2FC'] <= -LOG2FC_CUTOFF)]
    
    up_count = len(up_genes)
    down_count = len(down_genes)
    total = up_count + down_count
    
    results.append({
        "Dataset": name,
        "Up": up_count,
        "Down": down_count,
        "Total": total,
        "Up_Pct": (up_count / total) * 100 if total > 0 else 0,
        "Down_Pct": (down_count / total) * 100 if total > 0 else 0
    })

df_res = pd.DataFrame(results)
print("Calculated Proportions:")
print(df_res.to_string(index=False))

# Plot configuration
fig, ax = plt.subplots(figsize=(2.0, 1.5), dpi=300)

# Colors matching the previous plots: 
# Purple for Upregulated, Slate Grey for Downregulated
color_up = '#9e68e4'
color_down = '#94a3b8'

# X positions for the bars
x = np.arange(len(df_res))
bar_width = 0.55

# Plot stacked bars (Down at the bottom, Up stacked on top of Down)
bars_down = ax.bar(x, df_res["Down_Pct"], color=color_down, label="Downregulated", width=bar_width, edgecolor='none')
bars_up = ax.bar(x, df_res["Up_Pct"], bottom=df_res["Down_Pct"], color=color_up, label="Upregulated", width=bar_width, edgecolor='none')

# Add annotations inside the bar segments
for i in range(len(df_res)):
    down_pct = df_res.loc[i, "Down_Pct"]
    up_pct = df_res.loc[i, "Up_Pct"]
    down_val = df_res.loc[i, "Down"]
    up_val = df_res.loc[i, "Up"]
    
    # Downregulated label centered inside bottom segment
    ax.text(i, down_pct / 2.0, f"{down_pct:.1f}%\n({down_val:,})", 
            ha='center', va='center', fontsize=4.0, fontweight='bold', color='#ffffff', linespacing=1.0)
    
    # Upregulated label centered inside top segment
    ax.text(i, down_pct + (up_pct / 2.0), f"{up_pct:.1f}%\n({up_val:,})", 
            ha='center', va='center', fontsize=4.0, fontweight='bold', color='#ffffff', linespacing=1.0)

# Customizing Axes
ax.set_ylabel("DEG Proportion (%)", fontsize=5.5, fontweight='bold', color='#000000', labelpad=2)
ax.set_xticks(x)
ax.set_xticklabels(df_res["Dataset"], fontsize=5.0, fontweight='bold', color='#000000')
ax.tick_params(axis='x', pad=2)
ax.tick_params(axis='y', labelsize=4.5, colors='#000000', pad=1)
ax.set_ylim(0, 100)

# Borders
for spine in ['top', 'bottom', 'left', 'right']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#000000')
    ax.spines[spine].set_linewidth(0.4)

# Legend placement at the top
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False, fontsize=4.8, handletextpad=0.2, columnspacing=0.8)

# Subtle horizontal grid
ax.grid(axis='y', linestyle=':', color='#000000', alpha=0.15, linewidth=0.25)

plt.tight_layout(pad=0.2)

# Save Plot
output_plot_path = os.path.join(working_dir, "deg_proportions_plot.png")
plt.savefig(output_plot_path, dpi=600, bbox_inches='tight')
plt.close()
print(f"Stacked proportion plot saved to {output_plot_path}")
