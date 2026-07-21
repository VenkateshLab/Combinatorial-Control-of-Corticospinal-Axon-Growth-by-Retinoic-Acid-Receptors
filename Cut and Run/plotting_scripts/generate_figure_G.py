import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# Define paths
input_dir = '../input_data'
output_dir = '../output_figures'

rara_path = os.path.join(input_dir, 'Peaks_Annotated_RARA_vs_Input_p1e4.csv')
rarg_path = os.path.join(input_dir, 'Peaks_Annotated_RARG_vs_Input_p1e4.csv')

rara_deg_path = '../input_data/DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv'
rarg_deg_path = '../input_data/DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv'
combo_deg_path = '../input_data/DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv'

# Load data
df_rara_peaks = pd.read_csv(rara_path)
df_rarg_peaks = pd.read_csv(rarg_path)

df_rara_deg = pd.read_csv(rara_deg_path)
df_rarg_deg = pd.read_csv(rarg_deg_path)
df_combo_deg = pd.read_csv(combo_deg_path)

# Extract unique bound gene lists
rara_bound = set(df_rara_peaks['SYMBOL'].dropna().unique())
rarg_bound = set(df_rarg_peaks['SYMBOL'].dropna().unique())
combo_bound = rara_bound.union(rarg_bound)

# Helper function to compute percentages and Fisher exact p-value
def get_stats(df_deg, bound_genes):
    total_bg = len(df_deg)
    bound_bg = len(df_deg[df_deg['gene'].isin(bound_genes)])
    bg_pct = (bound_bg / total_bg) * 100
    
    # Significant DEGs (unadjusted p < 0.05)
    sig_degs = df_deg[df_deg['p_val'] < 0.05]
    sig_up = sig_degs[sig_degs['avg_log2FC'] > 0]
    sig_down = sig_degs[sig_degs['avg_log2FC'] < 0]
    
    # UP DEGs
    up_total = len(sig_up)
    up_bound = len(sig_up[sig_up['gene'].isin(bound_genes)])
    up_pct = (up_bound / up_total) * 100 if up_total > 0 else 0.0
    
    # Fisher test for UP
    bg_other_total = total_bg - up_total
    bg_other_bound = bound_bg - up_bound
    table_up = [[up_bound, up_total - up_bound], [bg_other_bound, bg_other_total - bg_other_bound]]
    p_up = stats.fisher_exact(table_up, alternative='greater')[1]
    
    # DOWN DEGs
    down_total = len(sig_down)
    down_bound = len(sig_down[sig_down['gene'].isin(bound_genes)])
    down_pct = (down_bound / down_total) * 100 if down_total > 0 else 0.0
    
    # Fisher test for DOWN
    bg_other_total_dn = total_bg - down_total
    bg_other_bound_dn = bound_bg - down_bound
    table_dn = [[down_bound, down_total - down_bound], [bg_other_bound_dn, bg_other_total_dn - bg_other_bound_dn]]
    p_dn = stats.fisher_exact(table_dn, alternative='greater')[1]
    
    return {
        'bg_pct': bg_pct,
        'up_pct': up_pct,
        'down_pct': down_pct,
        'p_up': p_up,
        'p_dn': p_dn
    }

# Compute stats for RARA, RARG, and Combo
rara_stats = get_stats(df_rara_deg, rara_bound)
rarg_stats = get_stats(df_rarg_deg, rarg_bound)
combo_stats = get_stats(df_combo_deg, combo_bound)

# Prepare plotting data
categories = ['RARA Up', 'RARA Down', 'RARG Up', 'RARG Down', 'Combo Up', 'Combo Down']
deg_pcts = [
    rara_stats['up_pct'], rara_stats['down_pct'],
    rarg_stats['up_pct'], rarg_stats['down_pct'],
    combo_stats['up_pct'], combo_stats['down_pct']
]
p_values = [
    rara_stats['p_up'], rara_stats['p_dn'],
    rarg_stats['p_up'], rarg_stats['p_dn'],
    combo_stats['p_up'], combo_stats['p_dn']
]

# Set up matplotlib figure
fig, ax = plt.subplots(figsize=(2.2, 2.0), dpi=300)

x = np.arange(len(categories))
width = 0.55

# Custom color scheme: exact user picks (face & edge)
colors_face = {
    'RARA': '#899cbf',
    'RARG': '#da9497',
    'Combo': '#b5a4cd'
}
colors_edge = {
    'RARA': '#677fac',
    'RARG': '#cf7579',
    'Combo': '#9c89b8'
}

# Determine colors for each bar
bar_groups = ['RARA', 'RARA', 'RARG', 'RARG', 'Combo', 'Combo']

# Plot single bar per category representing DEGs with exact custom face and edge colors
rects = []
for i in range(len(categories)):
    group = bar_groups[i]
    r = ax.bar(x[i], deg_pcts[i], width, 
               facecolor=colors_face[group], 
               edgecolor=colors_edge[group], 
               linewidth=0.5, zorder=2)
    rects.append(r[0])

# Customize axes
ax.set_ylabel('% of DEGs directly bound', fontsize=7, labelpad=2)
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=90, ha='center', fontsize=5.5)
ax.set_ylim(0, 9.0)
ax.tick_params(axis='both', which='major', labelsize=5.5, pad=1, length=1.5)

# Add exact percentage text on top of each bar
for i, rect in enumerate(rects):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., height + 0.15, f"{height:.1f}%",
            ha='center', va='bottom', fontsize=3.8, color='#1a202c', weight='bold')
    
    # Add statistical significance star on top of the DEG bar if p < 0.05
    if p_values[i] < 0.05:
        ax.text(rect.get_x() + rect.get_width()/2., height + 0.9, '*',
                ha='center', va='center', fontsize=8, color='#cf7579', weight='bold') # Highlight with RARG outline color

# Closed box style
for spine in ['top', 'right', 'left', 'bottom']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#1a202c')
    ax.spines[spine].set_linewidth(0.6)

# Create custom legend to display at the bottom (below the 90-degree rotated labels)
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=colors_face['RARA'], edgecolor=colors_edge['RARA'], linewidth=0.5, label='RARA'),
    Patch(facecolor=colors_face['RARG'], edgecolor=colors_edge['RARG'], linewidth=0.5, label='RARG'),
    Patch(facecolor=colors_face['Combo'], edgecolor=colors_edge['Combo'], linewidth=0.5, label='Combo')
]
ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.45), 
          ncol=3, fontsize=5.5, title='DEGs', title_fontsize=5.5, frameon=False, 
          handlelength=0.8, handleheight=0.8, columnspacing=1.0)

# Grid lines
ax.grid(axis='y', linestyle=':', alpha=0.3, linewidth=0.4, zorder=1)

# Adjust layout and save
plt.savefig(os.path.join(output_dir, 'FigureX_G.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Figure X.G updated with custom colors.")
