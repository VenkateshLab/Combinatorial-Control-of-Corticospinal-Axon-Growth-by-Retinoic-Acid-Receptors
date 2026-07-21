import os
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define paths
input_dir = '../input_data'
output_dir = '../output_figures'

rara_path = os.path.join(input_dir, 'Peaks_Annotated_RARA_vs_Input_p1e4.csv')
rarg_path = os.path.join(input_dir, 'Peaks_Annotated_RARG_vs_Input_p1e4.csv')
genes_path = '../input_data/overlaps_p_1e-4_genes.txt'

# Load data
df_rara = pd.read_csv(rara_path)
df_rarg = pd.read_csv(rarg_path)

with open(genes_path, 'r') as f:
    raw_genes = [line.strip() for line in f if line.strip()]

# Filter out Rik and Mir genes
genes = []
for g in raw_genes:
    gl = g.lower()
    if gl.endswith('rik') or gl.startswith('mir'):
        continue
    genes.append(g)

# Extract fold enrichment for co-bound genes
data = []
for g in genes:
    rara_sub = df_rara[df_rara['SYMBOL'] == g]
    rarg_sub = df_rarg[df_rarg['SYMBOL'] == g]
    
    rara_fe = rara_sub['fold_enrichment'].max() if len(rara_sub) > 0 else np.nan
    rarg_fe = rarg_sub['fold_enrichment'].max() if len(rarg_sub) > 0 else np.nan
    
    data.append({
        'gene': g,
        'rara_fe': rara_fe,
        'rarg_fe': rarg_fe,
        'rara_log2': np.log2(rara_fe),
        'rarg_log2': np.log2(rarg_fe),
        'max_val': max(rara_fe, rarg_fe)
    })

df_plot = pd.DataFrame(data).dropna()

# Sort genes by their maximum fold enrichment for a clean staircase effect
df_plot = df_plot.sort_values(by='max_val', ascending=True).reset_index(drop=True)

# Setup plot
fig, ax = plt.subplots(figsize=(2, 2.2), dpi=300)

# Colors: exact hex codes picked by the user
colors_face = {
    'RARA': '#899cbf',
    'RARG': '#da9497'
}
colors_edge = {
    'RARA': '#677fac',
    'RARG': '#cf7579'
}

# Find index of Rarg to highlight its row
rarg_idx = df_plot[df_plot['gene'] == 'Rarg'].index[0]

# Add a subtle light-pink background highlight strip behind Rarg (matching RARG color theme)
ax.axhspan(rarg_idx - 0.45, rarg_idx + 0.45, color='#ffebee', alpha=0.8, zorder=0, edgecolor='none')

# Plot dumbbells
for idx, row in df_plot.iterrows():
    # Draw horizontal connector line
    ax.plot([row['rara_log2'], row['rarg_log2']], [idx, idx], color='#cbd5e0', linewidth=0.8, zorder=1)
    
    # Draw RARA dot with exact face and edge colors
    ax.scatter(row['rara_log2'], idx, 
               facecolors=colors_face['RARA'], 
               edgecolors=colors_edge['RARA'], 
               linewidths=0.5, s=14, zorder=2)
    
    # Draw RARG dot with exact face and edge colors
    ax.scatter(row['rarg_log2'], idx, 
               facecolors=colors_face['RARG'], 
               edgecolors=colors_edge['RARG'], 
               linewidths=0.5, s=14, zorder=2)

# Customize axes
ax.set_xlabel('Fold enrichment', fontsize=7, labelpad=2)
ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot['gene'], fontsize=5.0)

# Apply italic styling for all genes, and bold + color highlight specifically for Rarg
for tick_label in ax.get_yticklabels():
    gene_name = tick_label.get_text()
    if gene_name == 'Rarg':
        tick_label.set_style('italic')
        tick_label.set_weight('bold')
        tick_label.set_color(colors_edge['RARG'])  # Highlight in Crimson Red outline color
    else:
        tick_label.set_style('italic')

ax.tick_params(axis='both', which='major', labelsize=5.5, pad=1, length=1.5)

# Define X-axis ticks (log2 values labeled as original values)
x_ticks_log = [0, 1, 2, 3, 4, 5, 6, 7]
x_ticks_labels = ['1', '2', '4', '8', '16', '32', '64', '128']
ax.set_xticks(x_ticks_log)
ax.set_xticklabels(x_ticks_labels, fontsize=5.5)
ax.set_xlim(-0.5, 7.5)

# Make it a clean closed box
for spine in ['top', 'right', 'left', 'bottom']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#1a202c')
    ax.spines[spine].set_linewidth(0.6)

# Add Legend at the bottom
legend_handles = [
    plt.Line2D([0], [0], marker='o', color='w', 
               markerfacecolor=colors_face['RARA'], 
               markeredgecolor=colors_edge['RARA'], 
               markeredgewidth=0.5, markersize=5, label='RARA'),
    plt.Line2D([0], [0], marker='o', color='w', 
               markerfacecolor=colors_face['RARG'], 
               markeredgecolor=colors_edge['RARG'], 
               markeredgewidth=0.5, markersize=5, label='RARG')
]
ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.16), 
          ncol=2, fontsize=5.5, frameon=False, handletextpad=0.2, columnspacing=1.0)

# Save figure
plt.savefig(os.path.join(output_dir, 'FigureX_E.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Figure X.E updated with exact custom colors.")
