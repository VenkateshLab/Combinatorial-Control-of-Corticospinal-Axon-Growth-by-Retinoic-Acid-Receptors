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

act_path = '../input_data/Combo_Direct_Activated.csv'
rep_path = '../input_data/Combo_Direct_Repressed.csv'

# Load data
df_rara = pd.read_csv(rara_path)
df_rarg = pd.read_csv(rarg_path)

df_act = pd.read_csv(act_path)
df_rep = pd.read_csv(rep_path)

# Filter out Rik and Mir genes
df_act = df_act[~df_act['gene'].str.lower().str.endswith('rik') & ~df_act['gene'].str.lower().str.startswith('mir')].copy()
df_rep = df_rep[~df_rep['gene'].str.lower().str.endswith('rik') & ~df_rep['gene'].str.lower().str.startswith('mir')].copy()

# Sort activated genes by avg_log2FC ascending
df_act = df_act.sort_values(by='avg_log2FC', ascending=True).reset_index(drop=True)

# Sort repressed genes by avg_log2FC ascending
df_rep = df_rep.sort_values(by='avg_log2FC', ascending=True).reset_index(drop=True)

# Combine them: repressed at bottom, activated at top
df_rep['Type'] = 'Repressed'
df_act['Type'] = 'Activated'
df_all = pd.concat([df_rep, df_act]).reset_index(drop=True)

# Lookup RARA and RARG peak strengths for each gene
plot_data = []
for idx, row in df_all.iterrows():
    g = row['gene']
    
    rara_sub = df_rara[df_rara['SYMBOL'] == g]
    rarg_sub = df_rarg[df_rarg['SYMBOL'] == g]
    
    # Use raw fold enrichment (0 if no peak present)
    rara_fe = rara_sub['fold_enrichment'].max() if len(rara_sub) > 0 else 0
    rarg_fe = rarg_sub['fold_enrichment'].max() if len(rarg_sub) > 0 else 0
    
    plot_data.append({
        'gene': g,
        'rara_fe': rara_fe,
        'rarg_fe': rarg_fe,
        'log2FC': row['avg_log2FC'],
        'type': row['Type']
    })

df_plot = pd.DataFrame(plot_data)

# Setup plot
fig, ax = plt.subplots(figsize=(2, 2.5), dpi=300)

# Colors: exact custom hex codes
colors_face = {
    'RARA': '#899cbf',
    'RARG': '#da9497'
}
colors_edge = {
    'RARA': '#677fac',
    'RARG': '#cf7579'
}

# Group annotation colors: Green for Activated, Slate Grey for Repressed
group_colors = {
    'Activated': '#2e7d32',
    'Repressed': '#718096'
}

# Find index of Rarg to highlight its row
rarg_rows = df_plot[df_plot['gene'] == 'Rarg']
if len(rarg_rows) > 0:
    rarg_idx = rarg_rows.index[0]
    # Add a subtle light-pink background highlight strip behind Rarg
    ax.axhspan(rarg_idx - 0.45, rarg_idx + 0.45, color='#ffebee', alpha=0.8, zorder=0, edgecolor='none')

# Plot butterfly bars using raw fold enrichment (linear scale)
for idx, row in df_plot.iterrows():
    text_color = group_colors[row['type']]
    fc_val = row['log2FC']
    sign = '+' if fc_val > 0 else ''
    fc_text = f"{sign}{fc_val:.2f}"
    
    # Left bar: RARA peak fold enrichment (plotted negative)
    if row['rara_fe'] > 0:
        ax.barh(idx, -row['rara_fe'], 
                facecolor=colors_face['RARA'], 
                edgecolor=colors_edge['RARA'], 
                linewidth=0.5, height=0.6, zorder=2)
        
        # Write Combo log2FC value at the end of the RARA bar (on the left)
        ax.text(-row['rara_fe'] - 2.5, idx, fc_text, ha='right', va='center', fontsize=3.8, color=text_color, weight='bold')
    
    # Right bar: RARG peak fold enrichment (plotted positive)
    if row['rarg_fe'] > 0:
        ax.barh(idx, row['rarg_fe'], 
                facecolor=colors_face['RARG'], 
                edgecolor=colors_edge['RARG'], 
                linewidth=0.5, height=0.6, zorder=2)
        
        # Write Combo log2FC value at the end of the RARG bar (on the right)
        ax.text(row['rarg_fe'] + 2.5, idx, fc_text, ha='left', va='center', fontsize=3.8, color=text_color, weight='bold')

# Draw vertical divider at x=0
ax.axvline(0, color='#1a202c', linewidth=0.8, zorder=3)

# Draw horizontal divider between Repressed and Activated sections
num_repressed = len(df_rep)
ax.axhline(num_repressed - 0.5, color='#7f7f7f', linestyle='--', linewidth=0.5, zorder=1)

# Annotate sections on the same side (left side at x = -75)
ax.text(-75, num_repressed - 0.4, 'Activated', fontsize=4.8, weight='bold', color=group_colors['Activated'], va='bottom', ha='left')
ax.text(-75, num_repressed - 0.6, 'Repressed', fontsize=4.8, weight='bold', color=group_colors['Repressed'], va='top', ha='left')

# Customize axes
ax.set_xlabel('Peak Fold Enrichment', fontsize=7, labelpad=2)
ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot['gene'], fontsize=5.0)

# Apply italic styling for all genes, and bold + color highlight specifically for Rarg
for tick_label in ax.get_yticklabels():
    gene_name = tick_label.get_text()
    if gene_name == 'Rarg':
        tick_label.set_style('italic')
        tick_label.set_weight('bold')
        tick_label.set_color(colors_edge['RARG'])
    else:
        tick_label.set_style('italic')

ax.tick_params(axis='both', which='major', labelsize=5.5, pad=1, length=1.5)

# Define symmetric X-axis ticks on linear scale (from -75 to 75)
x_ticks = [-75, -50, -25, 0, 25, 50, 75]
x_ticks_labels = ['75', '50', '25', '0', '25', '50', '75']
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_ticks_labels, fontsize=5.5)
# Set X-limit to (-87, 87) for symmetric plotting with space for the text labels
ax.set_xlim(-87, 87)

# Make it a clean closed box
for spine in ['top', 'right', 'left', 'bottom']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#1a202c')
    ax.spines[spine].set_linewidth(0.6)

# Add Legend at the bottom left/center of the plot
legend_handles = [
    plt.Rectangle((0,0),1,1, facecolor=colors_face['RARA'], edgecolor=colors_edge['RARA'], linewidth=0.5, label='RARA'),
    plt.Rectangle((0,0),1,1, facecolor=colors_face['RARG'], edgecolor=colors_edge['RARG'], linewidth=0.5, label='RARG')
]
ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.4, -0.14), 
          ncol=2, fontsize=5.5, frameon=False, handlelength=0.8, handleheight=0.8, columnspacing=1.0)

# Save figure
plt.savefig(os.path.join(output_dir, 'FigureX_F.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Figure X.F updated with exact custom colors.")
