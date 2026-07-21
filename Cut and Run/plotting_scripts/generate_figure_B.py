import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define paths
input_dir = '../input_data'
output_dir = '../output_figures'

rara_path = os.path.join(input_dir, 'Peaks_Annotated_RARA_vs_Input_p1e4.csv')
rarg_path = os.path.join(input_dir, 'Peaks_Annotated_RARG_vs_Input_p1e4.csv')

# Load files
df_rara = pd.read_csv(rara_path)
df_rarg = pd.read_csv(rarg_path)

# Helper function to group annotations
def group_annotation(ann):
    if not isinstance(ann, str):
        return 'Other'
    ann = ann.strip()
    if ann.startswith('Promoter'):
        return 'Promoter'
    elif ann.startswith('Exon'):
        return 'Exon'
    elif ann.startswith('Intron'):
        return 'Intron'
    elif 'UTR' in ann:
        return 'UTR'
    elif ann.startswith('Distal Intergenic'):
        return 'Distal Intergenic'
    elif ann.startswith('Downstream'):
        return 'Downstream'
    else:
        return 'Other'

df_rara['Group'] = df_rara['annotation'].apply(group_annotation)
df_rarg['Group'] = df_rarg['annotation'].apply(group_annotation)

# Categories list in stacking order
categories = ['Promoter', 'Exon', 'Intron', 'UTR', 'Distal Intergenic']

# Calculate counts and percentages
counts = {'RARA': {}, 'RARG': {}}
total_peaks = {'RARA': len(df_rara), 'RARG': len(df_rarg)}

for cat in categories:
    counts['RARA'][cat] = sum(df_rara['Group'] == cat)
    counts['RARG'][cat] = sum(df_rarg['Group'] == cat)

# Convert to percentages
percentages = {'RARA': [], 'RARG': []}

for cat in categories:
    percentages['RARA'].append((counts['RARA'][cat] / total_peaks['RARA']) * 100)
    percentages['RARG'].append((counts['RARG'][cat] / total_peaks['RARG']) * 100)

# Setup plot
fig, ax = plt.subplots(figsize=(2, 2), dpi=300)

# Colors: nature style palette
colors = {
    'Promoter': '#4C72B0',           # Muted Blue
    'Exon': '#55A868',               # Muted Green
    'Intron': '#C44E52',             # Muted Red
    'UTR': '#8172B3',                # Muted Purple
    'Distal Intergenic': '#CCB974'   # Muted Gold
}

# Plot bars
rara_bottom = 0
rarg_bottom = 0

bar_width = 0.5
x_positions = [0, 1]

for i, cat in enumerate(categories):
    rara_pct = percentages['RARA'][i]
    rarg_pct = percentages['RARG'][i]
    
    # Plot segment
    ax.bar(0, rara_pct, bottom=rara_bottom, width=bar_width, color=colors[cat])
    ax.bar(1, rarg_pct, bottom=rarg_bottom, width=bar_width, color=colors[cat])
    
    # Label RARA segment with percentage if height is reasonable (> 5%)
    if rara_pct > 5:
        ax.text(0, rara_bottom + rara_pct / 2, f"{rara_pct:.1f}%", 
                ha='center', va='center', color='white', fontsize=4.8, weight='bold')
        
    # Label RARG segment with percentage if height is reasonable (> 5%)
    if rarg_pct > 5:
        ax.text(1, rarg_bottom + rarg_pct / 2, f"{rarg_pct:.1f}%", 
                ha='center', va='center', color='white', fontsize=4.8, weight='bold')
        
    rara_bottom += rara_pct
    rarg_bottom += rarg_pct

# Customize axes
ax.set_ylabel('% of peaks', fontsize=7, labelpad=2)
ax.set_xticks(x_positions)
ax.set_xticklabels(["RARA", "RARG"], fontsize=7)
ax.set_ylim(0, 100)
ax.tick_params(axis='both', which='major', labelsize=6, pad=1, length=2)

# Make it a clean closed box (all spines visible)
for spine in ['top', 'right', 'left', 'bottom']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#1a202c')
    ax.spines[spine].set_linewidth(0.6)

# Add Legend at the bottom in 3 columns (creates 2 rows for 5 items) - slightly bigger font
legend_handles = [plt.Rectangle((0,0),1,1, color=colors[cat]) for cat in categories]
ax.legend(legend_handles, categories, loc='upper center', bbox_to_anchor=(0.5, -0.20),
          ncol=3, fontsize=5.5, frameon=False, handlelength=0.9, handleheight=0.9, 
          columnspacing=0.8, labelspacing=0.4)

# Save figure
plt.savefig(os.path.join(output_dir, 'FigureX_B.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Figure X.B updated with percentages and formatted ticks.")
