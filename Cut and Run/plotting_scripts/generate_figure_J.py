import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors

# Define paths
rara_peaks_csv = '../input_data/Peaks_Annotated_RARA_vs_Input_p1e4.csv'
rarg_peaks_csv = 'Gene_Ontology/Selected/Peaks_AnnotRIC_RARG_vs_Input_p1e4.csv' # Wait, let's make sure it's the correct path: Gene_Ontology/Selected/Peaks_Annotated_RARG_vs_Input_p1e4.csv
# In previous version it was: Gene_Ontology/Selected/Peaks_Annotated_RARG_vs_Input_p1e4.csv. Let's make sure we use that!
rara_peaks_csv = '../input_data/Peaks_Annotated_RARA_vs_Input_p1e4.csv'
rarg_peaks_csv = '../input_data/Peaks_Annotated_RARG_vs_Input_p1e4.csv'

output_dir = '../output_figures'
os.makedirs(output_dir, exist_ok=True)

# 1. Load peak data
df_rara = pd.read_csv(rara_peaks_csv)
df_rarg = pd.read_csv(rarg_peaks_csv)

# Build case-insensitive lookup maps
rara_fe = {str(k).upper(): v for k, v in df_rara.groupby('SYMBOL')['fold_enrichment'].max().to_dict().items()}
rarg_fe = {str(k).upper(): v for k, v in df_rarg.groupby('SYMBOL')['fold_enrichment'].max().to_dict().items()}

# Helper to format gene symbols to mouse standard (Italicized, Title-cased)
def to_mouse_casing(gene):
    if len(gene) <= 1:
        return gene.upper()
    return gene[0].upper() + gene[1:].lower()

# Define the 4 GO pathways and their genes (split into two lines for readability)
# Replacing Hindlimb morphogenesis & Bone morphogenesis with CNS/SCI-relevant terms:
# 1. Axon Guidance & Repair (or Axonal Regeneration)
# 2. Cell Differentiation (or Neuronal Differentiation)
# 3. Nervous System Development
# 4. ECM & Neuroprotection (or ECM & Cellular Repair)

pathways = [
    {
        'name': 'Axon Guidance\n& Repair',
        'genes': ['RARG', 'CPEB2', 'RPS29', 'NEMF', 'ZFAND5']
    },
    {
        'name': 'Cell\ndifferentiation',
        'genes': ['ALKBH5', 'NR2F1', 'NR2F6', 'INSM1', 'CHD7']
    },
    {
        'name': 'Nervous system\ndevelopment',
        'genes': ['DIP2B', 'ARID1B', 'PRKDC', 'ABR', 'ALOX5AP']
    },
    {
        'name': 'ECM &\nNeuroprotection',
        'genes': ['SKI', 'COMP', 'ATG9B', 'EXT1', 'MECOM']
    }
]

# 2. Render Split Horizontal Figure (7.5 x 1.3 inches, 500 DPI)
fig = plt.figure(figsize=(7.5, 1.3), dpi=500)
# 4 columns for pathways + 1 column for colorbar
gs = gridspec.GridSpec(1, 5, width_ratios=[5, 5, 5, 5, 0.6], wspace=0.28,
                       left=0.08, right=0.92, top=0.68, bottom=0.28)

# Custom soft rose colormap (transitions from a distinct soft gray for zero values to soft rose-pink for high values)
cmap_custom = mcolors.LinearSegmentedColormap.from_list('soft_rose', ['#e5e7eb', '#dca9ab', '#da9497'])

im = None

# Iterate through pathways and draw their respective subplots
for idx, p in enumerate(pathways):
    ax = fig.add_subplot(gs[0, idx])
    
    # Build data matrix for this pathway
    p_data = []
    p_genes_mouse = []
    for g in p['genes']:
        rara_val = rara_fe.get(g.upper(), 0.0)
        rarg_val = rarg_fe.get(g.upper(), 0.0)
        p_data.append([rara_val, rarg_val])
        p_genes_mouse.append(to_mouse_casing(g))
        
    p_matrix = np.array(p_data).T  # Shape: (2, n_genes)
    
    # Plot heatmap block
    im = ax.imshow(p_matrix, cmap=cmap_custom, aspect='auto', interpolation='nearest', vmin=0, vmax=95)
    
    # Show values inside cells
    for i in range(p_matrix.shape[0]):
        for j in range(p_matrix.shape[1]):
            val = p_matrix[i, j]
            # Since the colormap is soft and light, black text has excellent contrast throughout
            val_str = f"{val:.1f}" if val > 0 else "-"
            ax.text(j, i, val_str, ha='center', va='center', fontsize=4.8, color='black')
            
    # X-axis configuration (gene names rotated 45 degrees, aligned to the right)
    ax.set_xticks(np.arange(len(p_genes_mouse)))
    ax.set_xticklabels(p_genes_mouse, fontsize=5.5, rotation=45, ha='right',
                       fontstyle='italic', weight='semibold', color='black')
    
    # Y-axis configuration (only show labels on the first subplot)
    if idx == 0:
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['RARA', 'RARG'], fontsize=7, weight='bold', color='black')
    else:
        ax.set_yticks([])
        
    # Set pathway name as subplot title (enlarged to 7.0 and in solid black text)
    ax.set_title(p['name'], fontsize=7.0, weight='bold', pad=5, color='black')
    
    # Configure spines/borders
    ax.spines['top'].set_linewidth(0.4)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.spines['left'].set_linewidth(0.4)
    ax.spines['right'].set_linewidth(0.4)
    ax.tick_params(axis='both', which='both', length=0) # Hide ticks

# Colorbar in the last slot
cax = fig.add_subplot(gs[0, 4])
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=5, width=0.4, length=2, colors='black')
cbar.outline.set_linewidth(0.4)
cax.set_ylabel('Peak Fold Enrichment', fontsize=5.5, labelpad=2, color='black')

# Save figure
out_png = os.path.join(output_dir, 'FigureX_J.png')
plt.savefig(out_png, dpi=500, bbox_inches='tight', transparent=True)
plt.close()
print(f"Generated Split Horizontal Figure J Heatmap: {out_png}")

# Copy to brain folder
artifact_dir = '/home/manojkumar/.gemini/antigravity-cli/brain/986c8152-cb4b-460a-b37e-91341d120c0a/'
os.makedirs(artifact_dir, exist_ok=True)
shutil.copy(out_png, os.path.join(artifact_dir, 'FigureX_J.png'))
print("Figure X_J copied to brain artifacts folder as FigureX_J.png")
