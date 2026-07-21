import os
import gzip
import numpy as np
import matplotlib.pyplot as plt

# Define paths
input_dir = '../input_data'
output_dir = '../output_figures'
matrix_in = os.path.join(input_dir, 'tss_matrix.gz')

print(f"Reading matrix from {matrix_in}...")

# Load the matrix data
data_rows = []
with gzip.open(matrix_in, 'rt') as f:
    for line in f:
        if not line.startswith('@'):
            parts = line.strip().split('\t')
            data_rows.append([float(x) if x != 'nan' else 0.0 for x in parts[6:]])

data_arr = np.array(data_rows)
num_genes, total_bins = data_arr.shape
num_samples = 9
bins_per_sample = total_bins // num_samples

print(f"Loaded matrix: {num_genes} genes, {total_bins} total bins ({bins_per_sample} bins per sample)")

# Custom color coding matching user outlines:
# RARA = #677fac
# RARG = #cf7579
selected_indices = {
    'RARA': (4, '#677fac'),  # RARA Treated
    'RARG': (8, '#cf7579')   # RARG Treated
}

# Setup plot
fig, ax = plt.subplots(figsize=(2, 2), dpi=300)

# Generate X coordinates (-3.0 to 3.0 kb)
x_coords = np.linspace(-3.0, 3.0, bins_per_sample)

# Plot each sample's profile normalized to fraction of signal
for label, (idx, color) in selected_indices.items():
    start = idx * bins_per_sample
    end = (idx + 1) * bins_per_sample
    
    # Calculate average profile across all genes
    raw_profile = np.mean(data_arr[:, start:end], axis=0)
    
    # Normalize to fraction of total signal
    fraction_profile = raw_profile / np.sum(raw_profile)
    
    ax.plot(x_coords, fraction_profile, color=color, linewidth=1.0, label=label, zorder=3)

# Add vertical grey dashed line at x=0 (TSS) with thin width
ax.axvline(x=0, color='#7f7f7f', linestyle='--', linewidth=0.4, zorder=1)

# Customize axes
ax.set_xlabel('Distance from TSS (kb)', fontsize=7, labelpad=2)
ax.set_ylabel('Fraction of signal', fontsize=7, labelpad=2)
ax.set_xlim(-3.0, 3.0)
ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
ax.set_xticklabels(['-3', '-2', '-1', '0', '1', '2', '3'], fontsize=6)
ax.tick_params(axis='both', which='major', labelsize=6, pad=1, length=2)

# Make it a clean closed box
for spine in ['top', 'right', 'left', 'bottom']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#1a202c')
    ax.spines[spine].set_linewidth(0.6)

# Format y-axis ticks with scientific notation
ax.yaxis.get_offset_text().set_fontsize(5)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Add Legend at the bottom (clean, single row for 2 items)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2, 
          fontsize=6.5, frameon=False, handlelength=1.2, handleheight=0.8, 
          columnspacing=1.5, labelspacing=0.4)

# Save figure
plt.savefig(os.path.join(output_dir, 'FigureX_C.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Figure X.C updated with custom colors #677fac and #cf7579.")
