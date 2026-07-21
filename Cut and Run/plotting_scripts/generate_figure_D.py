import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Define paths
input_dir = '../input_data'
output_dir = '../output_figures'

rara_path = os.path.join(input_dir, 'Peaks_Annotated_RARA_vs_Input_p1e4.csv')
rarg_path = os.path.join(input_dir, 'Peaks_Annotated_RARG_vs_Input_p1e4.csv')

# Load files
df_rara = pd.read_csv(rara_path)
df_rarg = pd.read_csv(rarg_path)

rara_fe = df_rara['fold_enrichment'].dropna().values
rarg_fe = df_rarg['fold_enrichment'].dropna().values

# Perform Mann-Whitney U test
stat, p_val = stats.mannwhitneyu(rara_fe, rarg_fe, alternative='two-sided')
print(f"Mann-Whitney U p-value: {p_val:.4e}")

# Transform data to log2 scale
rara_fe_log = np.log2(rara_fe)
rarg_fe_log = np.log2(rarg_fe)

# Setup plot
fig, ax = plt.subplots(figsize=(2, 2), dpi=300)

# Colors: exact hex codes corrected by the user
colors_face = {
    'RARA': '#899cbf',
    'RARG': '#da9497'
}
colors_edge = {
    'RARA': '#677fac',
    'RARG': '#cf7579'
}

# 1. Plot Violins (width=0.7) on log2 values
vp = ax.violinplot([rara_fe_log, rarg_fe_log], positions=[0, 1], showmeans=False, showmedians=False, showextrema=False, widths=0.7)

# Customize violin bodies (set exact solid colors)
for i, body in enumerate(vp['bodies']):
    color_key = 'RARA' if i == 0 else 'RARG'
    body.set_facecolor(colors_face[color_key])
    body.set_edgecolor(colors_edge[color_key])
    body.set_linewidth(1.0)
    body.set_alpha(1.0)  # Solid rendering using the exact user-picked blended colors

# 2. Plot Boxplots inside the violins (width=0.14) on log2 values
bp = ax.boxplot([rara_fe_log, rarg_fe_log], positions=[0, 1], widths=0.14, 
                showfliers=False, patch_artist=True, zorder=3)

# Style boxplot elements
for i, box in enumerate(bp['boxes']):
    box.set_facecolor('#ffffff')  # White fill for box
    box.set_edgecolor('#1a202c')  # Dark charcoal border
    box.set_linewidth(0.8)
    box.set_alpha(0.8)

for whisker in bp['whiskers']:
    whisker.set_color('#1a202c')
    whisker.set_linewidth(0.8)

for cap in bp['caps']:
    cap.set_color('#1a202c')
    cap.set_linewidth(0.8)

for median in bp['medians']:
    median.set_color('#1a202c')
    median.set_linewidth(1.0)

# Customize axes
ax.set_ylabel('Fold enrichment', fontsize=7, labelpad=2)
ax.set_xticks([0, 1])
ax.set_xticklabels(['RARA', 'RARG'], fontsize=7)

# Define ticks in log2 space (starting from 0 up to 7)
y_ticks_log = [0, 1, 2, 3, 4, 5, 6, 7]
y_ticks_labels = ['1', '2', '4', '8', '16', '32', '64', '128']
ax.set_yticks(y_ticks_log)
ax.set_yticklabels(y_ticks_labels, fontsize=6)
ax.tick_params(axis='both', which='major', labelsize=6, pad=1, length=2)

# Set Y-limit from 0 to 9.2 in log2 space (gives headroom and prevents bottom clipping)
ax.set_ylim(0, 9.2)

# Make it a clean closed box
for spine in ['top', 'right', 'left', 'bottom']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#1a202c')
    ax.spines[spine].set_linewidth(0.6)

# 3. Add Significance Bracket in log2 space
# Bracket height
y_line = 7.7
y_tick = 7.2
# Draw horizontal line and vertical ticks
ax.plot([0, 0, 1, 1], [y_tick, y_line, y_line, y_tick], color='#4a5568', linewidth=0.6)

# Add p-value text above the bracket
p_text = f"***\np = {p_val:.2e}".replace('e-09', 'x10$^{-9}$')
ax.text(0.5, y_line + 0.1, p_text, ha='center', va='bottom', fontsize=5.5, weight='bold', color='#2d3748', linespacing=0.9)

# Save figure
plt.savefig(os.path.join(output_dir, 'FigureX_D.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Figure X.D saved with corrected solid custom colors.")
