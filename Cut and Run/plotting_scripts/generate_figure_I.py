import os
import math
import shutil
import matplotlib
import pyBigWig

# Override global Matplotlib line widths to make plot borders/lines half-thickness (0.5 pt)
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5

# Override pyGenomeTracks default figure margins
import pygenometracks.tracksClass
pygenometracks.tracksClass.DEFAULT_MARGINS['left'] = 0.08    # Adequate space on left to prevent label cropping
pygenometracks.tracksClass.DEFAULT_MARGINS['right'] = 0.99   # Let the plot occupy the full right edge of the screen
pygenometracks.tracksClass.DEFAULT_MARGINS['bottom'] = 0.02
pygenometracks.tracksClass.DEFAULT_MARGINS['top'] = 0.98

from pygenometracks.tracksClass import PlotTracks

# Define paths
rara_bw = '../input_data/RARA_Treated_rep1.bw'
rara_input_bw = '../input_data/RARA_Input_rep1.bw'
rarg_bw = '../input_data/RARG_Treated_rep1.bw'
rarg_input_bw = '../input_data/RARG_Input_rep1.bw'

gtf_path = '../input_data/genes_subset.gtf'
rara_motifs_bed = '../input_data/rara_motifs.bed9'
rarg_motifs_bed = '../input_data/rarg_motifs.bed9'

output_dir = '../output_figures'
os.makedirs(output_dir, exist_ok=True)

# 1. Define zoomed-in regions (2.5 Kb windows centered on peak summits or TSS)
zoomed_regions = [
    # (Gene, Coordinates)
    ('Rarg', 'chr15:102238250-102240750'),    # Co-bound summit ~102,239,500
    ('Nr2f1', 'chr13:78188500-78191000'),      # RARG peak summit ~78,189,743
    ('Nr2f6', 'chr8:71374650-71377150'),       # RARG peak summit ~71,375,902
    ('Lhx2', 'chr2:38338000-38340500'),        # Negative control centered on TSS ~38,339,281
    ('Npas3', 'chr12:53251350-53253850'),      # RARG peak summit ~53,252,600
    ('Nemf', 'chr12:69360000-69362500'),       # Co-bound peak summit ~69,361,210
    ('Rps29', 'chr12:69158300-69160800')       # Co-bound peak summit ~69,159,550
]

# 2. pyGenomeTracks Template (Omit title keys to remove labels on the left)
track_template = """
[RARA_Input]
file = {rara_input_bw}
height = 1.0
color = black
alpha = 1.0
min_value = 0
max_value = {max_limit}
show_data_range = True
summary_method = mean
type = fill
overlay_previous = no

[RARA_Treated]
file = {rara_bw}
color = #899cbf
alpha = 0.8
min_value = 0
max_value = {max_limit}
summary_method = mean
type = fill
overlay_previous = share-y

[spacer]

[RARG_Input]
file = {rarg_input_bw}
height = 1.0
color = black
alpha = 1.0
min_value = 0
max_value = {max_limit}
show_data_range = True
summary_method = mean
type = fill
overlay_previous = no

[RARG_Treated]
file = {rarg_bw}
color = #da9497
alpha = 0.8
min_value = 0
max_value = {max_limit}
summary_method = mean
type = fill
overlay_previous = share-y

[spacer]

[genes]
file = {gtf_path}
height = 0.5
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
labels = true
fontsize = 6
style = UCSC
gene_rows = 1
display = stacked
color = black

[spacer]

[RARA_Motifs]
file = {rara_motifs}
height = 0.25
file_type = bed
color = bed_rgb
border_color = none
labels = false
display = collapsed

[RARG_Motifs]
file = {rarg_motifs}
height = 0.25
file_type = bed
color = bed_rgb
border_color = none
labels = false
display = collapsed
"""

# Open BigWig files
bw_rara_t = pyBigWig.open(rara_bw)
bw_rara_i = pyBigWig.open(rara_input_bw)
bw_rarg_t = pyBigWig.open(rarg_bw)
bw_rarg_i = pyBigWig.open(rarg_input_bw)

# Generate zoomed plots
for name, region in zoomed_regions:
    chrom, coords = region.split(':')
    start, end = [int(x) for x in coords.split('-')]
    
    # Calculate max values in zoomed region
    def get_max_val(bw_t, bw_i, chrom, start, end):
        v_t = bw_t.stats(chrom, start, end, type='max')[0]
        v_i = bw_i.stats(chrom, start, end, type='max')[0]
        v_t = v_t if v_t is not None else 0.0
        v_i = v_i if v_i is not None else 0.0
        return max(v_t, v_i)
        
    max_rara = get_max_val(bw_rara_t, bw_rara_i, chrom, start, end)
    max_rarg = get_max_val(bw_rarg_t, bw_rarg_i, chrom, start, end)
    
    # Shared max limit
    overall_max = max(max_rara, max_rarg)
    max_limit = max(1.0, math.ceil(overall_max * 1.1))
    
    # Format INI file
    ini_file = f'../output_figures/tracks_{name}_zoomed.ini'
    with open(ini_file, 'w') as f:
        f.write(track_template.format(
            rara_bw=rara_bw,
            rara_input_bw=rara_input_bw,
            rarg_bw=rarg_bw,
            rarg_input_bw=rarg_input_bw,
            gtf_path=gtf_path,
            rara_motifs=rara_motifs_bed,
            rarg_motifs=rarg_motifs_bed,
            max_limit=max_limit
        ))
        
    # Render using the Python API of pyGenomeTracks to honor overrides
    out_png = f'../output_figures/FigureX_I_{name}.png'
    print(f"Generating label-free zoomed peak tracks for {name} ({region}) at 2.5x1.25 inches, 500 DPI...")
    
    # Initialize PlotTracks with:
    # Width: 2.5 inches = 6.35 cm
    # Height: 1.25 inches = 3.175 cm
    # FontSize: 5
    # DPI: 500
    # Track Label Width fraction: 0.08
    fig = PlotTracks(ini_file, 
                     fig_width=6.35, 
                     fig_height=3.175, 
                     fontsize=5, 
                     dpi=500,
                     track_label_width=0.08,
                     plot_regions=[(chrom, start, end)])
                     
    fig.plot(out_png, chrom, start, end)
    print(f"Successfully generated: {out_png}")

# Close BigWig files
bw_rara_t.close()
bw_rara_i.close()
bw_rarg_t.close()
bw_rarg_i.close()

# 3. Copy zoomed files to brain artifacts directory (disabled for self-contained package)
# artifact_dir = '/home/manojkumar/.gemini/antigravity-cli/brain/ecd96054-3de7-4ed8-bc31-e4fc2fa45bc5/zoomed_peaks'
# os.makedirs(artifact_dir, exist_ok=True)
# 
# for name, _ in zoomed_regions:
#     src_png = f'../output_figures/FigureX_I_{name}.png'
#     dst_png = os.path.join(artifact_dir, f'FigureX_I_{name}.png')
#     shutil.copy(src_png, dst_png)
    
# 4. Generate updated legend
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

fig_leg, ax_leg = plt.subplots(figsize=(8, 0.6))

legend_elements = [
    Patch(facecolor='#899cbf', label='RARA CUT&RUN'),
    Patch(facecolor='#da9497', label='RARG CUT&RUN'),
    Patch(facecolor='black', label='Input (rep1)'),
    Patch(facecolor='#899cbf', label='RARA JASPAR Motif'),
    Patch(facecolor='#da9497', label='RARG JASPAR Motif')
]
ax_leg.legend(handles=legend_elements, loc='center', ncol=3, frameon=False, fontsize=9.5)
ax_leg.axis('off')

legend_png = '../output_figures/FigureX_I_legend.png'
plt.savefig(legend_png, dpi=500, bbox_inches='tight', transparent=True)
plt.close()
print(f"Generated updated legend: {legend_png}")

# Copy legend to brain artifacts directory (disabled for self-contained package)
# shutil.copy(legend_png, os.path.join(artifact_dir, 'FigureX_I_legend.png'))
# print("Updated legend copied to brain artifacts directory.")

