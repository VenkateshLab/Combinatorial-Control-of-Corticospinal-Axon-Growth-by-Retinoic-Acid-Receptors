import os
import shutil
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

# Define paths
logo_dir = 'manuscript_images/output/logos_png/rara'
output_dir = '../output_figures'
os.makedirs(output_dir, exist_ok=True)

# Exact values and ordering requested by the user
rows_data = [
    {
        'name': 'Rarg',
        'id': 'MA0860.1',
        'rara_text': '16.1%',
        'rarg_text': '96.6%'
    },
    {
        'name': 'Rxra',
        'id': 'MA0512.2',
        'rara_text': '27.6%',
        'rarg_text': '13.7%'
    },
    {
        'name': 'RARA',
        'id': 'MA0730.1',
        'rara_text': '25.3%',
        'rarg_text': '14.3%'
    },
    {
        'name': 'RARA::RXRA',
        'id': 'MA0159.1',
        'rara_text': '93.1%',
        'rarg_text': '7.2%'
    }
]

# Render Figure (2.0 x 2.5 inches, 500 DPI)
fig = plt.figure(figsize=(2.0, 2.5), dpi=500)
# gridspec for precise layout alignment
gs = gridspec.GridSpec(5, 3, width_ratios=[1.25, 0.9, 0.9], wspace=0.08, hspace=0.1,
                       left=0.05, right=0.95, top=0.92, bottom=0.08)

# Table Headers
headers = ['Motif', 'RARA %', 'RARG %']
for col_idx, h in enumerate(headers):
    ax = fig.add_subplot(gs[0, col_idx])
    ax.axis('off')
    ax.text(0.5, 0.5, h, weight='bold', size=5, ha='center', va='center', color='#1f2937')

# Populate Rows
for row_idx, item in enumerate(rows_data):
    # Col 0: Logo & Name label
    ax_logo = fig.add_subplot(gs[row_idx+1, 0])
    ax_logo.axis('off')
    
    logo_path = os.path.join(logo_dir, f"jaspar_{item['id']}.png")
    if os.path.exists(logo_path):
        img = mpimg.imread(logo_path)
        ax_logo.imshow(img, aspect='equal')
    
    # Place text label below logo
    ax_logo.text(0.0, -0.22, f"{item['name']} ({item['id']})", 
                 size=3.5, ha='left', va='top', transform=ax_logo.transAxes, weight='semibold', color='#374151')
    
    # Col 1: RARA Stats
    ax_rara = fig.add_subplot(gs[row_idx+1, 1])
    ax_rara.axis('off')
    ax_rara.text(0.5, 0.5, item['rara_text'], size=5, ha='center', va='center', color='#111827')
    
    # Col 2: RARG Stats
    ax_rarg = fig.add_subplot(gs[row_idx+1, 2])
    ax_rarg.axis('off')
    ax_rarg.text(0.5, 0.5, item['rarg_text'], size=5, ha='center', va='center', color='#111827')

# Draw professional booktabs divider lines
line_top = plt.Line2D([0.05, 0.95], [0.93, 0.93], color='black', linewidth=0.8)
line_mid = plt.Line2D([0.05, 0.95], [0.77, 0.77], color='black', linewidth=0.4)
line_bot = plt.Line2D([0.05, 0.95], [0.08, 0.08], color='black', linewidth=0.8)
fig.add_artist(line_top)
fig.add_artist(line_mid)
fig.add_artist(line_bot)

# Save to output folder
out_png = os.path.join(output_dir, 'FigureX_H.png')
plt.savefig(out_png, dpi=500, bbox_inches='tight', transparent=True)
plt.close()
print(f"Generated Figure X.H table: {out_png}")

# Copy to brain artifacts folder
artifact_dir = '/home/manojkumar/.gemini/antigravity-cli/brain/ecd96054-3de7-4ed8-bc31-e4fc2fa45bc5/'
os.makedirs(artifact_dir, exist_ok=True)
shutil.copy(out_png, os.path.join(artifact_dir, 'FigureX_H.png'))
print("Figure X.H copied to brain artifacts folder as FigureX_H.png")
