import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_DIR = "../raw_data"
OUTPUT_DIR = "../figures"
P_CUTOFF = 0.05
FILES = {
    "RARA":  "DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv",
    "RARG":  "DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv",
    "Combo": "DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv",
}
range_keys = ["0.5 to 1.0", "1.0 to 2.0", "2.0 to 3.0", "> 3.0"]

up_colors   = ['#ebdffd', '#cca8fa', '#9e68e4', '#7030c0']
down_colors = ['#eceef1', '#cbd5e1', '#94a3b8', '#475569']

# ── Count DEGs from CSV ───────────────────────────────────────────────────────
up_data   = {r: [] for r in range_keys}
down_data = {r: [] for r in range_keys}

for fname in FILES.values():
    df  = pd.read_csv(os.path.join(INPUT_DIR, fname))
    sig = df[df["p_val"] < P_CUTOFF]
    u   = sig[sig["avg_log2FC"] >=  0.5]["avg_log2FC"]
    d   = sig[sig["avg_log2FC"] <= -0.5]["avg_log2FC"]

    up_data["0.5 to 1.0"].append(int(((u >= 0.5) & (u < 1.0)).sum()))
    up_data["1.0 to 2.0"].append(int(((u >= 1.0) & (u < 2.0)).sum()))
    up_data["2.0 to 3.0"].append(int(((u >= 2.0) & (u < 3.0)).sum()))
    up_data["> 3.0"].append(int((u >= 3.0).sum()))

    down_data["0.5 to 1.0"].append(int(((d <= -0.5) & (d > -1.0)).sum()))
    down_data["1.0 to 2.0"].append(int(((d <= -1.0) & (d > -2.0)).sum()))
    down_data["2.0 to 3.0"].append(int(((d <= -2.0) & (d > -3.0)).sum()))
    down_data["> 3.0"].append(int((d <= -3.0).sum()))

# ── Layout — identical y positions as original ────────────────────────────────
fig, ax = plt.subplots(figsize=(3.3, 3.0), dpi=300)

# bars: order top→bottom is >3, 2-3, 1-2, 0.5-1  (same as original)
rara_y  = np.array([4.7, 4.4, 4.1, 3.8])
rarg_y  = np.array([3.3, 3.0, 2.7, 2.4])
combo_y = np.array([1.9, 1.6, 1.3, 1.0])

bar_height = 0.25
gap        = 55    # centre gap — room for RARA/RARG/Combo labels

# ── Annotation helper ─────────────────────────────────────────────────────────
def annotate_bar(ax, y, val, direction):
    if val <= 0:
        return
    if direction == "up":
        if val >= 50:
            ax.text(gap + val / 2.0, y, str(val),
                    ha='center', va='center', fontsize=4.5,
                    fontweight='bold', color='#000000')
        else:
            ax.text(gap + val + 4, y, str(val),
                    ha='left', va='center', fontsize=4.5, color='#000000')
    else:
        if val >= 50:
            ax.text(-(gap + val / 2.0), y, str(val),
                    ha='center', va='center', fontsize=4.5,
                    fontweight='bold', color='#000000')
        else:
            ax.text(-(gap + val + 4), y, str(val),
                    ha='right', va='center', fontsize=4.5, color='#000000')

# ── Draw bars ─────────────────────────────────────────────────────────────────
group_ys = [rara_y, rarg_y, combo_y]

for r_idx, range_name in enumerate(range_keys):
    for g_idx, gy in enumerate(group_ys):
        y   = gy[r_idx]
        u   = up_data[range_name][g_idx]
        d   = down_data[range_name][g_idx]

        ax.barh(y,  u, left= gap, height=bar_height, color=up_colors[r_idx],   edgecolor='none')
        ax.barh(y, -d, left=-gap, height=bar_height, color=down_colors[r_idx], edgecolor='none')
        annotate_bar(ax, y, u, "up")
        annotate_bar(ax, y, d, "down")

# ── Centre labels — same positions as original ────────────────────────────────
for text, y in [("RARA", 4.25), ("RARG", 2.85), ("Combo", 1.45)]:
    ax.text(0, y, text, ha='center', va='center',
            fontsize=5.5, fontweight='bold', color='#000000',
            bbox=dict(facecolor='white', edgecolor='#000000',
                      boxstyle='round,pad=0.2', linewidth=0.5))

# ── Separator lines — same positions as original ──────────────────────────────
ax.axhline(y=3.55, color='#e5e5e5', linestyle='-', linewidth=0.4, zorder=0)
ax.axhline(y=2.15, color='#e5e5e5', linestyle='-', linewidth=0.4, zorder=0)

# ── Up / Down labels — outside the top spine (axes-fraction coords) ───────────
ax.text(0.25, 1.02, "◀ Down-regulated",
        ha='center', va='bottom', fontsize=5.5,
        color='#475569', fontweight='bold',
        transform=ax.transAxes, clip_on=False)
ax.text(0.75, 1.02, "Up-regulated ▶",
        ha='center', va='bottom', fontsize=5.5,
        color='#5a1fa0', fontweight='bold',
        transform=ax.transAxes, clip_on=False)

# ── Central dashed lines ──────────────────────────────────────────────────────
ax.axvline(x= gap, color='#000000', linestyle='--', linewidth=0.4)
ax.axvline(x=-gap, color='#000000', linestyle='--', linewidth=0.4)

# ── Y axis ────────────────────────────────────────────────────────────────────
ax.set_ylim(0.7, 5.0)
ax.set_yticks([])

# ── X axis ────────────────────────────────────────────────────────────────────
max_val = 160
step    = 50
# ticks at 0, 50, 100, 150 on each side; "0" shown only at ±gap
pos_ticks = np.arange(0, max_val + step, step)          # [0, 50, 100, 150]
x_ticks   = np.concatenate([-(pos_ticks[::-1] + gap), pos_ticks + gap])
x_labels  = [str(t) for t in np.concatenate([pos_ticks[::-1], pos_ticks])]
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels, fontsize=5.0, color='#000000')
ax.set_xlim(-(gap + max_val + 15), gap + max_val + 15)

ax.set_xlabel('Number of DEGs', fontsize=6.5, fontweight='bold',
              color='#000000', labelpad=2)

ax.grid(axis='x', linestyle=':', color='#000000', alpha=0.15, linewidth=0.3)

for spine in ['top', 'bottom', 'left', 'right']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_color('#000000')
    ax.spines[spine].set_linewidth(0.5)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_labels = [
    "Up [0.5-1]", "Up [1-2]", "Up [2-3]", "Up (>3)",
    "Down [-0.5 to -1]", "Down [-1 to -2]", "Down [-2 to -3]", "Down (<-3)"
]
handles = [mpatches.Patch(color=c) for c in up_colors + down_colors]
ax.legend(handles, legend_labels,
          loc='upper center', bbox_to_anchor=(0.5, -0.22),
          ncol=4, frameon=True, fontsize=4.5,
          handlelength=1.2, handletextpad=0.3,
          columnspacing=0.6, borderpad=0.3, edgecolor='#000000')

plt.tight_layout()

out = os.path.join(OUTPUT_DIR, "Figure_X_C_DEG_bidirectional_barplot.png")
plt.savefig(out, dpi=600, bbox_inches='tight')
plt.close()
print(f"Saved → {out}")
