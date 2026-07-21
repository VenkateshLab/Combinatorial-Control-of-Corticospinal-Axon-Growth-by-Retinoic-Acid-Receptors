import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_DIR = "../raw_data"
OUTPUT_DIR = "../figures"
P_CUTOFF  = 0.05
FC_CUTOFF = 0.5

FILES = {
    "RARA":  "DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv",
    "RARG":  "DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv",
    "Combo": "DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv",
}

# Colors
SAMPLE_COLORS = {
    "RARA":  "#E05C5C",   # coral-red
    "RARG":  "#5B9BD5",   # steel-blue
    "Combo": "#70AD47",   # leaf-green
}

COLOR_UP   = "#9e68e4"   # purple
COLOR_DOWN = "#94a3b8"   # slate grey

# ── Load gene sets ────────────────────────────────────────────────────────────
up_sets, down_sets = {}, {}

for label, fname in FILES.items():
    df  = pd.read_csv(os.path.join(INPUT_DIR, fname))
    sig = df[df["p_val"] < P_CUTOFF]
    up_sets[label]   = set(sig[sig["avg_log2FC"] >=  FC_CUTOFF]["gene"])
    down_sets[label] = set(sig[sig["avg_log2FC"] <= -FC_CUTOFF]["gene"])

set_names = list(FILES.keys())
n = len(set_names)

# ── Function to compute and sort intersections ────────────────────────────────
def get_sorted_intersections(sets):
    combos = []
    for r in range(1, n + 1):
        for subset in itertools.combinations(set_names, r):
            present = frozenset(subset)
            absent  = frozenset(set_names) - present
            
            genes = set.intersection(*[sets[s] for s in present])
            if absent:
                genes -= set.union(*[sets[s] for s in absent])
                
            cnt = len(genes)
            if cnt > 0:
                combos.append({
                    "present": present,
                    "count": cnt
                })
    combos.sort(key=lambda x: x["count"], reverse=True)
    return combos

up_combos = get_sorted_intersections(up_sets)
down_combos = get_sorted_intersections(down_sets)

# ── Create Figure (2.0 x 1.5 inches) ──────────────────────────────────────────
fig = plt.figure(figsize=(2.0, 1.5), dpi=300)

# Axes Coordinates
# [left, bottom, width, height]
ax_bar_up   = fig.add_axes([0.16, 0.46, 0.31, 0.38])
ax_mat_up   = fig.add_axes([0.16, 0.08, 0.31, 0.31])

ax_bar_down = fig.add_axes([0.65, 0.46, 0.31, 0.38])
ax_mat_down = fig.add_axes([0.65, 0.08, 0.31, 0.31])

# ── Helper to Plot Sub-UpSet ──────────────────────────────────────────────────
def plot_upset_sub(ax_bar, ax_mat, combos, color, title):
    n_bars = len(combos)
    xs = np.arange(n_bars)
    counts = [c["count"] for c in combos]
    
    # 1. Bar Chart
    bars = ax_bar.bar(xs, counts, color=color, width=0.58, edgecolor='none')
    for xi, cnt in enumerate(counts):
        ax_bar.text(xi, cnt + 2, str(cnt),
                    ha='center', va='bottom', fontsize=3.5, fontweight='bold', color='#000000')
                    
    ax_bar.set_xlim(-0.6, n_bars - 0.4)
    ax_bar.set_ylim(0, 120)  # Shared y-limit for direct comparison
    ax_bar.set_xticks([])
    ax_bar.set_yticks([0, 50, 100])
    ax_bar.set_ylabel("DEGs", fontsize=4.5, labelpad=1)
    ax_bar.tick_params(axis='y', labelsize=4.0, pad=1)
    
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.spines['bottom'].set_visible(True)
    ax_bar.spines['bottom'].set_color('#000000')
    ax_bar.spines['bottom'].set_linewidth(0.4)
    ax_bar.spines['left'].set_color('#000000')
    ax_bar.spines['left'].set_linewidth(0.4)
    ax_bar.yaxis.grid(True, linestyle=':', alpha=0.15, linewidth=0.25)
    ax_bar.set_axisbelow(True)
    ax_bar.set_title(title, fontsize=5.5, fontweight='bold', pad=3)

    # 2. Dot Matrix
    ax_mat.set_xlim(-0.6, n_bars - 0.4)
    ax_mat.set_ylim(-0.6, n - 0.4)
    ax_mat.set_xticks([])
    ax_mat.set_yticks([])
    ax_mat.spines[:].set_visible(False)
    
    for xi, c in enumerate(combos):
        present = c["present"]
        active_ys = []
        for yi, sname in enumerate(set_names):
            on = sname in present
            dot_color = SAMPLE_COLORS[sname] if on else '#dde0e5'
            ax_mat.scatter(xi, yi, s=12, zorder=4, color=dot_color, linewidths=0)
            if on:
                active_ys.append(yi)
        if len(active_ys) > 1:
            ax_mat.plot([xi, xi], [min(active_ys), max(active_ys)],
                        color='#444444', lw=0.8, zorder=3)

# ── Draw UP and DOWN Plots ───────────────────────────────────────────────────
plot_upset_sub(ax_bar_up, ax_mat_up, up_combos, COLOR_UP, "Upregulated")
plot_upset_sub(ax_bar_down, ax_mat_down, down_combos, COLOR_DOWN, "Downregulated")

# ── Draw Row Labels (Set Names) in the Gaps ───────────────────────────────────
mat_bottom = 0.08
mat_height = 0.31
ylim_range = (n - 0.4) - (-0.6)  # 3.0

for yi, sname in enumerate(set_names):
    fig_y = mat_bottom + (yi + 0.6) / ylim_range * mat_height
    # Labels for UP Matrix (placed on the left of it)
    fig.text(0.14, fig_y, sname,
             ha='right', va='center',
             fontsize=4.5, fontweight='bold',
             color=SAMPLE_COLORS[sname])
             
    # Labels for DOWN Matrix (placed on the left of it)
    fig.text(0.63, fig_y, sname,
             ha='right', va='center',
             fontsize=4.5, fontweight='bold',
             color=SAMPLE_COLORS[sname])

# Save Plot
output_plot_path = os.path.join(OUTPUT_DIR, "Figure_X_E_DEG_upset.png")
plt.savefig(output_plot_path, dpi=600, bbox_inches='tight', pad_inches=0.02)
plt.close()
print(f"Side-by-side UpSet plot saved to {output_plot_path}")
