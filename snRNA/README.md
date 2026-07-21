# Differential Expression and GO Enrichment Analysis Assets

This repository contains the raw data, execution scripts, and high-resolution figures for the differential gene expression (DEG) and Gene Ontology (GO) pathway enrichment analyses of **RARA**, **RARG**, and **Combo** versus **GFP** control.

## Directory Structure

github/snRNA/
├── raw_data/                 # Raw input tables and intermediate CSVs
├── scripts/                  # Python and R scripts to reproduce analysis
├── figures/                  # High-resolution (600 DPI) publication-ready plots
├── CellRanger.sh             # CellRanger alignment and counting docker script
└── README.md                 # Project summary and description
```

---

## Figures Description (Figure X Panels C to I)
*Note: Panels A and B are schematics not included in this directory.*

### [Figure X.C: DEG Bi-directional Barplot](figures/Figure_X_C_DEG_bidirectional_barplot.png)
* **Description**: Shows the total number of significant DEGs in RARA, RARG, and Combo grouped by fold-change ranges ($0.5\text{ to }1.0$, $1.0\text{ to }2.0$, $2.0\text{ to }3.0$, and $>3.0$).
* **Script**: `scripts/deg_bidirectional_plot.py`

### [Figure X.D: DEG Proportions Plot](figures/Figure_X_D_DEG_proportions.png)
* **Description**: A stacked percentage bar plot comparing the relative fractions of upregulated and downregulated DEGs.
* **Script**: `scripts/generate_deg_proportions_plot.py`

### [Figure X.E: DEG UpSet Intersections](figures/Figure_X_E_DEG_upset.png)
* **Description**: An UpSet plot displaying the intersections of DEGs between the RARA, RARG, and Combo groups.
* **Script**: `scripts/generate_combined_upset_plot.py`

### [Figure X.F: Common Genes GO Butterfly Plot](figures/Figure_X_F_Common_GO_butterfly.png)
* **Description**: A bilateral butterfly barplot for the 95 common genes (significant in all 3 groups), showing Fold Enrichment for key biological processes (UP on the right in purple, DOWN on the left in blue/grey).
* **Script**: `scripts/plot_all_3_go_butterfly.py`

### [Figure X.G: Common Genes Vertical Heatmap](figures/Figure_X_G_Common_vertical_heatmap.png)
* **Description**: A vertical split heatmap of the 95 common genes divided into four categories: *Similar*, *RARG > Combo*, *RARA > Combo*, and *Combo > Indiv*. Features row dendrograms and thin outer borders.
* **Script**: `scripts/generate_deg_heatmap.R`

### [Figure X.H: Unique Genes GO Butterfly Plot](figures/Figure_X_H_Unique_GO_butterfly.png)
* **Description**: A bilateral butterfly barplot for the unique (condition-specific) GO terms, grouped and colored by condition (RARA = light, RARG = medium, Combo = dark).
* **Script**: `scripts/plot_unique_go_butterfly.py`

### [Figure X.I: Unique Genes Horizontal Heatmap](figures/Figure_X_I_Unique_horizontal_heatmap.png)
* **Description**: A horizontal split heatmap of the 472 unique DEGs. Non-significant fold changes are masked to exactly 0 (White). Labels exactly 10 biologically relevant genes (5 UP, 5 DOWN related to regeneration/remodeling/metabolism) per category.
* **Script**: `scripts/generate_deg_heatmap_horizontal.R`

---

## Inputs and Raw Data

* **DEG CSVs**:
  * `DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv`
  * `DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv`
  * `DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv`
* **DAVID Reports**: Raw CSV tables exported from DAVID Functional Annotation Tool.
* **Intermediate Matrices**:
  * `common_genes_fold_changes.csv`: Matrix of raw fold change values for the 95 common genes.
  * `unique_genes_masked_fold_changes.csv`: Matrix of masked fold change values for the 472 unique genes.
