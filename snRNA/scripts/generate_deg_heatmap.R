library(ComplexHeatmap)
library(circlize)

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_DIR <- "../raw_data"
OUTPUT_DIR <- "../figures"
DAVID_DIR <- INPUT_DIR
P_CUTOFF  <- 0.05
FC_CUTOFF <- 0.5
TOLERANCE <- 0.20  # Tolerance for classifying "no much change" genes

FILES <- list(
  RARA  = file.path(BASE_DIR, "DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv"),
  RARG  = file.path(BASE_DIR, "DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv"),
  Combo = file.path(BASE_DIR, "DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv")
)

# ── Load Data ─────────────────────────────────────────────────────────────────
rara  <- read.csv(FILES$RARA)
rarg  <- read.csv(FILES$RARG)
combo <- read.csv(FILES$Combo)

# Filter out viral vector constructs (pAAV entries)
rara  <- rara[!grepl("^paav", rara$gene, ignore.case = TRUE), ]
rarg  <- rarg[!grepl("^paav", rarg$gene, ignore.case = TRUE), ]
combo <- combo[!grepl("^paav", combo$gene, ignore.case = TRUE), ]

# Identify significant DEGs in each group
rara$sig  <- rara$p_val < P_CUTOFF & abs(rara$avg_log2FC) >= FC_CUTOFF
rarg$sig  <- rarg$p_val < P_CUTOFF & abs(rarg$avg_log2FC) >= FC_CUTOFF
combo$sig <- combo$p_val < P_CUTOFF & abs(combo$avg_log2FC) >= FC_CUTOFF

# Get the union of all genes
all_genes <- unique(c(rara$gene, rarg$gene, combo$gene))

# Create mapping vectors for Fold Change and Significance
fc_rara_map  <- setNames(rara$avg_log2FC, rara$gene)
fc_rarg_map  <- setNames(rarg$avg_log2FC, rarg$gene)
fc_combo_map <- setNames(combo$avg_log2FC, combo$gene)

sig_rara_map  <- setNames(rara$sig, rara$gene)
sig_rarg_map  <- setNames(rarg$sig, rarg$gene)
sig_combo_map <- setNames(combo$sig, combo$gene)

# Build Master DataFrame
master <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
master$FC_RARA  <- ifelse(is.na(fc_rara_map[master$gene]), 0, fc_rara_map[master$gene])
master$FC_RARG  <- ifelse(is.na(fc_rarg_map[master$gene]), 0, fc_rarg_map[master$gene])
master$FC_Combo <- ifelse(is.na(fc_combo_map[master$gene]), 0, fc_combo_map[master$gene])

master$Sig_RARA  <- ifelse(is.na(sig_rara_map[master$gene]), FALSE, sig_rara_map[master$gene])
master$Sig_RARG  <- ifelse(is.na(sig_rarg_map[master$gene]), FALSE, sig_rarg_map[master$gene])
master$Sig_Combo <- ifelse(is.na(sig_combo_map[master$gene]), FALSE, sig_combo_map[master$gene])

# Filter to keep ONLY the common genes (significant in all 3 groups)
master_common <- master[master$Sig_RARA & master$Sig_RARG & master$Sig_Combo, ]

C <- abs(master_common$FC_Combo)
A <- abs(master_common$FC_RARA)
G <- abs(master_common$FC_RARG)

# Classify the 95 common genes into the 4 requested categories
master_common$category <- NA

# 1. Similar (no much change: all three conditions within tolerance of each other)
similar_idx <- (abs(C - A) <= TOLERANCE) & (abs(C - G) <= TOLERANCE) & (abs(A - G) <= TOLERANCE)
master_common$category[similar_idx] <- "Similar"

# 2. Combo higher than individual: C > A and C > G
master_common$category[is.na(master_common$category) & (C > A) & (C > G)] <- "Combo > Indiv"

# 3. High in RARA compared to Combo: A > C and A >= G
master_common$category[is.na(master_common$category) & (A > C) & (A >= G)] <- "RARA > Combo"

# 4. High in RARG compared to Combo: G > C and G > A
master_common$category[is.na(master_common$category) & (G > C) & (G > A)] <- "RARG > Combo"

# Handle any remaining genes that are intermediate (classify by their highest individual FC)
master_common$category[is.na(master_common$category) & (A > C)] <- "RARA > Combo"
master_common$category[is.na(master_common$category) & (G > C)] <- "RARG > Combo"
master_common$category[is.na(master_common$category)]           <- "Combo > Indiv"

print(paste("Total common DEGs plotted in heatmap:", nrow(master_common)))
print(table(master_common$category))

# Prepare matrix
mat <- as.matrix(master_common[, c("FC_RARA", "FC_RARG", "FC_Combo")])
rownames(mat) <- master_common$gene
colnames(mat) <- c("RARA", "RARG", "Combo")

# Order the categories top-to-bottom exactly as requested:
# 1. Similar (Shared in all 3, no changes)
# 2. RARG > Combo (RARG higher)
# 3. RARA > Combo (RARA higher)
# 4. Combo > Indiv (Combo higher)
cat_levels <- c("Similar", "RARG > Combo", "RARA > Combo", "Combo > Indiv")
category   <- factor(master_common$category, levels = cat_levels)

# ── Colors ────────────────────────────────────────────────────────────────────
# Themed Colormap: Slate Grey (#94a3b8) for Downregulated, White (#ffffff) for Neutral, Purple (#9e68e4) for Upregulated
col_fun <- colorRamp2(c(-3, 0, 3), c("#94a3b8", "#ffffff", "#9e68e4"))

# ── Select 2 Genes per Category Programmatically ──────────────────────────────
# We select the top 2 genes with the highest absolute Combo FC in each category
genes_to_mark <- c()
for(cat in cat_levels) {
  cat_df <- master_common[master_common$category == cat, ]
  cat_df_sorted <- cat_df[order(abs(cat_df$FC_Combo), decreasing = TRUE), ]
  top_2_genes <- head(cat_df_sorted$gene, 2)
  genes_to_mark <- c(genes_to_mark, top_2_genes)
}

print("Selected genes to label (2 per category):")
print(genes_to_mark)

row_idx <- match(genes_to_mark, rownames(mat))

# Create Row Annotation with short link width for narrow layout
ha_row <- rowAnnotation(
  link = anno_mark(
    at = row_idx, 
    labels = genes_to_mark, 
    labels_gp = gpar(fontsize = 4.2, fontface = "bold.italic", col = "#1a1a1a"),
    link_width = unit(2.0, "mm"),
    padding = unit(0.5, "mm")
  )
)

# ── Heatmap ───────────────────────────────────────────────────────────────────
ht <- Heatmap(
  mat,
  name = "log2(FC)",
  col = col_fun,
  row_split = category,
  cluster_columns = FALSE,
  column_order = c("RARA", "RARG", "Combo"),
  show_row_names = FALSE,
  right_annotation = ha_row,
  
  # Thin outer border only (no lines between individual cells)
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.25),
  rect_gp = gpar(col = NA), # remove cell borders
  width = unit(7.0, "mm"), # 2.33 mm per column (total 7 mm for the 3 columns)
  
  # Row titles (categories) - rotated 90 degrees and styled very small to fit narrow layout
  row_title_gp = gpar(fontsize = 5.0, fontface = "bold", col = "#1a1a1a"),
  row_title_rot = 90,
  row_title_side = "left",
  
  # Column names styled tiny and rotated 90 degrees
  column_names_gp = gpar(fontsize = 5.0, fontface = "bold", col = "#1a1a1a"),
  column_names_rot = 90,
  column_names_centered = TRUE,
  
  # Clustering params
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_dend_width = unit(3.0, "mm"),  # Compact dendrogram tree
  row_dend_gp = gpar(lwd = 0.25),
  
  # DO NOT cluster the row slices (removes the global dotted dendrogram tree connecting categories)
  cluster_row_slices = FALSE,
  
  # Legend parameters (horizontal legend to save width)
  heatmap_legend_param = list(
    title = "log2(FC)",
    direction = "horizontal",
    title_gp = gpar(fontsize = 5.0, fontface = "bold"),
    labels_gp = gpar(fontsize = 4.5),
    grid_height = unit(2, "mm"),
    legend_width = unit(1.8, "cm")
  )
)

# Save Plot (exactly 1.2" width x 3.5" height, high-DPI publication ready)
out_path <- file.path(OUTPUT_DIR, "Figure_X_G_Common_vertical_heatmap.png")
png(out_path, width = 1.2, height = 3.5, units = "in", res = 600)
# Draw the heatmap with legend at the bottom to conserve width
draw(
  ht, 
  heatmap_legend_side = "bottom",
  padding = unit(c(1, 1, 1.2, 1), "mm") # slight extra padding at bottom for rotated column names
)
dev.off()

print(paste("Heatmap saved successfully to:", out_path))
