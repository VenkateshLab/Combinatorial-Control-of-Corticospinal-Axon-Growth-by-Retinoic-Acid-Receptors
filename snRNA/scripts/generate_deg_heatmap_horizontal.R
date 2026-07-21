library(ComplexHeatmap)
library(circlize)

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_DIR <- "../raw_data"
OUTPUT_DIR <- "../figures"
DAVID_DIR <- INPUT_DIR
P_CUTOFF  <- 0.05
FC_CUTOFF <- 0.5

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

# Classify into the 3 unique categories (significant in ONLY one condition) with counts in labels
master$category <- NA
master$category[master$Sig_RARA & !master$Sig_RARG & !master$Sig_Combo] <- "RARA Unique (n=136)"
master$category[master$Sig_RARG & !master$Sig_RARA & !master$Sig_Combo] <- "RARG Unique (n=163)"
master$category[master$Sig_Combo & !master$Sig_RARA & !master$Sig_RARG] <- "Combo Unique (n=173)"

# Filter to keep only these 3 unique categories
master_unique <- master[!is.na(master$category), ]

# Apply Significance Mask: set fold change to 0 if the gene is not significant in that condition
master_unique$FC_RARA  <- ifelse(master_unique$Sig_RARA, master_unique$FC_RARA, 0)
master_unique$FC_RARG  <- ifelse(master_unique$Sig_RARG, master_unique$FC_RARG, 0)
master_unique$FC_Combo <- ifelse(master_unique$Sig_Combo, master_unique$FC_Combo, 0)

print(paste("Total unique DEGs plotted in horizontal masked heatmap:", nrow(master_unique)))
print(table(master_unique$category))

# Save the master CSV file containing the masked unique genes fold change matrix
master_csv_path <- file.path(INPUT_DIR, "unique_genes_masked_fold_changes.csv")
write.csv(
  master_unique[, c("gene", "FC_RARA", "FC_RARG", "FC_Combo", "category")],
  master_csv_path,
  row.names = FALSE
)
print(paste("Saved master CSV file to:", master_csv_path))

# Prepare matrix (Rows = genes, Columns = conditions)
mat <- as.matrix(master_unique[, c("FC_RARA", "FC_RARG", "FC_Combo")])
rownames(mat) <- master_unique$gene
colnames(mat) <- c("RARA", "RARG", "Combo")

# Transpose matrix for horizontal heatmap (Rows = conditions, Columns = genes)
mat_t <- t(mat)
colnames(mat_t) <- master_unique$gene
rownames(mat_t) <- c("RARA", "RARG", "Combo")

# Order the column categories (slices): RARA Unique, RARG Unique, Combo Unique
cat_levels <- c("RARA Unique (n=136)", "RARG Unique (n=163)", "Combo Unique (n=173)")
category   <- factor(master_unique$category, levels = cat_levels)

# ── Colors ────────────────────────────────────────────────────────────────────
# Themed Colormap: Slate Grey (#94a3b8) for Downregulated, White (#ffffff) for Neutral, Purple (#9e68e4) for Upregulated
col_fun <- colorRamp2(c(-3, 0, 3), c("#94a3b8", "#ffffff", "#9e68e4"))

# ── Select 30 Biologically Relevant Genes (Regeneration, Remodeling, Metabolism) ──
genes_to_mark <- c(
  # RARA Unique (n=136): 5 UP, 5 DOWN
  "Elp3", "Pebp1", "Sos1", "Dnm1l", "Mterf3",          # UP
  "Wnt7b", "Lrig1", "Tmem106b", "Cdk6", "Pecr",         # DOWN
  
  # RARG Unique (n=163): 5 UP, 5 DOWN
  "Hras", "Cdc42", "Aars", "Slc25a44", "Adamts18",      # UP
  "Cdc42ep3", "Sema5a", "Dctn6", "Mapk8ip3", "Tmtc4",   # DOWN
  
  # Combo Unique (n=173): 5 UP, 5 DOWN
  "Clu", "Cisd1", "Syngap1", "Pak1ip1", "Apopt1",       # UP
  "Maoa", "Rhobtb3", "Map3k10", "Stambp", "Zbtb26"      # DOWN
)

print("Selected biologically relevant genes to label (10 per unique category):")
print(genes_to_mark)

col_idx <- match(genes_to_mark, colnames(mat_t))

# Create Column Annotation for the bottom of the columns
ha_col <- HeatmapAnnotation(
  link = anno_mark(
    at = col_idx, 
    labels = genes_to_mark, 
    side = "bottom",
    labels_gp = gpar(fontsize = 4.0, fontface = "bold.italic", col = "#1a1a1a"),
    link_height = unit(2.5, "mm"),
    padding = unit(0.4, "mm")
  )
)

# ── Heatmap ───────────────────────────────────────────────────────────────────
ht <- Heatmap(
  mat_t,
  name = "log2(FC)",
  col = col_fun,
  
  # Column split and order
  column_split = category,
  cluster_rows = FALSE,
  row_order = c("RARA", "RARG", "Combo"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  
  # Bottom annotation for column (gene) labels
  bottom_annotation = ha_col,
  
  # Thin outer border only (no lines between individual cells)
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.25),
  rect_gp = gpar(col = NA), # remove cell borders
  
  # Dimensions to fit in 6.8" x 1.25"
  height = unit(11.0, "mm"),   # 3 rows (RARA, RARG, Combo) fit in 11 mm vertical space
  width = unit(13.5, "cm"),   # 472 columns fit in 13.5 cm horizontal space
  
  # Row name styling
  row_names_gp = gpar(fontsize = 5.5, fontface = "bold", col = "#1a1a1a"),
  row_names_side = "left",
  
  # Column titles (categories) styled small
  column_title_gp = gpar(fontsize = 5.0, fontface = "bold", col = "#1a1a1a"),
  column_title_side = "top",
  
  # Column clustering (Show the column hierarchical clustering tree/dendrogram!)
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  column_dend_height = unit(4.0, "mm"),  # Compact dendrogram tree at the top
  column_dend_gp = gpar(lwd = 0.25),
  cluster_column_slices = FALSE,          # Do not cluster the slices
  
  # Legend parameters
  heatmap_legend_param = list(
    title = "log2(FC)",
    title_gp = gpar(fontsize = 5.0, fontface = "bold"),
    labels_gp = gpar(fontsize = 4.5),
    grid_width = unit(2, "mm"),
    legend_height = unit(1.8, "cm")
  )
)

# Save Plot (exactly 6.8" width x 1.25" height, high-DPI publication ready)
out_path <- file.path(OUTPUT_DIR, "Figure_X_I_Unique_horizontal_heatmap.png")
png(out_path, width = 6.8, height = 1.25, units = "in", res = 600)
draw(
  ht, 
  padding = unit(c(1, 2, 5, 2), "mm")
)
dev.off()

print(paste("Horizontal masked heatmap saved successfully to:", out_path))
