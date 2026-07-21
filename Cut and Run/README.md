# RARA and RARG CUT&RUN Analysis & Plotting Pipeline

This repository contains the complete computational workflow and figure generation code for Retinoic Acid Receptor Alpha (RARA) and Retinoic Acid Receptor Gamma (RARG) CUT&RUN profiling in mouse cell lines. 

The pipeline spans from raw sequencing reads (`.fastq.gz`) quality control, genome alignment, peak calling, and peak annotations, up to reproducing the final publication-ready figures (**Figure X, Panels B to J**).

---

## Directory Structure

```directory
github/Cut and Run/
├── README.md                          # This documentation file
├── cut_and_run_analysis_pipeline.sh   # Bash pipeline script (Raw data -> Annotated Peaks)
├── annotate_peaks.R                  # R script for ChIPseeker peak annotation
├── input_data/                        # Self-contained source datasets
│   ├── Peaks_Annotated_RARA_vs_Input_p1e4.csv
│   ├── Peaks_Annotated_RARG_vs_Input_p1e4.csv
│   ├── RARA_CutandRun_GO.csv
│   ├── RARG_CutandRun_GO.csv
│   ├── Combo_Direct_Activated.csv
│   ├── Combo_Direct_Repressed.csv
│   ├── DEG_Neuronal_RARA_vs_GFP_minpct5_noLogFC.csv
│   ├── DEG_Neuronal_RARG_vs_GFP_minpct5_noLogFC.csv
│   ├── DEG_Neuronal_RARA_RARG_vs_GFP_minpct5_noLogFC.csv
│   ├── overlaps_p_1e-4_genes.txt
│   ├── tss_matrix.gz
│   ├── genes_subset.gtf               # Subsetted mm10 RefGene GTF (7 plotted genes)
│   ├── RARA_Treated_rep1.bw           # BigWig track files
│   ├── RARA_Input_rep1.bw
│   ├── RARG_Treated_rep1.bw
│   ├── RARG_Input_rep1.bw
│   ├── rara_motifs.bed9               # Motif match tracks
│   ├── rarg_motifs.bed9
│   └── logos_png/                     # JASPAR Motif logo PNGs
│       ├── jaspar_MA0159.1.png
│       ├── jaspar_MA0512.2.png
│       ├── jaspar_MA0730.1.png
│       └── jaspar_MA0860.1.png
├── output_figures/                    # Plotted high-resolution images
│   ├── FigureX_B.png
│   ├── FigureX_C.png
│   ├── FigureX_D.png
│   ├── FigureX_E.png
│   ├── FigureX_F.png
│   ├── FigureX_G.png
│   ├── FigureX_H.png
│   ├── FigureX_I_*.png                # 7 genomic screenshots + 1 legend
│   ├── FigureX_J.png                  # GO-grouped horizontal heatmap
│   └── figure_legend.md               # Detailed Figure X panel legend
└── plotting_scripts/                  # Python figure-rendering scripts
    ├── generate_figure_B.py
    ├── generate_figure_C.py
    ├── generate_figure_D.py
    ├── generate_figure_E.py
    ├── generate_figure_F.py
    ├── generate_figure_G.py
    ├── generate_figure_H.py
    ├── generate_figure_I.py
    └── generate_figure_J.py
```

---

## Part 1: Primary Data Processing Pipeline

All primary data processing steps are documented and compiled in the executable shell script `cut_and_run_analysis_pipeline.sh`. The workflow consists of the following steps:

1. **Quality Control on Raw Reads (`FastQC v0.11.9`):** Assesses sequencing quality of raw paired-end fastq files.
2. **Adapter Trimming and Read Filtering (`fastp v0.23.2`):** Trims sequencing adapters, filters out reads with low-quality bases (Phred < 20), and filters out short fragments (< 25bp).
3. **Genome Alignment (`Bowtie2 v2.4.4`):** Align trimmed paired-end reads to the reference mouse genome (`mm10`) index. Optimized parameters for CUT&RUN:
   ```bash
   bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700
   ```
4. **BAM Processing (`Samtools v1.13`):** Converts SAM alignments to BAM format, sorts alignments by genomic coordinates, and indexes the final sorted BAM files.
5. **Peak Calling (`MACS2 v2.2.7.1`):** Identifies RARA and RARG enriched peaks relative to the shared Input control sample in paired-end mode with a p-value threshold of $p < 10^{-4}$:
   ```bash
   macs2 callpeak -t treatment.bam -c input.bam -f BAMPE -g mm -p 1e-4 -n RARA_vs_Input
   ```
6. **Peak Genomic Feature Annotation (`ChIPseeker v1.28.3`):** Annotates called peaks with nearest genomic features (Promoter, Intron, Exon, UTR, Distal Intergenic) relative to the mm10 refGene transcript model database. Executed using:
   ```bash
   Rscript annotate_peaks.R peaks.xls annotated_peaks.csv
   ```

---

## Part 2: Reproducing the Figures (Panels B to J)

### Python Dependencies

To run the plotting scripts and regenerate the final figures, ensure the following Python packages are installed:

```bash
pip install pandas numpy matplotlib scipy pyBigWig pyGenomeTracks
```

*Note: `pyGenomeTracks` is required for Panel I track snapshots.*

### How to Run

To regenerate a specific panel, run the corresponding python script inside the `plotting_scripts/` directory. All scripts load data relatively from `../input_data/` and output figures directly into `../output_figures/`.

For example:
```bash
cd plotting_scripts
python3 generate_figure_B.py
python3 generate_figure_J.py
```

### Figure-to-Code Mapping:

| Panel | Figure Filename | Description | Code File |
| :--- | :--- | :--- | :--- |
| **B** | `FigureX_B.png` | Peak Genomic Distribution Bar Plot | `generate_figure_B.py` |
| **C** | `FigureX_C.png` | TSS Profile Metagene Line Plot | `generate_figure_C.py` |
| **D** | `FigureX_D.png` | Box + Violin Peak Fold Enrichment | `generate_figure_D.py` |
| **E** | `FigureX_E.png` | Co-bound Peak Fold Enrichment | `generate_figure_E.py` |
| **F** | `FigureX_F.png` | Peak Enrichment & snRNA Expression | `generate_figure_F.py` |
| **G** | `FigureX_G.png` | direct binding percentage in DEGs | `generate_figure_G.py` |
| **H** | `FigureX_H.png` | JASPAR Motif Occupancy Table | `generate_figure_H.py` |
| **I** | `FigureX_I_*.png` | Zoomed locus snapshots & Legend | `generate_figure_I.py` |
| **J** | `FigureX_J.png` | GO-grouped Peak Fold Enrichment Heatmap | `generate_figure_J.py` |

---

## Figure Legend

Detailed descriptions of each panel can be found inside the text catalog at **`github/output_figures/figure_legend.md`**.
