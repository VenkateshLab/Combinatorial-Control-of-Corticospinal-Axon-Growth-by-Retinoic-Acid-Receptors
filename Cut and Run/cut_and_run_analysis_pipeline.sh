#!/usr/bin/env bash

# ==============================================================================
# CUT&RUN Bioinformatics Analysis Pipeline Script
# ==============================================================================
# This script outlines the complete computational workflow used for processing 
# RARA and RARG CUT&RUN data, including quality control, adapter trimming, 
# alignment to the reference genome, peak calling, and peak genomic annotation.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status

# Define directory variables
RAW_DATA_DIR="raw_data"
QC_DIR="qc_reports"
TRIM_DIR="trimmed_data"
ALIGN_DIR="alignments"
PEAK_DIR="peaks_output"

# Reference Genome Index path
BOWTIE2_INDEX="reference/mm10" # Path to Bowtie2 index for mouse (mm10)

echo "Starting CUT&RUN Analysis Pipeline..."

# ------------------------------------------------------------------------------
# STEP 1: Quality Control on Raw Fastq Files
# ------------------------------------------------------------------------------
echo "Step 1: Running FastQC on raw reads..."
mkdir -p ${QC_DIR}
fastqc -o ${QC_DIR} ${RAW_DATA_DIR}/*.fastq.gz

# ------------------------------------------------------------------------------
# STEP 2: Adapter Trimming and Read Filtering (fastp)
# ------------------------------------------------------------------------------
echo "Step 2: Running fastp for adapter trimming and low-quality filtering..."
mkdir -p ${TRIM_DIR}

# Example processing RARA treatment sample
fastp \
    -i ${RAW_DATA_DIR}/RARA_treated_R1.fastq.gz \
    -I ${RAW_DATA_DIR}/RARA_treated_R2.fastq.gz \
    -o ${TRIM_DIR}/RARA_treated_trimmed_R1.fastq.gz \
    -O ${TRIM_DIR}/RARA_treated_trimmed_R2.fastq.gz \
    -h ${QC_DIR}/RARA_treated_fastp.html \
    -j ${QC_DIR}/RARA_treated_fastp.json \
    --detect_adapter_for_pe \
    --qualified_quality_filter --qualified_phred 20 --length_required 25

# Example processing RARG treatment sample
fastp \
    -i ${RAW_DATA_DIR}/RARG_treated_R1.fastq.gz \
    -I ${RAW_DATA_DIR}/RARG_treated_R2.fastq.gz \
    -o ${TRIM_DIR}/RARG_treated_trimmed_R1.fastq.gz \
    -O ${TRIM_DIR}/RARG_treated_trimmed_R2.fastq.gz \
    -h ${QC_DIR}/RARG_treated_fastp.html \
    -j ${QC_DIR}/RARG_treated_fastp.json \
    --detect_adapter_for_pe \
    --qualified_quality_filter --qualified_phred 20 --length_required 25

# Example processing Input control sample
fastp \
    -i ${RAW_DATA_DIR}/Input_R1.fastq.gz \
    -I ${RAW_DATA_DIR}/Input_R2.fastq.gz \
    -o ${TRIM_DIR}/Input_trimmed_R1.fastq.gz \
    -O ${TRIM_DIR}/Input_trimmed_R2.fastq.gz \
    -h ${QC_DIR}/Input_fastp.html \
    -j ${QC_DIR}/Input_fastp.json \
    --detect_adapter_for_pe \
    --qualified_quality_filter --qualified_phred 20 --length_required 25

# ------------------------------------------------------------------------------
# STEP 3: Post-Trimming Quality Control
# ------------------------------------------------------------------------------
echo "Step 3: Running post-trimming FastQC..."
fastqc -o ${QC_DIR} ${TRIM_DIR}/*.fastq.gz

# ------------------------------------------------------------------------------
# STEP 4: Alignment to Reference Genome (Bowtie2)
# ------------------------------------------------------------------------------
# Bowtie2 settings optimized for CUT&RUN:
# --local --very-sensitive-local: allows soft clipping of adapters
# --no-mixed --no-discordant: forces paired-end alignment of fragments
# -I 10 -X 700: sets fragment size ranges (10bp to 700bp)
# ------------------------------------------------------------------------------
echo "Step 4: Aligning trimmed reads using Bowtie2 and converting to sorted BAM..."
mkdir -p ${ALIGN_DIR}

for sample in RARA_treated RARG_treated Input; do
    echo "Aligning sample: ${sample}..."
    bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 \
        -I 10 -X 700 -x ${BOWTIE2_INDEX} \
        -1 ${TRIM_DIR}/${sample}_trimmed_R1.fastq.gz \
        -2 ${TRIM_DIR}/${sample}_trimmed_R2.fastq.gz \
        -S ${ALIGN_DIR}/${sample}.sam
        
    # Convert to BAM, sort, and index
    samtools view -bS ${ALIGN_DIR}/${sample}.sam | samtools sort -o ${ALIGN_DIR}/${sample}_sorted.bam
    samtools index ${ALIGN_DIR}/${sample}_sorted.bam
    
    # Remove large intermediate SAM file
    rm ${ALIGN_DIR}/${sample}.sam
done

# ------------------------------------------------------------------------------
# STEP 5: Peak Calling (MACS2)
# ------------------------------------------------------------------------------
# Calls peaks for RARA and RARG against Input background control.
# Uses -f BAMPE flag for paired-end mode (which reads fragment sizes directly).
# Set p-value threshold to 1e-4 (-p 1e-4) to call enriched regions.
# ------------------------------------------------------------------------------
echo "Step 5: Calling peaks against Input using MACS2..."
mkdir -p ${PEAK_DIR}

# Call RARA peaks
macs2 callpeak \
    -t ${ALIGN_DIR}/RARA_treated_sorted.bam \
    -c ${ALIGN_DIR}/Input_sorted.bam \
    -f BAMPE -g mm -p 1e-4 \
    -n RARA_vs_Input --outdir ${PEAK_DIR}

# Call RARG peaks
macs2 callpeak \
    -t ${ALIGN_DIR}/RARG_treated_sorted.bam \
    -c ${ALIGN_DIR}/Input_sorted.bam \
    -f BAMPE -g mm -p 1e-4 \
    -n RARG_vs_Input --outdir ${PEAK_DIR}

# ------------------------------------------------------------------------------
# STEP 6: Genomic Feature Annotation (ChIPseeker)
# ------------------------------------------------------------------------------
# Annotate the peaks with promoters, introns, exons, distal intergenic regions 
# using ChIPseeker in R.
# ------------------------------------------------------------------------------
echo "Step 6: Annotating peaks using ChIPseeker in R..."
Rscript annotate_peaks.R ${PEAK_DIR}/RARA_vs_Input_peaks.xls ${PEAK_DIR}/Peaks_Annotated_RARA_vs_Input_p1e4.csv
Rscript annotate_peaks.R ${PEAK_DIR}/RARG_vs_Input_peaks.xls ${PEAK_DIR}/Peaks_Annotated_RARG_vs_Input_p1e4.csv

echo "CUT&RUN Analysis Pipeline completed successfully!"
