#!/bin/bash

# Check if Docker is installed
if ! command -v docker &> /dev/null
then
    echo "Docker is not installed. Please install Docker before running this script."
    exit 1
fi

# Reference genome of mm10 can be downloaded from "https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads"
# Default parameters
REFERENCE_GENOME="Reference_genome_transgene/Transgene_Reference/"
FASTQ_DIR="Raw_data"
LOCAL_MEM=200
LOCAL_CORES=40
IMAGE_NAME="litd/docker-cellranger"
FORCE_CELLS=5000

# Function to run Cell Ranger
run_cellranger() {
    SAMPLE_ID=$1
    echo "Running Cell Ranger for sample: $SAMPLE_ID"

    docker run --rm -it -v "$PWD":/data -w /data $IMAGE_NAME \
        cellranger count \
        --id="$SAMPLE_ID" \
        --transcriptome="$REFERENCE_GENOME" \
        --fastqs="$FASTQ_DIR/" \
        --sample="$SAMPLE_ID" \
        --create-bam=true \
        --localcores=$LOCAL_CORES \
        --localmem=$LOCAL_MEM \
        --include-introns=true \
        --force-cells=$FORCE_CELLS

    echo "Processing for $SAMPLE_ID completed!"
}

# Check for sample arguments
if [ $# -eq 0 ]; then
    echo "Usage: $0 sample1 [sample2 sample3 ...]"
    exit 1
fi

# Run for each sample provided
for SAMPLE in "$@"
do
    run_cellranger "$SAMPLE"
done

echo "All samples have been processed."
