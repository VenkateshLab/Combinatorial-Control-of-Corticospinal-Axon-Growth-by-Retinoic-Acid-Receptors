#!/usr/bin/env Rscript

# Load necessary libraries
suppressMessages(library(ChIPseeker))
suppressMessages(library(GenomicRanges))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
suppressMessages(library(data.table))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript annotate_peaks.R inputfile.xls outputfile.csv")
}

input_file <- args[1]
output_file <- args[2]

# Read the file, manually filtering out lines that start with '#'
lines <- readLines(input_file)
data_lines <- lines[!grepl("^#", lines)]

temp_file <- tempfile()
writeLines(data_lines, temp_file)

peak_data <- fread(temp_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Convert to GenomicRanges object
peaks <- makeGRangesFromDataFrame(df = peak_data, 
                                  keep.extra.columns = TRUE,
                                  seqnames.field = "chr", 
                                  start.field = "start", 
                                  end.field = "end")

# Annotate peaks
peakAnno <- annotatePeak(peaks, 
                         TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                         annoDb = "org.Mm.eg.db")

# Convert annotation to data frame
peakAnno_df <- as.data.frame(peakAnno)

# Save to CSV
write.csv(peakAnno_df, file = output_file, row.names = FALSE)
cat("Done: ", output_file, "\n")
