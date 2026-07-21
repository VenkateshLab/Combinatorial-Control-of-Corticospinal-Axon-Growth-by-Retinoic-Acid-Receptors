**Combinatorial Control of Corticospinal Axon Growth by Retinoic Acid Receptors**  
This repository contains computational workflows supporting the study **“Combinatorial Control of Corticospinal Axon Growth by Retinoic Acid Receptors.”** The analyses examine how retinoic acid receptor alpha ( **RARA**), retinoic acid receptor gamma ( **RARG**), and their combined activity influence corticospinal neuron gene regulation and axon-growth-associated programs.  
Two complementary sequencing approaches are included:  
- **CUT&RUN** profiling to identify genomic regions bound by RARA and RARG.  
- **Single-nucleus RNA sequencing (snRNA-seq)** to characterize transcriptional responses to RARA, RARG, and combined RARA/RARG expression relative to GFP controls.  
Together, these datasets connect receptor occupancy with changes in gene expression to identify candidate direct targets and biological pathways involved in corticospinal axon growth.  
**Repository structure**  
.  
 ├── README.md  
 ├── Methods_Draft.md  
 ├── Cut and Run/  
 │   ├── README.md  
 │   ├── cut_and_run_analysis_pipeline.sh  
 │   └── annotate_peaks.R  
 └── snRNA/  
     ├── README.md  
     └── CellRanger.sh  
   
**Analysis workflows**  
**CUT&RUN**  
The CUT&RUN workflow processes paired-end sequencing data for RARA, RARG, and input-control samples. It includes:  
1. Raw-read quality control with FastQC.  
2. Adapter trimming and read filtering with fastp.  
3. Alignment to the mouse mm10 genome with Bowtie2.  
4. BAM sorting and indexing with SAMtools.  
5. Peak calling against input controls with MACS2.  
6. Peak annotation with ChIPseeker in R.  
Detailed instructions and parameters are provided in the [CUT&RUN README.](Cut%20and%20Run/README.md "Cut%20and%20Run/README.md")  
**Single-nucleus RNA sequencing**  
The snRNA-seq workflow uses a custom mouse reference containing AAV-GFP, AAV-RARA, and AAV-RARG transgene sequences. The repository includes a Docker-based Cell Ranger script for generating gene–barcode count matrices from FASTQ files. Downstream analyses compare RARA, RARG, and combined RARA/RARG conditions with GFP controls and include differential-expression, gene-overlap, heatmap, and Gene Ontology analyses.  
Detailed analysis and figure descriptions are provided in the [snRNA-seq README.](snRNA/README.md "snRNA/README.md")  
**Getting started**  
Clone the repository and enter its directory:  
git clone <https://github.com/VenkateshLab/Combinatorial-Control-of-Corticospinal-Axon-Growth-by-Retinoic-Acid-Receptors>  
 cd <rCombinatorial-Control-of-Corticospinal-Axon-Growth-by-Retinoic-Acid-Receptors>  
   
Before running either workflow, update the input-data and reference-genome paths in the relevant shell script. The scripts assume paired-end FASTQ input and access to the indicated command-line tools.  
Run the CUT&RUN workflow from its directory:  
cd "Cut and Run"  
 bash cut_and_run_analysis_pipeline.sh  
   
Run Cell Ranger for one or more snRNA-seq samples from the snRNA directory:  
cd snRNA  
 bash CellRanger.sh <sample_id> [additional_sample_ids]  
   
The Cell Ranger script currently requests 40 CPU cores and 200 GB of memory. Adjust LOCAL_CORES, LOCAL_MEM, FORCE_CELLS, the Docker image, and input/reference paths to match the computing environment and experiment.  
**Main software requirements**  
- FastQC  
- fastp  
- Bowtie2  
- SAMtools  
- MACS2  
- R with Bioconductor packages ChIPseeker, TxDb.Mmusculus.UCSC.mm10.knownGene, and org.Mm.eg.db  
- Docker and 10x Genomics Cell Ranger  
- Python and R packages required by the downstream plotting scripts described in the module READMEs  
Software versions and analysis parameters should be kept consistent with those documented in the individual workflows.  
**Data and reproducibility**  
Large raw sequencing files and reference-genome resources are not stored in this repository. Place data in the locations expected by each workflow, or revise the directory variables in the scripts. The analysis scripts are intended to document the computational workflow and may require adaptation to local file names, sample identifiers, and computing resources.  
The current manuscript-oriented methods summary is available in [Methods_Draft.md. Items in square brackets are placeholders that should be completed before publication.](Methods_Draft.md "Methods_Draft.md")  
**Citation**  
If you use this code or associated data, please cite the accompanying publication. Full citation details will be added when available.  
**Contact**  
For questions about the analysis or data, please contact the repository authors through the corresponding GitHub repository.  
