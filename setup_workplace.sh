#!/bin/bash

# Define the root folder
PIPELINE_DIR="short-pipe1"

# Create the folder structure
echo "Creating folder structure for the pipeline..."

mkdir -p $PIPELINE_DIR/{raw_reads,results/{fastqc,trimmed_reads,kraken2,viral_reads,assemblies,centrifuge,deepvirfinder,virsorter2,reports},scripts,config}

# Add placeholders or example files to guide usage
echo "Adding example files..."

# Placeholder for raw reads
touch $PIPELINE_DIR/raw_reads/README.txt
echo "Place raw paired-end FASTQ files here, named as <sample>_R1.fastq.gz and <sample>_R2.fastq.gz." > $PIPELINE_DIR/raw_reads/README.txt

# Placeholder for config.yaml
cat <<EOL > $PIPELINE_DIR/config/config.yaml
samples: ["sample1", "sample2"]  # List your sample names here
threads: 50                     # Default number of threads
kraken2_db: "/path/to/kraken2/db"  # Specify Kraken2 database path
centrifuge_db: "/path/to/centrifuge/db"  # Specify Centrifuge database path
virsorter2_db: "/path/to/virsorter2/db"  # Specify VirSorter2 database path
EOL

# Placeholder for the report generation script
cat <<'EOL' > $PIPELINE_DIR/scripts/generate_report.py
import pandas as pd

# Placeholder script for report generation
# Replace this with your implementation
print("Report generation script goes here.")
EOL

# Instructions for pipeline usage
cat <<EOL > $PIPELINE_DIR/README.txt
Pipeline Directory Structure:
-----------------------------
- raw_reads/          # Place your raw FASTQ files here
- results/            # All pipeline outputs will be stored here
  - fastqc/           # FastQC reports
  - trimmed_reads/    # Trimmed reads from Trim Galore
  - kraken2/          # Kraken2 outputs
  - viral_reads/      # Viral reads filtered by KrakenTools
  - assemblies/       # Contigs assembled by SPAdes
  - centrifuge/       # Centrifuge taxonomical classification results
  - deepvirfinder/    # DeepVirFinder predictions
  - virsorter2/       # VirSorter2 outputs
  - reports/          # Final reports for each sample
- config/             # Configuration files (e.g., config.yaml)
- scripts/            # Custom scripts (e.g., generate_report.py)

Steps to Run the Pipeline:
--------------------------
1. Place raw reads in 'raw_reads/' folder, named as '<sample>_R1.fastq.gz' and '<sample>_R2.fastq.gz'.
2. Update 'config/config.yaml' with your sample names, database paths, and thread settings.
3. Activate the pipeline environment:
   mamba activate viral_pipeline
4. Run the pipeline with:
   snakemake --cores <number_of_cores>

EOL

echo "Pipeline folder structure is ready in '$PIPELINE_DIR'."

