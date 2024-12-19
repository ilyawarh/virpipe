# short-pipe-v1
short reads, standard scheme, nt database (IMG/VR)

Pipeline Directory Structure:
-----------------------------
- raw_reads/            
- results/              
  - fastqc/             
  - trimmed_reads/      
  - kraken2/            
  - viral_reads/        
  - assemblies/         
  - centrifuge/         
  - deepvirfinder/      
  - virsorter2/         
  - reports/            
- config/
  - config.yaml
  - deepvirfinder_env.yaml
  - virsorter2_env.yaml               
- scripts/
  - generate_report.py              
- tools/
  - KrakenTools/
  - DeepVirFinder/
  - VirSorter2/   	                
- Snakefile

Additional Files:
--------------------------
- setup_environment.sh - Bash script to install and activate the env (including third party tools)
- IMGvr_to_centrifuge_DB.py - Python script to transform IMG/VR database files as usable input for "centrifuge-build"  

Steps to Run the Pipeline:
--------------------------
1. Place raw reads in 'raw_reads/' folder, named as '<sample>_R1_001.fastq.gz' and '<sample>_R2_001.fastq.gz'.
2. Update 'config/config.yaml' with your sample names and database paths.
3. Activate the pipeline environment:
   `conda activate virpipe_short1`
4. Run the pipeline with:
   `snakemake --cores <number_of_cores> --use-conda`
   - you may use `--conda-frontend mamba` if your mamba version is below 2.*

  
<img src="https://github.com/user-attachments/assets/0b349893-86c5-4486-8860-b05c79730d67" width="450" heigt="700" />
