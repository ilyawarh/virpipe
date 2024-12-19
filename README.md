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
- Snakefile

Additional Files:
--------------------------
- setup_environment.sh - Bash script to install the env (including third party tools)
- setup_workplace.sh - Bash script to expand this workplace
- IMGvr_to_centrifuge_DB.py - Python script to transform IMG/VR database files to usable input for "centrifuge-build"  

Steps to Run the Pipeline:
--------------------------
1. Place raw reads in 'raw_reads/' folder, named as '<sample>_R1_001.fastq.gz' and '<sample>_R2_001.fastq.gz'.
2. Update 'config/config.yaml' with your sample names, database paths, and thread settings.
3. Activate the pipeline environment:
   conda activate virpipe_short1
4. Run the pipeline with:
   snakemake --cores <number_of_cores> --use-conda 
