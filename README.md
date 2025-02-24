# short-pipe-v2
short reads, standard scheme, nt+nr databases (Base Cemtrifuge + Kaiju Indexes)

Pipeline Directory Structure:
-----------------------------
- raw_reads/            
- results/              
  - fastqc/             
  - trimmed_reads/      
  - centrifuge/                    
  - assemblies/                  
  - kaiju/
  - genomad     
  - virsorter2/         
  - reports/            
- config/
  - config.yaml
  - genomad.yaml
  - virsorter2_env.yaml               
- scripts/
  - generate_report2.py              
- tools/
  - VirSorter2/                   
- Snakefile
- PE_Snakefile
- SE_Snakefile

Additional Files:
--------------------------
- setup_environment.sh - Bash script to install the env (including third party tools)

Steps to Run the Pipeline:
--------------------------
1. Place raw reads in 'raw_reads/' folder, named as '<sample>_R1_001.fastq.gz' and '<sample>_R2_001.fastq.gz'.
2. Update 'config/config.yaml' with your sample names and database paths.
3. Activate the pipeline environment:
   `conda activate virpipe_short2`
4. Run the pipeline with:
   `snakemake --cores <number_of_cores> --use-conda`
   - If your mamba version is below 2.* you may use `--conda-frontend mamba` to speed up the subenvs building;
   - You can run single-end library analysis by specifying `--snakefile SE_Snakefile`; default paired-end sequences analysis is copied to PE_Snakefile so you can modify actual Snakefile; 

To be added:
--------------------------
- IMGvr_to_Kaiju.py - python script to transform IMG/VR database as usable inpute for Kaiju custom database
- report_summary.py - python script to make a comprehensive summary for large dataset analysis (one .tsv file summing up {sample}_reports.tsv and MultiQC html)


## DAG
<img src="https://github.com/user-attachments/assets/434f5f39-8d63-4b37-83ba-a635c537d4b7" width="400" height="700" />
