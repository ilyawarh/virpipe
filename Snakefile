configfile: "config/config.yaml"

SAMPLES = config["samples"]
RESULTS_DIR = "results"
#THREADS = config["threads"]


rule all:
	input:
	  expand(f"{RESULTS_DIR}/reports/{{sample}}_report.tsv", sample=SAMPLES)
#	  expand(f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_R1_fastqc.html", sample=SAMPLES)


rule fastqc_before:
    input:
        r1="raw_reads/{sample}_R1_001.fastq.gz",
        r2="raw_reads/{sample}_R2_001.fastq.gz"
    output:
        html=[f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_R1_fastqc.html", f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_R2_fastqc.html"],
        zip=[f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_R1_fastqc.zip", f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_R2_fastqc.zip"]
    shell:
        "fastqc {input.r1} {input.r2} -o {RESULTS_DIR}/fastqc/ -t 8"

rule trim_galore:
    input:
        r1="raw_reads/{sample}_R1_001.fastq.gz",
        r2="raw_reads/{sample}_R2_001.fastq.gz"
    output:
        r1_trimmed=f"{RESULTS_DIR}/trimmed_reads/{{sample}}_R1_val_1.fq.gz",
        r2_trimmed=f"{RESULTS_DIR}/trimmed_reads/{{sample}}_R2_val_2.fq.gz"
    shell:
        "trim_galore --paired --length 21 -j 16 -o {RESULTS_DIR}/trimmed_reads/ {input.r1} {input.r2}"

rule kraken2:
    input:
        r1=rules.trim_galore.output.r1_trimmed,
        r2=rules.trim_galore.output.r2_trimmed
    output:
        report=f"{RESULTS_DIR}/kraken2/{{sample}}/{{sample}}_kraken2.report",
        classified=f"{RESULTS_DIR}/kraken2/{{sample}}/{{sample}}_classified.fastq"
    threads: workflow.cores
    shell:
        "kraken2 --paired --db {config[kraken2_db]} --threads {threads} --confidence 0.3 "
        "--report {output.report} --output {output.classified} {input.r1} {input.r2}"

rule extract_viral_reads:
    input:
        report=rules.kraken2.output.report,
        classified=rules.kraken2.output.classified
    output:
        r1_viral=f"{RESULTS_DIR}/viral_reads/{{sample}}_R1_viral.fastq",
        r2_viral=f"{RESULTS_DIR}/viral_reads/{{sample}}_R2_viral.fastq"
    shell:
        "tools/KrakenTools/extract_kraken_reads.py -k {input.classified} -1 {RESULTS_DIR}/trimmed_reads/{wildcards.sample}_R1_val_1.fq.gz -2 {RESULTS_DIR}/trimmed_reads/{wildcards.sample}_R2_val_2.fq.gz -r {input.report} --include-children "
        "--exclude --taxid 2 2759 -o {output.r1_viral} -o2 {output.r2_viral} --fastq-output"

rule metaSPAdes:
    input:
        r1_viral=rules.extract_viral_reads.output.r1_viral,
        r2_viral=rules.extract_viral_reads.output.r2_viral
    output:
        assembly=f"{RESULTS_DIR}/assemblies/{{sample}}/contigs.fasta"
    threads: workflow.cores
    shell:
        "spades.py --meta -1 {input.r1_viral} -2 {input.r2_viral} -o {RESULTS_DIR}/assemblies/{wildcards.sample}/ "
        "--threads {threads} --only-assembler"

rule centrifuge:
    input:
        assembly=rules.metaSPAdes.output.assembly
    output:
        report=f"{RESULTS_DIR}/centrifuge/{{sample}}/{{sample}}_centrifuge_summary.tsv",
	reads=f"{RESULTS_DIR}/centrifuge/{{sample}}/{{sample}}_centrifuge_reads.tsv"
    threads: workflow.cores
    shell:
        "centrifuge -x {config[centrifuge_db]} -U {input.assembly} --report-file {output.report}  -S {output.reads} --threads {threads} -f"

rule deepvirfinder:
    input:
        assembly=rules.metaSPAdes.output.assembly
    output: scores=f"{RESULTS_DIR}/deepvirfinder/{{sample}}/contigs.fasta_gt1bp_dvfpred.txt"
    threads: workflow.cores
    conda: "config/deepvirfinder_env.yaml"
    shell:
        "python tools/DeepVirFinder/dvf.py -i {input.assembly} -o {RESULTS_DIR}/deepvirfinder/{wildcards.sample} -c {threads}"

rule virsorter2:
    input:
        assembly=rules.metaSPAdes.output.assembly
    output: file=f"{RESULTS_DIR}/virsorter2/{{sample}}/final-viral-score.tsv"
    threads: workflow.cores
    conda: "config/virsorter2_env.yaml"
    shell: "virsorter config --init-source --db-dir={config[virsorter2_db]} && virsorter run -i {input.assembly} -w {RESULTS_DIR}/virsorter2/{wildcards.sample} --db-dir {config[virsorter2_db]} --min-score 0.5 -j {threads} --include-groups RNA"

rule generate_report:
    input:
        centrifuge=f"{RESULTS_DIR}/centrifuge/{{sample}}/{{sample}}_centrifuge_reads.tsv",
        deepvirfinder=f"{RESULTS_DIR}/deepvirfinder/{{sample}}/contigs.fasta_gt1bp_dvfpred.txt",
        virsorter2=f"{RESULTS_DIR}/virsorter2/{{sample}}/final-viral-score.tsv"
    output:
        report=f"{RESULTS_DIR}/reports/{{sample}}_report.tsv"
    script:
        "scripts/generate_report.py"
