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
        r1_trimmed=f"{RESULTS_DIR}/trimmed_reads/{{sample}}_R1_001_val_1.fq.gz",
        r2_trimmed=f"{RESULTS_DIR}/trimmed_reads/{{sample}}_R2_001_val_2.fq.gz"
    threads: workflow.cores  
    shell:
        "trim_galore --paired --length 21 -j {threads} -o {RESULTS_DIR}/trimmed_reads/ {input.r1} {input.r2}"

rule centrifuge:
    input:
        r1_trimmed=rules.trim_galore.output.r1_trimmed,
	r2_trimmed=rules.trim_galore.output.r2_trimmed
    output:
        reads1=f"{RESULTS_DIR}/centrifuge/{{sample}}/un-conc-mate.1.fastq",
        reads2=f"{RESULTS_DIR}/centrifuge/{{sample}}/un-conc-mate.2.fastq",
#        report=f"{RESULTS_DIR}/centrifuge/{{sample}}/{{sample}}_centrifuge_summary.tsv",
	report=f"{RESULTS_DIR}/centrifuge/{{sample}}/{{sample}}_centrifuge_reads.tsv"
    threads: workflow.cores
    shell:
        "centrifuge -x {config[centrifuge_origin_db]} -1 {input.r1_trimmed} -2 {input.r2_trimmed}  -S {output.report} -q --threads {threads} "
        "--un-conc {RESULTS_DIR}/centrifuge/{wildcards.sample}/ --min-hitlen 30 --host-taxids 9606 --exclude-taxids 10239,2157 && for i in {RESULTS_DIR}/centrifuge/{wildcards.sample}/un-conc-mate*; "
        "do mv $i $i.fastq; done "


rule metaSPAdes:
    input:
        r1_viral=rules.centrifuge.output.reads1,
        r2_viral=rules.centrifuge.output.reads2
    output:
        assembly=f"{RESULTS_DIR}/assemblies/{{sample}}/contigs.fasta"
    threads: workflow.cores
    shell:
        "spades.py --meta -1 {input.r1_viral} -2 {input.r2_viral} -o {RESULTS_DIR}/assemblies/{wildcards.sample}/ "
        "--threads {threads} --only-assembler"

rule kaiju:
    input:
        assembly=rules.metaSPAdes.output.assembly
    output:
        preport=f"{RESULTS_DIR}/kaiju/{{sample}}/{{sample}}_kaiju_report.tsv",
        summary=f"{RESULTS_DIR}/kaiju/{{sample}}/{{sample}}_kaiju_summary.tsv"
    threads: workflow.cores
    shell:
        "kaiju -t {config[kaiju_dmps]}/nodes.dmp -f {config[kaiju_db]} -i {input.assembly} -z {threads} -o {output.preport} "
        "&& kaiju-addTaxonNames -t {config[kaiju_dmps]}/nodes.dmp -n {config[kaiju_dmps]}/names.dmp -i {output.preport} -o {output.summary} "  


rule genomad:
    input:
        assembly=rules.metaSPAdes.output.assembly
    output: 
        scores=f"{RESULTS_DIR}/genomad/{{sample}}/contigs_summary/contigs_virus_summary.tsv"
#    threads: workflow.cores
    conda: 
        "config/genomad.yaml"
    shell:
        "genomad end-to-end --cleanup {input.assembly} {RESULTS_DIR}/genomad/{wildcards.sample} {config[genomad_db]} "

rule virsorter2:
    input:
        assembly=rules.metaSPAdes.output.assembly
    output: 
        file=f"{RESULTS_DIR}/virsorter2/{{sample}}/final-viral-score.tsv"
    threads: workflow.cores
    conda: 
        "config/virsorter2_env.yaml"
    shell: "virsorter config --init-source --db-dir={config[virsorter2_db]} && virsorter run -i {input.assembly} -w {RESULTS_DIR}/virsorter2/{wildcards.sample} --db-dir {config[virsorter2_db]} --min-score 0.5 -j {threads} --include-groups 'RNA,ssDNA,NCLDV' --provirus-off"

rule generate_report:
    input:
        kaiju=f"{RESULTS_DIR}/kaiju/{{sample}}/{{sample}}_kaiju_summary.tsv",
        virsorter2=f"{RESULTS_DIR}/virsorter2/{{sample}}/final-viral-score.tsv",
        genomad=f"{RESULTS_DIR}/genomad/{{sample}}/contigs_summary/contigs_virus_summary.tsv"
    output:
        report=f"{RESULTS_DIR}/reports/{{sample}}_report.tsv"
    script:
        "scripts/generate_report2.py"
