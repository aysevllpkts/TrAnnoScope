# Quality check for Illumina short reads. 

import glob

# Define config file:
configfile: "config/config.yaml"

# Create wildcards for Illumina files, SAMPLE refers the names, FRR is the forward and reverse reads.
RNASEQ_INDIR = config["RNASEQ_INDIR"]
SAMPLE, FRR = glob_wildcards(RNASEQ_INDIR + "/{sample}_{frr}.fq.gz")

# OUTPUT project name
#PROJECT = config["PROJECT"]

# Output directory path 
OUTDIR = "output_" + config["PROJECT"]

localrules: all, MultiQC

rule all: 
    "Target files for the workflow."
    input:
        OUTDIR + "/illumina/qc/rawMultiQC"

rule FastQC: 
    "Create quality report for raw Illumina reads."
    input:
        rawdata = RNASEQ_INDIR+"/{sample}_{frr}.fq.gz"
    output:
        OUTDIR + "/illumina/qc/rawFastQC/{sample}_{frr}_fastqc.html",
        OUTDIR + "/illumina/qc/rawFastQC/{sample}_{frr}_fastqc.zip"
    log:
        "logs/illumina_qc/rawFastQC_{sample}_{frr}.log"
    conda:
        "../envs/qc.yaml"
    params:
        outdir = OUTDIR + "/illumina/qc/rawFastQC"
    threads:
        config["fastqc"]["threads"]
    resources:
        mem_mb = config["fastqc"]["memory"]
    benchmark:
        "benchmarks/illumina_qc/rawFastQC_{sample}_{frr}.benchmark.txt"
    priority: 50
    shell:
        """
        echo FastQC is starting for {wildcards.sample}_{wildcards.frr}.fq.gz > {log}
        fastqc -f fastq -o {params.outdir} {input.rawdata} &>> {log}    
        echo DONE... >> {log}
        """
        
rule MultiQC:
    "Generate MultiQC report for raw sample FASTQC reports."
    input: 
        expand(OUTDIR + "/illumina/qc/rawFastQC/{sample}_{frr}_fastqc.zip", sample=SAMPLE, frr=["1","2"])
    output:
        out=directory(OUTDIR + "/illumina/qc/rawMultiQC")
    log:
        "logs/illumina_qc/rawMultiQC.log"
    conda:
        "../envs/qc.yaml"
    threads:
        1
    resources:
        mem_mb=100
    benchmark:
        "benchmarks/illumina_qc/rawMultiQC.benchmark.txt"
    priority: 40
    shell:
        """
        # Create output directory fro MultiQC 
        mkdir -p {output.out} &> {log}

        # Run MultiQC for raw FastQC output
        multiqc -f -o {output.out} {input} -n raw_multiqc_report.html &>> {log} 
        """
