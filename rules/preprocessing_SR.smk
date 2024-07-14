# Filtering and trimming of Illumina short reads. 

from snakemake.utils import min_version
# Set minimum Snakemake version
min_version('6.0')

import glob

# Define config file:
configfile: "config/config.yaml"

# Create wildcards for Illumina files, SAMPLE refers the names, FRR is the forward and reverse reads.
RNASEQ_INDIR = config["RNASEQ_INDIR"]
SAMPLE, FRR = glob_wildcards(RNASEQ_INDIR + "/{sample}_{frr}.fq.gz")

# OUTDIR path
OUTDIR = "output_" + config["PROJECT"]

localrules: all, trimmedMultiQC

rule all: 
    "Target files for the workflow."
    input:
        OUTDIR + "/illumina/trimmedMultiQC"


rule fastqscreen:
    "Run FastQScreen on each raw Illumina short reads."             # fastqscreen processes PE reads separately  
    input:                                                          # because of it, it can be different number of reads each PE
        read=RNASEQ_INDIR+"/{sample}_{frr}.fq.gz"
    output:
        outfile=OUTDIR + "/illumina/filtered_fastqscreen/{sample}_{frr}.tagged_filter.fq.gz"
    conda:
        "../envs/qc.yaml"
    params:
        aligner="bowtie2",
        config_file=config["fastqscreen"]["conf"],
        outdir=OUTDIR + "/illumina/filtered_fastqscreen",
        uncompressed = "{sample}_{frr}.tagged_filter.fastq"
    log:
        "logs/illumina/fastqscreen_{sample}_{frr}.log"
    threads:
        config["fastqscreen"]["threads"]
    resources:
        mem_mb = config["fastqscreen"]["memory"]
    benchmark:
        "benchmarks/illumina_qc/fastqscreen_{sample}_{frr}.benchmark.txt"
    shell:
        """
        # Create output directory for FastQ Screen 
        mkdir -p {params.outdir} &> {log}

        # Compressed fq.gz - fastqscreen does not work with compressed files
        gunzip -c {input.read} > {params.outdir}/{wildcards.sample}_{wildcards.frr}.fq
        # Run FastQ Screen for raw reads
        fastq_screen {params.outdir}/{wildcards.sample}_{wildcards.frr}.fq --threads {threads} --conf {params.config_file} --aligner {params.aligner} --nohits --outdir {params.outdir} &>> {log}
        # Compressed file again
        gzip -c {params.outdir}/{params.uncompressed} > {output.outfile}
        # Remove temporary uncompressed files 
        rm {params.outdir}/{wildcards.sample}_{wildcards.frr}.fq {params.outdir}/{params.uncompressed}
        """


if config["END"] == "PE":
    rule fastp_PE:                                                         # In order to prevent ERROR for different number of reads in PE
        input:                                                          # version of the fastp is IMPORTANT. I am using 0.23.1 (0.23.4 do not run diffferent read number files!)    
            F=OUTDIR + "/illumina/filtered_fastqscreen/{sample}_1.tagged_filter.fq.gz",
            R=OUTDIR + "/illumina/filtered_fastqscreen/{sample}_2.tagged_filter.fq.gz"
        output:
            #directory("preprocessing/illumina/trimmed_fastp"),
            trimmed_F=OUTDIR + "/illumina/trimmed_fastp/{sample}_1.trimmed.fq.gz",
            trimmed_R=OUTDIR + "/illumina/trimmed_fastp/{sample}_2.trimmed.fq.gz",
            json=OUTDIR + "/illumina/trimmed_fastp/{sample}.json",
            html=OUTDIR + "/illumina/trimmed_fastp/{sample}.html"
        conda:
            "../envs/qc.yaml"
        params:
            params=config["fastp"]["parameters"],
            adapter=config["fastp"]["adapter"],
            outdir = OUTDIR + "/illumina/trimmed_fastp"
        log:
            "logs/illumina/fastp_PE_{sample}.log" 
        threads:
            config["fastp"]["threads"]
        resources:
            mem_mb = config["fastp"]["memory"]
        benchmark:
            "benchmarks/illumina/fastp_PE_{sample}.benchmark.txt"
        shell:
            """
            # Create the output directory for fastp if it doesn't exist
            mkdir -p {params.outdir} &> {log}

            # Run fastp for filtered reads
            fastp -i {input.F} -I {input.R} -o {output.trimmed_F} -O {output.trimmed_R} --thread {threads} --adapter_fasta {params.adapter} --json {output.json} --html {output.html} {params.params} 2> {log} 
            """
elif config["END"] == "SE":
    rule fastp_SE:                                                         # In order to prevent ERROR for different number of reads in PE
        input:                                                          # version of the fastp is IMPORTANT. I am using 0.23.1 (0.23.4 do not run diffferent read number files!)    
            F=OUTDIR + "/illumina/filtered_fastqscreen/{sample}_1.tagged_filter.fastq.gz",
        output:
            #directory("preprocessing/illumina/trimmed_fastp"),
            trimmed_F=OUTDIR + "/illumina/trimmed_fastp/{sample}_1.trimmed.fq.gz",
            json=OUTDIR + "/illumina/trimmed_fastp/{sample}.json",
            html=OUTDIR + "/illumina/trimmed_fastp/{sample}.html"
        conda:
            "../envs/qc.yaml"
        params:
            params = config["fastp"]["parameters"],
            adapter = config["fastp"]["adapter"],
            outdir = OUTDIR + "/illumina/trimmed_fastp"
        log:
            "logs/illumina/fastp_SE_{sample}.fastp.log" 
        threads:
            config["fastp"]["threads"]
        resources:
            mem_mb = config["fastp"]["memory"]
        benchmark:
            "benchmarks/illumina/fastp_SE_{sample}.benchmark.txt"
        shell:
            """
            # Create the output directory for fastp if it doesn't exist
            mkdir -p {params.outdir} &> {log}

            # Run fastp for filtered reads
            fastp -i {input.F} -o {output.trimmed_F} --thread {threads} --adapter_fasta {params.adapter} --json {output.json} --html {output.html} {params.params} &>> {log} 
            """
else:
    print("Choose PE or SE in config file for END variable..")

rule trimmedFastQC:
    input:
        trimmed=OUTDIR + "/illumina/trimmed_fastp/{sample}_{frr}.trimmed.fq.gz"
    output:
        OUTDIR + "/illumina/trimmedFastQC/{sample}_{frr}.trimmed_fastqc.html",
        OUTDIR + "/illumina/trimmedFastQC/{sample}_{frr}.trimmed_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    params:
        outdir=OUTDIR + "/illumina/trimmedFastQC"
    threads:
        config["fastqc"]["threads"]
    resources:
        mem_mb = config["fastqc"]["memory"]
    log:
        "logs/illumina_qc/trimmedFastQC_{sample}_{frr}.log"
    shell:
        """
        # Create directory for FastQC output
        mkdir -p {params.outdir} &> {log}

        # Run FastQC for trimmed reads
        fastqc -f fastq -o {params.outdir} {input.trimmed} -t {threads} %>> {log}
        """

rule trimmedMultiQC:
    input:
        expand(OUTDIR + "/illumina/trimmedFastQC/{sample}_{frr}.trimmed_fastqc.zip", sample=SAMPLE, frr=FRR) 
    output:
        out=directory(OUTDIR + "/illumina/trimmedMultiQC")
    conda:
        "../envs/qc.yaml"
    threads:
        1
    resources:  
        mem_mb = 100
    log:
        "logs/illumina_qc/trimmed_MultiQC.log"
    shell:
        """
        # Create a directory for MultiQC output
        mkdir -p {output.out} &> {log}

        # Run MultiQC in trimmed FastQC outputs
        multiqc -o {output.out} {input} -n trimmed_multiqc_report.html &>> {log}
        """


