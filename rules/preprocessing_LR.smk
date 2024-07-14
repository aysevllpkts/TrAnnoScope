# Aysevil Pektas
# 12-03-2024
# isoseq3 preprocessing snakefile

import glob
import os 

######## NOTES #########
# Add --latency-wait 
# If you change something in environment.yaml, it will rerun all the pre-obtained files.
# core-dump can be due to file naming problems. 
# for lima order is input primer OUTDIR in the command
# do not use same variable name for input and OUTDIR. (for example, fl cannot be used as a variable in OUTDIR in one rule and input in another)
########################

######## RULES #########
# subreads to ccs
# ccs to fl
# fl to flnc
# merge flnc bam files
# flnc to cluster
########################

# Define config file:
configfile: "config/config.yaml"

## Define wildcards:
ISOSEQ_DATA = config["ISOSEQ_INDIR"]
IDS,=glob_wildcards(ISOSEQ_DATA+"/{id}.subreads.bam")

# Obtain unique sample names:
SAMPLES, MOVIE = glob_wildcards(ISOSEQ_DATA+"/{sample}_{movie}.subreads.bam") # movie is not necessary.
SAMPLES = set(SAMPLES)

# OUTDIR directory for this SMK.
OUTDIR = "output_" + config["PROJECT"] + "/pacbio"

# Function to add suffix for primers that used in lima (FL)
"""
def add_suffix(primer_file):
    fasta=open(primer_file)
    head=[]
    for line in fasta:
        if line[0] == '>':
            line=line.strip()
            head.append(line[1:len(line)])
    return "--".join(head)

"""

# For CCS2FL rule, we need suffix for the file
def add_suffix(primer_file):
    with open(primer_file, 'r') as fasta:
        head = [line.strip()[1:] for line in fasta if line.startswith('>')]
    return "--".join(head)

primer_suffix = add_suffix(config["primers"])


rule all:   
    input: 
        expand(OUTDIR + "/ccs/{id}.ccs.bam", id=IDS),
        expand(OUTDIR + "/fl/{id}.fl.{primer_suffix}.bam", id=IDS, primer_suffix=primer_suffix),
        expand(OUTDIR + "/flnc/{id}.flnc.bam", id=IDS),
        expand(OUTDIR + "/flnc/{sample}.flnc.fofn", sample=SAMPLES),
        expand(OUTDIR + "/clustered/{sample}.clustered.bam", sample=SAMPLES)

# Rules:
rule subreads2ccs:
    input:
        subreads = ISOSEQ_DATA +"/{id}.subreads.bam" 
    output:
        ccs = OUTDIR + "/ccs/{id}.ccs.bam",
        report = OUTDIR + "/ccs/{id}.ccs_report.txt"
    log: 
        "logs/pacbio/subreads2ccs.{id}.log"
    params:
        min_rq=config["subreads2ccs"]["min_rq"],
        min_passes=config["subreads2ccs"]["min_passes"]
    threads:
        config["subreads2ccs"]["threads"]
    conda:
        "../envs/pacbio.yaml"
    resources: 
        mem_mb = config["subreads2ccs"]["memory"]
    shell:
        """
        ccs --min-passes {params.min_passes} --min-rq {params.min_rq} --report-file {output.report} --log-file {log} -j {threads} {input} {output.ccs}
        """

rule ccs2fl:
    input:
        ccs = OUTDIR + "/ccs/{id}.ccs.bam",
        primers = config["primers"]
    output:
        fl = OUTDIR + "/fl/{id}.fl.{primer_suffix}.bam"
    log:
        "logs/pacbio/ccs2fl.{id}.{primer_suffix}.log"
    conda:
        "../envs/pacbio.yaml"
    params:
        outfile = OUTDIR + "/fl/{id}.fl.bam"
    resources: 
        mem_mb = config["ccs2fl"]["memory"]
    threads:
        config["ccs2fl"]["threads"]
    shell:
        """
        lima --isoseq --log-file {log} {input.ccs} {input.primers} {params.outfile} --num-threads {threads}
        """

rule fl2flnc:
    input:
        fl = expand(OUTDIR + "/fl/{id}.fl.{primer_suffix}.bam", id = IDS, primer_suffix = primer_suffix)
    output:
        flnc = OUTDIR + "/flnc/{id}.flnc.bam"
    log:
        "logs/pacbio/fl2flnc_isoseq3_refine.{id}.log"
    conda:
        "../envs/pacbio.yaml"
    params: 
        indir = OUTDIR + "/fl",
        primers = config["primers"] 
    threads:
        config["fl2flnc"]["threads"]
    resources: 
        mem_mb = config["fl2flnc"]["memory"]
    shell:
        """
        isoseq3 refine --require-polya --log-file {log} {params.indir}/{wildcards.id}.fl.{primer_suffix}.bam {params.primers} {output.flnc} -j {threads}
        """

rule create_dataset:
    input:
        file = expand(OUTDIR + "/flnc/{id}.flnc.bam", id=IDS),
    output:
        dataset = OUTDIR + "/flnc/{sample}.flnc.fofn"
    log:
        "logs/pacbio/create_dataset.{sample}.log"
    params:
        pattern = "{sample}_*.flnc.bam",
        outfile = "{sample}.flnc.fofn",
        dir = OUTDIR + "/flnc"
    threads:
        1
    resources:
        mem_mb = 100
    shell:
        """
        cd {params.dir} 
        ls {params.pattern} > {params.outfile} 
        cd ../../../..
        """
# Here you go
if config["flnc2cluster"]["cluster_type"] == "cluster2":
    rule flnc2cluster2:
        input:
            dataset = OUTDIR + "/flnc/{sample}.flnc.fofn"
        output:
            clustered = OUTDIR + "/clustered/{sample}.clustered.bam",
            fasta = OUTDIR + "/clustered/{sample}.clustered.fasta"
        log:
            "logs/pacbio/cluster2.{sample}.log"
        conda:
            "../envs/pacbio.yaml"
        threads:
            config["flnc2cluster"]["threads"]
        resources: 
            mem_mb = config["flnc2cluster"]["memory"]
        shell:  
            """ 
            echo "isoseq3 cluster2 was used for clustering." > {log} 
            isoseq3 cluster2 --log-file {log} {input.dataset} {output.clustered} -j {threads} >> {log}
            bamtools convert -format fasta -in {output.clustered} -out {output.fasta} >> {log}
            """

elif config["flnc2cluster"]["cluster_type"] == "cluster":
    rule flnc2cluster:
        input:
            dataset = OUTDIR + "/flnc/{sample}.flnc.fofn"
        output:
            clustered = OUTDIR + "/clustered/{sample}.clustered.bam"
        params:
            fasta = OUTDIR + "/clustered/{sample}.clustered.hq.fasta.gz"
        threads:
            config["flnc2cluster"]["threads"]
        resources: 
            mem_mb = config["flnc2cluster"]["memory"]    
        shell:
            """
            echo "isoseq3 cluster was used for clustering." > {log}
            isoseq3 cluster --log-file {log} --use-qvs {input.dataset} {output.clustered} -j {threads} >> {log}
            gzip -d {params.fasta} >> {log}
            """
