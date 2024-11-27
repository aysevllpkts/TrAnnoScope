# Aysevil Pektas
# 23-03-2024
# Error correction by short reads - OPTIONAL rule!

import os
import glob

# Define config file:
configfile: "config/config.yaml"

##------------------------- Define your input file ---------------------------##

# Define output folder:
OUTDIR = "output_" + config["PROJECT"] 

"""
def get_samples(directory, pattern):
    files = glob.glob(os.path.join(directory, pattern))
    return [re.match(r"(.*)\..*\.fasta", os.path.basename(f)).group(1) for f in files if re.match(r"(.*)\..*\.fasta", os.path.basename(f))]
""" 

# Check the 'proceed' variable from config
if config["error_correction"]["proceed"] == "yes":
    print("Proceeding with error correction module...")

    # if you already have preprocessed (High Quality) pacbio reads
    if config["clean_long_reads"] != "":
        if os.path.isdir(config["clean_long_reads"]):
            print("Your LR input defined in the config file will be used for error correction step.")
            print("Only your desired input files in the directory should contain .fasta extension")
            ISOSEQ_DIR = config['clean_long_reads']                                          # path
            fasta_files = glob.glob(os.path.join(config["clean_long_reads"] , "*.fasta"))    # list of input files
            LR_SAMPLES = [os.path.basename(path).split('.')[0] for path in fasta_files]
            #LR_SAMPLES = get_samples(config['clean_long_reads'], "*.fasta")                  # wildcards
            LR_SAMPLES = set(sorted(LR_SAMPLES))
            # rest of the suffix of the input file.
            LR_file_suffix = '.'.join(os.path.basename(fasta_files[0]).split('.')[1:]) # Get file name except wildcard part
        else:
            print(f"Error: The path for LR provided in the config file does not exist: {config['clean_long_reads']}")
            LR_SAMPLES = None  # Or handle the error accordingly
    elif config["clean_long_reads"] == "": #nd os.path.isdir(OUTDIR + "/pacbio/remove_contaminants/clean_reads"):
        print("LR input files created by contamnination removal step will be used for error correction step.")
        ISOSEQ_DIR = OUTDIR + "/pacbio/remove_contaminants/clean_reads"                                                                                              # path
        LR_SAMPLES, = glob_wildcards(OUTDIR + "/pacbio/remove_contaminants/clean_reads/{sample}.clean.fasta")                           # wildcards
        LR_SAMPLES = set(sorted(LR_SAMPLES))
        fasta_files = [os.path.join(OUTDIR, f"pacbio/remove_contaminants/clean_reads/{sample}.clean.fasta") for sample in LR_SAMPLES]   # list of input files
        LR_file_suffix = 'clean.fasta'

    else:
        print("You don't have a input file. Define the path of your input directory correctly or run the previous steps of the workflow.")
        LR_SAMPLES = None 


    # if you already have clean short reads
    if config["clean_short_reads"] != "":
        if os.path.isdir(config["clean_short_reads"]):
            print("\nYour SR input defined in the config file will be used for contamination removal step.") 
            print("Only your desired input files in the directory should contain .fq.gz extension")      
            RNASEQ_DIR = config['clean_short_reads']                                        # path
            fq_files = glob.glob(os.path.join(config["clean_short_reads"] , "*.fq.gz"))     # list of files
            SR_SAMPLES = [os.path.basename(path).split('_')[0] for path in fq_files]        # wildcards SAMPLE 
            SR_SAMPLES = set(sorted(SR_SAMPLES))
            SR_file_suffix = fq_files[0].split('.', 1)[1]                                   # file suffix, everything except <sample>
        else:
            print(f"Error: The path provided in the config file does not exist: {config['clean_short_reads']}")
            SR_SAMPLES = None 
    elif config["clean_short_reads"] == "" and os.path.isdir(OUTDIR + "/illumina/trimmed_fastp"):
        print("SR Input files created by RNASEQ preprocessing step will be used for contamination removal step.")
        RNASEQ_DIR = OUTDIR + "/illumina/trimmed_fastp"                                                      # path
        fq_files = glob.glob(os.path.join(OUTDIR + "/illumina/trimmed_fastp/", "*.fq.gz"))                   # list of files
        SR_SAMPLES, FRR = glob_wildcards(OUTDIR + "/illumina/trimmed_fastp/{sample}_{frr}.trimmed.fq.gz")    # wildcards SAMPLE and FR
        SR_SAMPLES = set(sorted(SR_SAMPLES))
        SR_file_suffix = "trimmed.fq.gz"                                                                     # file suffix, everything except SAMPLE and FR
    else:
        print("You don't have SR input files. Define the path of your input directory correctly or run the previous steps of the workflow.")
        SR_SAMPLES = None 

    # Handling the else case: Ensure that Snakemake will stop execution if the transcriptome is not defined
    if SR_SAMPLES is None or LR_SAMPLES is None:
        print(SR_SAMPLES)
        print(LR_SAMPLES)
        raise ValueError("Input files are not correctly defined. Please provide the path in the config file or ensure the previous steps are completed.")
    elif SR_SAMPLES != LR_SAMPLES:
        print(f"Input LR directory: {ISOSEQ_DIR} and SR directory: {RNASEQ_DIR}")
        print(SR_SAMPLES)
        print(LR_SAMPLES)
        raise ValueError("Wildcards for short are long reads do not overlaps, Please ensure the files have same sample wildcards!")
    else:
        print(f"Input LR directory: {ISOSEQ_DIR} and SR directory: {RNASEQ_DIR}")
        print(f"Found LR samples: {LR_SAMPLES} and SR samples {SR_SAMPLES}")
        print(f"Input LR files: {fasta_files} and SR files {fq_files}")
        print(f"LR file suffix: {LR_file_suffix} and SR file suffix: {SR_file_suffix}")


    ##--------------------------------------------------------------------------------##

    rule all:
        input:
            expand(OUTDIR + "/pacbio/error_correction/{sample}_comp_msbwt.npy", sample = LR_SAMPLES),
            expand(OUTDIR + "/pacbio/error_correction/{sample}.clean.corrected.fasta", sample = LR_SAMPLES),

    rule bwt:
        input: 
            fq_files
        output:
            OUTDIR + "/pacbio/error_correction/{sample}_comp_msbwt.npy"
        log:
            "logs/error_correction/bwt.{sample}.log"
        conda:
            "../envs/fmlrc.yaml"
        params:
            indir = RNASEQ_DIR
        threads:
            config["error_correction"]["threads"]   
        resources:
            mem_mb=config["error_correction"]["memory"]
        benchmark:
            "benchmarks/pacbio/error_correction/bwt_{sample}.benchmark.txt"
        shell:
            """
            gunzip -c {params.indir}/{wildcards.sample}_?.*.fq.gz \
                | awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN \
                | fmlrc-convert {output}
            """

    rule fmlrc:
        input:
            bwt = OUTDIR + "/pacbio/error_correction/{sample}_comp_msbwt.npy",
            long_reads = ISOSEQ_DIR + "/{sample}." + LR_file_suffix
        output:
            OUTDIR + "/pacbio/error_correction/{sample}.clean.corrected.fasta"
        log:
            "logs/error_correction/fmlrc.{sample}.log"
        conda:
            "../envs/fmlrc.yaml"
        params:
            parameters = config["fmlrc_parameters"]
        threads:
            config["error_correction"]["threads"] 
        resources:
            mem_mb=config["error_correction"]["memory"]
        benchmark:
            "benchmarks/pacbio/error_correction/fmlrc_{sample}.benchmark.txt"
        shell:
            """
            fmlrc -p {threads} {params.parameters} {input.bwt} {input.long_reads} {output}
            """
else:
    print("Error correction module is skipped.")
