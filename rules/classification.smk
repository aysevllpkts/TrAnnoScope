# Aysevil Pektas
# 24-03-2024
# CD-HIT-EST

import glob
import os

configfile: "config/config.yaml"

##------------------------- Define your input file ---------------------------##

# Define output folder:
OUTDIR = "output_" + config["PROJECT"] 

# if you already have preprocessed (High Quality) pacbio reads
if config["input_for_classification"] != "":
    if os.path.isdir(config["input_for_classification"]):
        print("Your LR input defined in the config file will be used for classification step.")
        print("Only your desired input files in the directory should contain .fasta extension")
        ISOSEQ_DIR = config['input_for_classification']                                          # path
        fasta_files = glob.glob(os.path.join(config["input_for_classification"] , "*.fasta"))    # list of input files
        LR_SAMPLES = [os.path.basename(path).split('.')[0] for path in fasta_files]
        #LR_SAMPLES = get_samples(config['input_for_classification'], "*.fasta")                  # wildcards
    else:
        print(f"Error: The path for LR provided in the config file does not exist: {config['clean_long_reads']}")
        LR_SAMPLES = None  # Or handle the error accordingly
elif config["input_for_classification"] == "" and os.path.isdir(OUTDIR + "/pacbio/error_correction"):
    print("LR input files created by error correction step will be used for classification step.")
    ISOSEQ_DIR = OUTDIR + "/pacbio/error_correction"                                                                                              # path
    LR_SAMPLES, = glob_wildcards(OUTDIR + "/pacbio/error_correction/{sample}.clean.corrected.fasta")                           # wildcards
    fasta_files = [os.path.join(OUTDIR, f"pacbio/error_correction/{sample}.clean.corrected.fasta") for sample in LR_SAMPLES]   # list of input files
else:
    print("You don't have a input file. Define the path of your input directory correctly or run the previous steps of the workflow.")
    LR_SAMPLES = None 

# rest of the suffix of the input file.
LR_file_suffix = '.'.join(os.path.basename(fasta_files[0]).split('.')[1:]) # Get file name except wildcard part

# Handling the else case: Ensure that Snakemake will stop execution if the transcriptome is not defined
if LR_SAMPLES is None:
    raise ValueError("Input files are not correctly defined. Please provide the path in the config file or ensure the previous steps are completed.")
else:
    print(f"Input LR directory: {ISOSEQ_DIR}")
    print(f"Found LR samples: {LR_SAMPLES}")
    print(f"Input LR files: {fasta_files}")
    print(f"LR file suffix: {LR_file_suffix}")

##------------------------------------------------------------------------------##

optional_input = [] 
if config["cd-hit-est"]["proceed"] == "yes":
    optional_input.append(expand(OUTDIR + "/pacbio/classification/cdhit/{sample}_clean_corrected_cdhit.fasta", sample = LR_SAMPLES))
    optional_input.append(OUTDIR + "/pacbio/classification/evigene/merged_clean_corrected.fasta")
else: 
    optional_input.append(OUTDIR + "/pacbio/classification/evigene/merged_clean_corrected.fasta")
#print(optional_input)

rule all:   
    input:
        optional_input,
        #"preprocessing/pacbio/evigene/okayset/merged_clean_corrected_cdhit.okay.mrna"
        OUTDIR + "/pacbio/classification/evigene/evigene.done"

if config["cd-hit-est"]["proceed"] == "yes":
    rule cd_hit_est:
        input:  
            ISOSEQ_DIR + "/{sample}." + LR_file_suffix
        output:
            OUTDIR + "/pacbio/classification/cdhit/{sample}_clean_corrected_cdhit.fasta"
        log:
            "logs/classification/cd_hit_est.{sample}.log"
        conda: 
            "../envs/cd-hit-est.yaml"
        params:
            similarity = config["cd-hit-est"]["similarity"],
            extra_parameters = config["cd-hit-est"]["extra_parameters"]
        threads:
            config["cd-hit-est"]["threads"]
        resources:
            mem_mb = config["cd-hit-est"]["memory"]
        benchmark:
            "benchmarks/pacbio/classification/cdhitest_{sample}.benchmark.txt"
        shell:
            """
            cd-hit-est -i {input} -o {output} -c {params.similarity} -T {threads} -M {resources.mem_mb} {params.extra_parameters} &> {log}
            """

    rule merge_sample:
        input:
            expand(OUTDIR + "/pacbio/classification/cdhit/{sample}_clean_corrected_cdhit.fasta", sample = LR_SAMPLES)
        output: 
            OUTDIR + "/pacbio/classification/evigene/merged_clean_corrected.fasta" # "." change to "_" in naming, Evigene has problems with "." naming. 
        threads:
            1
        resources:
            mem_mb = 100 
        shell:
            """
            cat {input} >> {output}
            """
else: 
    rule merge_sample:
        input:
            ISOSEQ_DIR + "/{sample}." + LR_file_suffix
        output: 
            OUTDIR + "/pacbio/classification/evigene/merged_clean_corrected.fasta"
        threads:
            1
        resources:
            mem_mb = 100 
        shell:
            """
            cat {input} >> {output}
            """    

rule evigene: 
    input:
        OUTDIR + "/pacbio/classification/evigene/merged_clean_corrected.fasta" 
    output:
        touch(OUTDIR + "/pacbio/classification/evigene/evigene.done")
    log:
        "logs/classification/evigene_tr2aacds4.log"
    conda:
        "../envs/evigene.yaml"
    params:
        fasta = "merged_clean_corrected.fasta",
        minaa = config["evigene"]["min_aa"],
        extra = config["evigene"]["extra"],
        path = OUTDIR + "/pacbio/classification/evigene"
    threads:
        config["evigene"]["threads"]
    resources:
        mem_mb=config["evigene"]["memory"]
    benchmark:
        "benchmarks/pacbio/classification/evigene.benchmark.txt"
    shell: 
        """
        mkdir -p {params.path}
        cd {params.path}
        $EVIGENEHOME/scripts/prot/tr2aacds4.pl -cdnaseq {params.fasta} -NCPU={threads} -MAXMEM={resources.mem_mb} -MINAA={params.minaa} -logfile ../../../../{log} {params.extra} &>> ../../../../{log} 
        """
