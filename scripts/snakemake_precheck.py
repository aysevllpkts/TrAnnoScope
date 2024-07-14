# Main scripts to carry out the snakemake rules of FLTransAnnot
from snakemake.utils import min_version
# Set minimum Snakemake version
min_version('7.0')

import yaml
import os
import time

with open("configs/config.yaml", "r") as config_file:
    config = yaml.load(config_file , Loader=yaml.FullLoader)

# Parameters to control pipeline
project = config["PROJECT"]

## Do you need to do quality control?
qc = config["illumina_QC"]
print("Is quality control conda packages required?\n If you run this before, no need to run again", qc)

if qc:
    # Double check that the user really wants to do QC instead of forgetting to change the param after doing QC
    print("Are you sure to install conda packages for illumina short reads?\n If yes, type 'y'; if not, type 'n' and set 'QC' to 'no' in the config file")
    qc_2nd = input()
    if qc_2nd == "y":
        print("Start install conda environment for Illumina short reads snakemake run!")
        os.system("nice -5 snakemake -s rules/quality_SR.smk --use-conda --conda-create-envs-only 2>&1 | tee logs/illumina_qc/log_quality_control_SR_check.txt")
        print("Installing QC packages for Illumina short reads is done!")
        #os._exit(0)
    else: 
        print("Skipping installation of QC conda environment for quality check Illumina reads!")

trim = config["illumina_TRIMMED"]

if trim:
    if qc == "no":
       os.system("nice -5 snakemake -s rules/preprocessing_SR.smk --use-conda --conda-create-envs-only 2>&1 | tee logs/illumina_qc/log_preprocessing_SR_check.txt") 
    else:
        print("No need to install additional conda packages, there were installed in QC step\n For fastq_screen, you need genomes and indexes for contamination check and conf file")
        print("Do you want to download pre-indexed Bowtie2 genomes and configuration file?\n If you already have Genomes folder and conf file, you can skip this step!")
        install_input = input()
        if install_input == "y":
            print("Downloading is starting...")
            os.system("nice -5 fastq_screen --get_genomes --outdir configs/")
            print("Downloading is done!")
            #os._exit(0)
        else:
            print("No need to install qc.yaml environment")

print("You can add custom pre-indexed Bowtie2 genomes for your organims. Specifically for the removal of mitocondrial fragments")
print("You have to download the genome and index it with bowtie2 manually and add it to the fastq-screen CONF file!. Look at INSTRUCTIONS.txt for more details.)

process_pacbio = config["pacbio_PROCESSING"]

if process_pacbio:
    print("Installing pacbio.yaml for preprocessing PacBio reads")
    os.system("nice -5 snakemake -s rules/preprocessing_LR.smk --use-conda --conda-create-envs-only 2>&1 | tee logs/pacbio/log_preprocessing_LR_check.txt")
else:
    print("No need to install qc.yaml environment")

