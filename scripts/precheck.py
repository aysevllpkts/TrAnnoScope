# Main scripts to carry out the snakemake rules of FLTransAnnot
from snakemake.utils import min_version
import yaml
import os
import sys
import argparse
import subprocess
import logging

# Set minimum Snakemake version
min_version('7.0')

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_config(config_path):
    try:
        with open(config_path, "r") as config_file:
            return yaml.load(config_file, Loader=yaml.FullLoader)
    except FileNotFoundError:
        logging.error(f"Configuration file {config_path} not found.")
        sys.exit(1)
    except yaml.YAMLError as exc:
        logging.error(f"Error parsing YAML configuration: {exc}")
        sys.exit(1)

"""
def install_requirements(rule, log_file, command):
    command = f"nice -5 snakemake -s {rule} --use-conda --conda-create-envs-only"
    try:
        with open(f"logs/{log_file}", "w") as log:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print(output.strip())
                    log.write(output)
            rc = process.poll()
        if rc != 0:
            logging.error(f"Snakemake step {rule} installation failed. Check logs/{log_file} for details.")
            sys.exit(rc)
    except Exception as e:
        logging.error(f"An error occurred while running Snakemake: {e}")
        sys.exit(1)
"""

# Argument parsing
parser = argparse.ArgumentParser(description="Install requirements for FLTransAnnot pipeline")
subparsers = parser.add_subparsers(dest="step", metavar="STEP", help="Which step you want to run")

# Common arguments for all steps
common_parser = argparse.ArgumentParser(add_help=False)
common_parser.add_argument("-c", "--config", default="configs/config.yaml", help="Specify the path for snakemake config file (default: configs/config.yaml)")

# Subparser for each step
subparsers.add_parser("qc_rnaseq", parents=[common_parser], help="Quality Control for Illumina short reads")
subparsers.add_parser("preprocessing_rnaseq", parents=[common_parser], help="Filtering and Trimming of Illumina short reads")
subparsers.add_parser("preprocessing_pacbio", parents=[common_parser], help="Processing PacBio long reads")
subparsers.add_parser("remove_contaminants", parents=[common_parser], help="Contamination removal of PacBio long reads")
subparsers.add_parser("error_correction", parents=[common_parser], help="Error correction of PacBio long reads")
subparsers.add_parser("classification", parents=[common_parser], help="Clustering of PacBio long reads")
subparsers.add_parser("annotation", parents=[common_parser], help="Annotation of PacBio long reads")

args = parser.parse_args()

if not args.step:
    parser.print_help()
    sys.exit(1)

logging.info(f"You will install requirements for {args.step} step")

# Load configuration file
config = load_config(args.config)


if args.step == "qc_rnaseq":
    print("Do you want to install conda requirements for qc_rnaseq?")
    install_input = input()
    if install_input == "y":
        print("Installing package requirements for qc_rnaseq...")
        os.system("nice -5 snakemake -s rules/quality_SR.smk --use-conda --conda-create-envs-only 2>&1 | tee logs/precheck/log_quality_SR_check.txt")
        print("Requirements was installed.")  

if args.step == "preprocessing_rnaseq":
    print("Do you want to install conda requirements for preprocessing_rnaseq?\n if youhave already run qc_rnaseq, no need to install them again!")
    install_input = input()
    if install_input == "y":
        os.system("nice -5 snakemake -s rules/preprocessing_SR.smk --use-conda --conda-create-envs-only 2>&1 | tee logs/precheck/log_preprocessing_SR_check.txt") 
        print("Conda requirements for preprocessing_rnaseq was installed.")
    print("Do you want to download pre-indexed Bowtie2 genomes and configuration file?\n If you already have Genomes folder and conf file, you can skip this step!")
    install_input = input()
    if install_input == "y":
        print("Downloading is starting for FastQScreen Genome Database...")
        os.system("nice -5 fastq_screen --get_genomes --outdir data/fastqscreen 2>>&1 | tee logs/precheck/log_preprocessing_SR_check.txt")
        print("Downloading is done!")

if args.step == "preprocessing_pacbio":
    print("Do you want to install conda requirements for preprocessing_pacbio?")
    install_input = input()
    if install_input == "y":
        print("Installing package requirements for preprocessing_pacbio...")
        os.system("nice -5 snakemake -s rules/preprocessing_LR.smk --use-conda --conda-create-envs-only 2>&1 | tee logs/precheck/log_preprocessing_LR_check.txt") 
        print("Requirements was installed.")

if args.step == "remove_contaminants":
    print("Do you want to install conda requirements for remove_contaminants?")
    install_input = input()
    if install_input == "y":
        print("Installing package requirements for remove_contaminants...")
        os.system("nice -5 snakemake -s rules/remove_contaminants.smk --use-conda --conda-create-envs-only 2>&1 | tee logs/precheck/log_remove_contaminants_check.txt") 
        print("Requirements was installed.")


