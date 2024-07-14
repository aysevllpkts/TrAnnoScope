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

def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            logging.info(f"Created directory: {directory}")
        except Exception as e:
            logging.error(f"Failed to create directory {directory}: {e}")
            sys.exit(1)

def run_command(command, log_file):
    try:
        with open(log_file, "w") as log:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    logging.info(output.strip())
                    log.write(output)
            rc = process.poll()
        if rc != 0:
            logging.error(f"Command failed. Check {log_file} for details.")
            sys.exit(rc)
    except Exception as e:
        logging.error(f"An error occurred while running command: {e}")
        sys.exit(1)

def run_snakemake(rule, log_file, cores):
    command = f"nice -5 snakemake -s {rule} --latency-wait 30 --cores {cores} --use-conda"
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
            logging.error(f"Snakemake step {rule} failed. Check logs/{log_file} for details.")
            sys.exit(rc)
    except Exception as e:
        logging.error(f"An error occurred while running Snakemake: {e}")
        sys.exit(1)

def install_requirements(step, rule, log_file):
    logging.info(f"Do you want to install conda requirements for {step}? (y/n)")
    install_input = input().strip().lower()
    if install_input == 'y':
        logging.info(f"Installing package requirements for {step}...")
        create_directory_if_not_exists(os.path.dirname(log_file))
        run_command(f"nice -5 snakemake -s {rule} --use-conda --conda-create-envs-only", log_file)
        logging.info("Requirements installed.")

def download_genomes():
    logging.info("Do you want to download pre-indexed Bowtie2 genomes and configuration file? If you already have Genomes folder and conf file, you can skip this step! (y/n)")
    install_input = input().strip().lower()
    if install_input == 'y':
        logging.info("Downloading FastQScreen Genome Database...")
        create_directory_if_not_exists("logs/precheck")
        run_command("nice -5 fastq_screen --get_genomes --outdir data/fastqscreen", "logs/precheck/log_download_genomes.txt")
        logging.info("Downloading completed.")

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

# Dictionary to map steps to Snakemake rules and log files
steps = {
    "qc_rnaseq": ("rules/quality_SR.smk", "logs/precheck/log_quality_SR_check.txt"),
    "preprocessing_rnaseq": ("rules/preprocessing_SR.smk", "logs/precheck/log_preprocessing_SR_check.txt"),
    "preprocessing_pacbio": ("rules/preprocessing_LR.smk", "logs/precheck/log_preprocessing_LR_check.txt"),
    "remove_contaminants": ("rules/contamination.smk", "logs/precheck/log_remove_contaminants_check.txt"),
    "error_correction": ("rules/error_correction.smk", "logs/precheck/log_error_correction_check.txt"),
    "classification": ("rules/classification.smk", "logs/precheck/log_clustering_check.txt"),
    "annotation": ("rules/annotation.smk", "logs/precheck/log_annotation_check.txt"),
}

# Execute the step
if args.step in steps:
    rule, log_file = steps[args.step]
    install_requirements(args.step, rule, log_file)
    
    # Additional genome download step for preprocessing_rnaseq
    if args.step == "preprocessing_rnaseq":
        download_genomes()
    # Additional installations for contamination removal
    if args.step == "remove_contaminants":
        # install taxdump
        logging.info(f"Do you want to install taxdump for blobtools in {args.step}? (y/n)")
        install_input = input().strip().lower()        
        if install_input == "y":
            print("Installation is starting...")
            create_directory_if_not_exists(os.path.dirname("configs/taxdump"))   
            run_command("nice -5 mkdir configs/taxdump", "logs/precheck/log_download_taxdump.txt")
            run_command("nice -5 wget -qO- https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf - -C configs/taxdump" , "logs/precheck/log_download_taxdump.txt")
        else:
            print("Skip installing taxdump, you should include path to config file!")
        # do not download busco data, change code in snakefile
    
    if args.step == "annotation":
        logging.info(f"Do you want to install Trinotate-v4.0.2 for {args.step}? (y/n)")
        install_input = input().strip().lower()

        if install_input == "y":
            if os.path.isfile("configs/Trinotate-Trinotate-v4.0.2/Trinotate"):
                print(f"Trinotate-v4.0.2 already exists. Skipping installation.")  

            else:
                print("Installation is starting for Trinotate-v4.0.2...")
                run_command("nice -5 wget -qO- https://github.com/Trinotate/Trinotate/archive/refs/tags/Trinotate-v4.0.2.tar.gz | tar xzf - -C configs", "logs/precheck/log_download_trinotate.txt")
                print("Trinotate-v4-0-2 was installed successfully to configs folder!")
            
            if os.path.isfile("data/TRINOTATE_DATA_DIR"):
                print(f"TRINOTATE_DATA_DIR folder already exists. Skipping installation database.")
            
            else:
                print("TRINOTATE_DATA_DIR installation is starting... YOU HAVE TO INSTALL TRINOTATE FIRST and required conda packages for snakemake!")
                #install_requirements("annotation", "rules/trinity_init.smk", "logs/precheck/log_install_trinotate.txt")
                #run_snakemake("rules/trinity_init.smk", "logs/precheck/log_install_trinotate.txt", 1)
                os.system("nice -5 snakemake -s rules/trinotate_init.smk --use-conda --cores 4 2>&1 | tee logs/precheck/log_remove_contaminants_check.txt") 
                print("TRINOTATE_DATA_DIR was installed successfully to data folder!")

        else:   
            print("Skipping installation. You should define Trinotate and TRINOTATE_DATA_DIR folders in config file!")
     
     if args.step == "annotation":
        logging.info(f"Activating signalp6")
        run_command("nice -5")
            
else:
    logging.error("Invalid step specified. Please choose a valid step.")
    sys.exit(1)
