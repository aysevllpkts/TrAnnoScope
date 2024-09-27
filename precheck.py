# A python script for managing the pre-requisites and steps involved in running TrAnnoScope pipeline.

from snakemake.utils import min_version
import yaml
import os
import sys
import argparse
import subprocess
import logging

# Set minimum Snakemake version
min_version('6.0')

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_config(config_path):
    """Loads the Snakemake configuration file specified by the --config argument."""
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
    """Creates the directory if it does not exist."""
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            logging.info(f"Created directory: {directory}")
        except Exception as e:
            logging.error(f"Failed to create directory {directory}: {e}")
            sys.exit(1)

def run_command(command, log_file):
    """Executes a system command and logs the output and errors."""
    try:
        with open(log_file, "w") as log:
            logging.info(f"Running command: {command}")
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

def run_snakemake(rule, config, log_file, cores):
    """Runs a Snakemake rule, logging the output and errors to a file."""
    command = f"nice -5 snakemake -s {rule} --configfile {config} --latency-wait 30 --cores {cores} --use-conda --nolock"
    run_command(command, log_file)

def prompt_user(question):
    """Prompts the user with a question and returns their response."""
    logging.info(question)
    return input().strip().lower()

def install_requirements(log_file):
    """Checks if the user wants to install conda requirements for TrAnnoScope. It will install all the conda requirements if the user chooses to."""
    if prompt_user(f"Do you want to install conda requirements for TrAnnoScope? (y/n)") == 'y':
        logging.info(f"Installing package requirements...")
        create_directory_if_not_exists(os.path.dirname(log_file))
        run_command(f"nice -5 snakemake -s rules/conda_env.smk --use-conda --nolock --conda-create-envs-only --cores 1", log_file)
        logging.info("Requirements installed.")

def download_fastqscreen_genomes():
    """Downloads pre-indexed Bowtie2 genomes and the configuration file if the user chooses to do so."""
    if prompt_user("Do you want to download pre-indexed Bowtie2 genomes and configuration file for FastQScreen? If you already have Genomes folder and conf file, you can skip this step! (y/n)") == 'y':
        logging.info("Downloading FastQScreen Genome Database...")
        create_directory_if_not_exists("logs/precheck")
        create_directory_if_not_exists("resources/fastqscreen")
        run_command("nice -5 fastq_screen --get_genomes --outdir resources/fastqscreen", "logs/precheck/log_download_genomes.txt")
        run_command("nice -5 cp resources/fastqscreen/FastQ_Screen_Genomes/fastq_screen.conf resources/fastqscreen/", "logs/precheck/log_copied_fastqscreen_conf.txt")
        logging.info("Downloading completed.")

def download_busco_dataset():
    """Downloads BUSCO databases if the user chooses to do so."""
    if prompt_user("Do you want to download BUSCO dataset? If you already have busco_downloads folder, you can skip this step! (y/n)") == "y":
        logging.info("Downloading BUSCO Dataset...")
        create_directory_if_not_exists("logs/precheck")
        #create_directory_if_not_exists("resources/busco_downloads")
        run_command("busco --list-datasets", "logs/precheck/log_busco_list_datasets.txt")
        dataset_input = prompt_user("Which dataset do you want to download?")
        run_command(f"cd resources/", "logs/precheck/log_copied_busco_downloads.txt")
        run_command(f"cd resources/; busco --download {dataset_input}", "logs/precheck/log_busco_download.txt")
        logging.info("BUSCO Dataset is downloaded.")


def main():
    parser = argparse.ArgumentParser(description="Install requirements for TrAnnoScope pipeline")
    subparsers = parser.add_subparsers(dest="step", metavar="STEPS", help="Which step do you want to run")

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("-c", "--config", required=True, help="Specify the path for snakemake config file (default: configs/config.yaml)")

    subparsers.add_parser("all", parents=[common_parser], help="Run full analysis in TrAnnoScope pipeline")
    subparsers.add_parser("qc_rnaseq", parents=[common_parser], help="Quality Control for short reads")
    subparsers.add_parser("preprocessing_rnaseq", parents=[common_parser], help="Filtering and Trimming of short reads")
    subparsers.add_parser("preprocessing_pacbio", parents=[common_parser], help="Processing PacBio long reads")
    subparsers.add_parser("remove_contaminants", parents=[common_parser], help="Contamination removal of long reads")
    subparsers.add_parser("error_correction", parents=[common_parser], help="Error correction of long reads")
    subparsers.add_parser("classification", parents=[common_parser], help="Clustering of long reads")
    subparsers.add_parser("annotation", parents=[common_parser], help="Annotation of long reads")
    subparsers.add_parser("quality_assessment", parents=[common_parser], help="Quality assessment of transcriptome")

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
        "quality_assessment": ("rules/quality_assessment.smk", "logs/precheck/log_quality_assessment_check.txt")
    }

    # Execute the steps
    if args.step in steps:
        log_file = "logs/precheck/log_conda_env.txt"
        install_requirements(log_file)

        # Additional check for FastQScreen Genome download for RNAseq preprocessing step.
        if args.step == "preprocessing_rnaseq":
            download_fastqscreen_genomes()

        # Download BUSCO Dataset
        if args.step == "remove_contaminants" or args.step == "quality_assessment":
            rule, log_file = steps[args.step]
            download_busco_dataset()

        # Additional check for Contamination Removal step. 
        if args.step == "remove_contaminants":
            if prompt_user(f"Do you want to install taxdump for blobtools in {args.step}? (y/n)") == "y":
                logging.info("Installing taxdump...")
                create_directory_if_not_exists("resources/taxdump")
                #run_command("nice -5 mkdir resources/taxdump", "logs/precheck/log_download_taxdump.txt")
                run_command("nice -5 wget -qO- https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf - -C resources/taxdump", "logs/precheck/log_download_taxdump.txt")
            else:
                logging.info("Skipping taxdump installation. You should include the path to config file!")

        # Additional Trinotate and related Database check for annotation step.
        if args.step == "annotation":
            if prompt_user(f"Do you want to install Trinotate-v4.0.2 for {args.step}? (y/n)") == "y":
                if not os.path.isfile("resources/Trinotate-Trinotate-v4.0.2/Trinotate"):
                    logging.info("Installing Trinotate-v4.0.2...")
                    run_command("nice -5 wget -qO- https://github.com/Trinotate/Trinotate/archive/refs/tags/Trinotate-v4.0.2.tar.gz | tar xzf - -C resources", "logs/precheck/log_download_trinotate.txt")
                    logging.info("Trinotate-v4.0.2 installed successfully.")
                else:
                    logging.info("Trinotate-v4.0.2 already exists. Skipping installation.")
                
                if not os.path.isdir("resources/TRINOTATE_DATA_DIR") or not os.path.isdir("resources/TRINOTATE_DATA_DIR/EGGNOG_DATA_DIR"):
                    logging.info("Installing TRINOTATE_DATA_DIR/EGGNOG_DATA_DIR...")
                    run_snakemake("rules/trinotate_init.smk", args.config, "logs/precheck/log_TRINOTATE_DATA_DIR_check.txt", 2)
                    logging.info("TRINOTATE_DATA_DIR installed successfully.")
                else:
                    logging.info("TRINOTATE_DATA_DIR/EGGNOG_DATA_DIR already exists. Skipping installation.")
            else:
                logging.info("Skipping Trinotate installation. You should define Trinotate and TRINOTATE_DATA_DIR folders in config file!")


    elif args.step == "all":
        
        # Load all necessary conda environments
        log_file = "logs/precheck/log_conda_env.txt"
        install_requirements(log_file)

        # Additional check for FastQScreen Genome download for RNAseq preprocessing step.
        download_fastqscreen_genomes() # user input defined in the function

        # Download BUSCO Dataset - user input defined in the function
        download_busco_dataset()

        # Additional check for Contamination Removal step. 
        if prompt_user(f"Do you want to install taxdump for blobtools in {args.step}? (y/n)") == "y":
            logging.info("Installing taxdump...")
            create_directory_if_not_exists("resources/taxdump")
            #run_command("nice -5 mkdir resources/taxdump", "logs/precheck/log_download_taxdump.txt")
            run_command("nice -5 wget -qO- https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf - -C resources/taxdump", "logs/precheck/log_download_taxdump.txt")
        else:
            logging.info("Skipping taxdump installation. You should include the path to config file!")

        # Additional Trinotate and related Database check for annotation step.
        if prompt_user(f"Do you want to install Trinotate-v4.0.2 for {args.step}? (y/n)") == "y":
            if not os.path.isfile("resources/Trinotate-Trinotate-v4.0.2/Trinotate"):
                logging.info("Installing Trinotate-v4.0.2...")
                run_command("nice -5 wget -qO- https://github.com/Trinotate/Trinotate/archive/refs/tags/Trinotate-v4.0.2.tar.gz | tar xzf - -C resources", "logs/precheck/log_download_trinotate.txt")
                logging.info("Trinotate-v4.0.2 installed successfully.")
            else:
                logging.info("Trinotate-v4.0.2 already exists. Skipping installation.")
            
            if not os.path.isdir("resources/TRINOTATE_DATA_DIR") or not os.path.isdir("resources/TRINOTATE_DATA_DIR/EGGNOG_DATA_DIR"):
                logging.info("Installing TRINOTATE_DATA_DIR/EGGNOG_DATA_DIR...")
                run_snakemake("rules/trinotate_init.smk", args.config, "logs/precheck/log_TRINOTATE_DATA_DIR_check.txt", 2)
                logging.info("TRINOTATE_DATA_DIR installed successfully.")
            else:
                logging.info("TRINOTATE_DATA_DIR/EGGNOG_DATA_DIR already exists. Skipping installation.")
        else:
            logging.info("Skipping Trinotate installation. You should define Trinotate and TRINOTATE_DATA_DIR folders in config file!")

    else:
        logging.error("Invalid step specified. Please choose a valid step.")
        sys.exit(1)

if __name__ == "__main__":
    main()

