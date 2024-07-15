# Main script to carry out the Snakemake rules of TrAnnoScope
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
    """
    Load the Snakemake configuration from a YAML file.

    Args:
        config_path (str): Path to the configuration file.

    Returns:
        dict: Configuration dictionary.
    """
    try:
        with open(config_path, "r") as config_file:
            return yaml.load(config_file, Loader=yaml.FullLoader)
    except FileNotFoundError:
        logging.error(f"Configuration file {config_path} not found.")
        sys.exit(1)
    except yaml.YAMLError as exc:
        logging.error(f"Error parsing YAML configuration: {exc}")
        sys.exit(1)

def run_snakemake(rule, config, log_file, cores):
    """
    Execute a Snakemake rule.

    Args:
        rule (str): Snakemake rule file.
        config (str): Path to the configuration file for Snakemake.
        log_file (str): Log file name.
        cores (int): Number of cores to use.
    """
    #command = f"nice -5 snakemake -s {rule} --latency-wait 30 --cores {cores} --configfile {config} --use-conda --nolock"
    command = f"nice -5 snakemake -s {rule} --latency-wait 30 --cores {cores} --configfile {config} --use-conda --nolock --rerun-incomplete --reason --show-failed-logs --keep-going --printshellcmds"
    try:
        os.makedirs("logs", exist_ok=True)
        with open(f"logs/{log_file}", "w") as log:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for output in process.stdout:
                if output:
                    print(output.strip())
                    log.write(output)
            process.wait()
        if process.returncode != 0:
            logging.error(f"Snakemake step {rule} failed. Check logs/{log_file} for details.")
            sys.exit(process.returncode)
    except Exception as e:
        logging.error(f"An error occurred while running Snakemake: {e}")
        sys.exit(1)

def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Running TrAnnoScope pipeline")
    parser.add_argument("-l", "--loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level (default: INFO)")

    subparsers = parser.add_subparsers(dest="step", metavar="STEP", help="Which step you want to run")

    # Common arguments for all steps
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("-c", "--config", required=True, help="Specify the path for Snakemake config file")
    common_parser.add_argument("-t", "--cores", type=int, default=1, help="Specify the number of cores you want to use (default: 1)")

    # Subparser for each step
    subparsers.add_parser("qc_rnaseq", parents=[common_parser], help="Quality Control for Illumina short reads")
    subparsers.add_parser("preprocessing_rnaseq", parents=[common_parser], help="Filtering and Trimming of Illumina short reads")
    subparsers.add_parser("preprocessing_pacbio", parents=[common_parser], help="Processing PacBio long reads")
    subparsers.add_parser("remove_contaminants", parents=[common_parser], help="Contamination removal of PacBio long reads")
    subparsers.add_parser("error_correction", parents=[common_parser], help="Error correction of PacBio long reads")
    subparsers.add_parser("classification", parents=[common_parser], help="Clustering of PacBio long reads")
    subparsers.add_parser("annotation", parents=[common_parser], help="Annotation of PacBio long reads")
    subparsers.add_parser("quality_assessment", parents=[common_parser], help="Quality assessment of the transcriptome for nucleotide and protein sequences"),
    subparsers.add_parser("all", parents=[common_parser], help="Run all the steps")

    args = parser.parse_args()

    if not args.step:
        parser.print_help()
        sys.exit(1)

    return args

def main():
    args = parse_args()

    logging.getLogger().setLevel(args.loglevel)

    logging.info(f"You will run {args.step} step")

    # Load configuration file
    config = load_config(args.config)

    # Dictionary to map steps to Snakemake rules and log files
    steps = {
        "qc_rnaseq": ("rules/quality_SR.smk", "log_quality_control_SR.txt", "Start Quality Control for Illumina short reads!", "Quality control is done!\n Please check the report and decide whether trimming is needed."),
        "preprocessing_rnaseq": ("rules/preprocessing_SR.smk", "log_trim.txt", "Start Filtering and Trimming of Illumina short reads!", "Filtering and Trimming of Illumina short reads are done!"),
        "preprocessing_pacbio": ("rules/preprocessing_LR.smk", "log_processing_LR.txt", "Start Processing PacBio long reads!", "Processing PacBio long reads is done!\n Now, you have HQ FL reads."),
        "remove_contaminants": ("rules/contamination.smk", "log_remove_contaminants.txt", "Start contamination removal of PacBio long reads!", "Contamination removal is done!\n You can proceed to the clustering/classification step."),
        "error_correction": ("rules/error_correction.smk", "log_error_correction.txt", "Start error correction of PacBio long reads!", "Error correction is done!\n You can proceed to the clustering/classification step."),
        "classification": ("rules/classification.smk", "log_clustering.txt", "Start clustering of PacBio long reads!", "Clustering is done!\n Now, you have TrAnnoScope results."),
        "annotation": ("rules/annotation.smk", "log_annotation.txt", "Start annotation of PacBio long reads!", "Annotation is done!\n Now, you have TrAnnoScope results."),
        "quality_assessment": ("rules/quality_assessment.smk", "log_quality_assessment.txt", "Start quality assessment for the transcriptome!", "Quality assessment is done!\n Now, you have TrAnnoScope results.")
    }

    # Execute the step
    if args.step in steps:
        rule, log_file, start_message, end_message = steps[args.step]
        logging.info(start_message)
        run_snakemake(rule, args.config, log_file, args.cores)
        logging.info(end_message)
    
    elif args.step == "all":
        for step in step: 
            rule, log_file, start_message, end_message = steps[step]
            logging.info(start_message)
            run_snakemake(rule, log_file, args.jobs, args.account_name, args.partition, args.time)
            logging.info(end_message)
        
        """
        rule, log_file, start_message, end_message = steps["qc_rnaseq"]
        logging.info(start_message)
        run_snakemake("rules/quality_SR.smk", args.config, log_file, args.cores)
        logging.info(end_message)

        rule, log_file, start_message, end_message = steps["preprocessing_rnaseq"]
        logging.info(start_message)
        run_snakemake("rules/preprocessing_SR.smk", args.config, log_file, args.cores)
        logging.info(end_message)

        rule, log_file, start_message, end_message = steps["preprocessing_pacbio"]
        logging.info(start_message)
        run_snakemake("rules/preprocessing_LR.smk", args.config, log_file, args.cores)
        logging.info(end_message)

        rule, log_file, start_message, end_message = steps["remove_contaminants"]
        logging.info(start_message)
        run_snakemake("rules/contamination.smk", args.config, log_file, args.cores)
        logging.info(end_message)

        rule, log_file, start_message, end_message = steps["error_correction"]
        logging.info(start_message)
        run_snakemake("rules/error_correction.smk", args.config, log_file, args.cores)
        logging.info(end_message)

        rule, log_file, start_message, end_message = steps["classification"]
        logging.info(start_message)
        run_snakemake("rules/classification.smk", args.config, log_file, args.cores)
        logging.info(end_message)

        rule, log_file, start_message, end_message = steps["annotation"]
        logging.info(start_message)
        run_snakemake("rules/annotation.smk", args.config, log_file, args.cores)
        logging.info(end_message)

        rule, log_file, start_message, end_message = steps["quality_assessment"]
        logging.info(start_message)
        run_snakemake("rules/quality_assessment.smk", args.config, log_file, args.cores)
        logging.info(end_message)
        """

    else:
        logging.error("Invalid step specified. Please choose a valid step.")
        sys.exit(1)

    
    

if __name__ == "__main__":
    main()

