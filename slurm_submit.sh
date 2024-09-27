#!/bin/bash
#SBATCH -J snakemake
#SBATCH --time=12:00:00
#SBATCH --partition short
#SBATCH -N 1
#SBATCH -n 2 
#SBATCH --mem=2G
#SBATCH --account MolGen
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aysevilpektas@mbg.au.dk

# Function to display help
show_help() {
    echo "Usage: $0 STEP [-A ACCOUNT_NAME]"
    echo
    echo "Run the TrAnnoScope pipeline with the specified step and optionally with the account."
    echo
    echo "Arguments:"
    echo "  STEP       The step to execute. Available steps are:"
    echo "             qc_rnaseq             - Quality Control for short reads"
    echo "             preprocessing_rnaseq  - Filtering and Trimming of short reads"
    echo "             preprocessing_pacbio  - Processing PacBio long reads"
    echo "             remove_contaminants   - Contamination removal of long reads"
    echo "             error_correction      - Error correction of long reads"
    echo "             classification        - Clustering of long reads"
    echo "             annotation            - Annotation of long reads"
    echo "             quality_assessment    - Quality assessment of the transcriptome"
    echo "             all                   - Run all the steps"
    echo
    echo "Options:"
    echo "  -A ACCOUNT_NAME  Specify the SLURM account name for your project (default: False)"
    echo "  -h, --help  Show this help message and exit"
    echo
    echo "Example:"
    echo "  $0 qc_rnaseq -A <account_name>"
}

declare -A steps
steps["qc_rnaseq"]="rules/quality_SR.smk|log_quality_control_SR.txt|Starting Quality Control for short reads!|Quality control is done! Please check the report and decide whether trimming is needed."
steps["preprocessing_rnaseq"]="rules/preprocessing_SR.smk|log_trim.txt|Starting Filtering and Trimming of short reads!|Filtering and Trimming of short reads are done!"
steps["preprocessing_pacbio"]="rules/preprocessing_LR.smk|log_processing_LR.txt|Starting Processing PacBio long reads!|Processing PacBio long reads is done! Now, you have HQ FL reads."
steps["remove_contaminants"]="rules/contamination.smk|log_remove_contaminants.txt|Starting contamination removal of long reads!|Contamination removal is done! You can proceed to the clustering/classification step."
steps["error_correction"]="rules/error_correction.smk|log_error_correction.txt|Starting error correction of long reads!|Error correction is done! You can proceed to the clustering/classification step."
steps["classification"]="rules/classification.smk|log_clustering.txt|Starting clustering of long reads!|Clustering is done! Now, you have TrAnnoScope results."
steps["annotation"]="rules/annotation.smk|log_annotation.txt|Starting annotation of long reads!|Annotation is done! Now, you have TrAnnoScope results."
steps["quality_assessment"]="rules/quality_assessment.smk|log_quality_assessment.txt|Starting quality assessment for the transcriptome!|Quality assessment is done! Now, you have TrAnnoScope results."

# Define the order of steps
step_order=("qc_rnaseq" "preprocessing_rnaseq" "preprocessing_pacbio" "remove_contaminants" "error_correction" "classification" "annotation" "quality_assessment")


# Check if help is requested
if [[ $1 == "-h" || $1 == "--help" ]]; then
    show_help
    exit 0
fi

# Check if the required arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Error: Missing arguments."
    show_help
    exit 1
fi

STEP=$1
USE_ACCOUNT=false  # Default is not to use the account


# Shift the arguments to process the optional -A flag
shift
while getopts "A" opt; do
    case $opt in
        A) USE_ACCOUNT=true ;;
        *) show_help; exit 1 ;;
    esac
done

run_step() {
    local step=$1
    IFS='|' read -r -a step_details <<< "${steps[$step]}"
    echo "${step_details[2]}"
    if [ "$USE_ACCOUNT" = true ]; then
        snakemake_cmd="snakemake -s ${step_details[0]} --configfile config/config.yaml --use-conda -j 500 \
            --cluster 'sbatch -A {cluster.account_name} --time={cluster.time} -p {cluster.partition} --mem={resources.mem_mb} \
            -c {threads} -o {cluster.output} -e {cluster.error}' --cluster-config config/slurm_config.yaml \
            --latency-wait=60 --reason --show-failed-logs --keep-going --printshellcmds --rerun-incomplete --nolock" 
    else
        snakemake_cmd="snakemake -s ${step_details[0]} --configfile config/config.yaml --use-conda -j 100 \
            --cluster 'sbatch --time={cluster.time} -p {cluster.partition} --mem={resources.mem_mb} \
            -c {threads} -o {cluster.output} -e {cluster.error}' --cluster-config config/slurm_config.yaml \
            --latency-wait=60 --reason --show-failed-logs --keep-going --printshellcmds --rerun-incomplete --nolock"
    fi

    # Run the Snakemake command
    echo "Running: $snakemake_cmd"
    eval $snakemake_cmd
    echo "${step_details[3]}"
}

# Run all steps in the specified order 
if [ "$STEP" == "all" ]; then
    for step in "${step_order[@]}"; do
        run_step $step
    done
else
    run_step $STEP
fi
