# TrAnnoScope

**TrAnnoScope** is a comprehensive workflow designed to process PacBio and Illumina reads to generate full-length transcripts and annotate them effectively.

## Features
- **Quality Control**: Performs QC on Illumina short reads.
- **Preprocessing**: Filters and trims Illumina reads, processes PacBio long reads
- **Contaminants Removal** Removes contamination from PacBio long reads.
- **Error Correction**: Corrects errors in PacBio reads by Illumina reads.
- **Classification**: Cluster and Classify PacBio reads.
- **Annotation**: Annotates the full-length transcripts.
- **Quality Assessment**: Assesses the quality of the transcriptome.

## Usage

### Step-by-Step Guide

1. **Quality Control for Illumina Reads**:
    ```bash
    ./slurm_submit.sh qc_rnaseq [-A <account_name>]
    ```

2. **Filtering and Trimming Illumina Reads**:
    ```bash
    ./slurm_submit.sh preprocessing_rnaseq [-A <account_name>]
    ```

3. **Processing PacBio Reads**:
    ```bash
    ./slurm_submit.sh preprocessing_pacbio [-A <account_name>]
    ```

4. **Removing Contaminants**:
    ```bash
    ./slurm_submit.sh remove_contaminants [-A <account_name>]
    ```

5. **Error Correction**:
    ```bash
    ./slurm_submit.sh error_correction [-A <account_name>]
    ```

6. **Clustering PacBio Reads**:
    ```bash
    ./slurm_submit.sh classification [-A <account_name>]
    ```

7. **Annotating PacBio Reads**:
    ```bash
    ./slurm_submit.sh annotation [-A <account_name>]
    ```

8. **Quality Assessment**:
    ```bash
    ./slurm_submit.sh quality_assessment [-A <account_name>]
    ```

9. **Running All Steps**:
    ```bash
    ./slurm_submit.sh all -A [-A <account_name>]
    ```

### Test Data

To verify the workflow, you can use the provided test data located in the `test_data/` directory.

1. **Run the workflow with the test data**:


## License
TrAnnoScope is licensed under the [MIT License](LICENSE).

## Contributions
Contributions are welcome. Please fork the repository and submit a pull request.

## Contact
For any questions or issues, please contact [aysevilpektas@mbg.au.dk](mailto:aysevilpektas@mbg.au.dk).
