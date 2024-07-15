# TrAnnoScope

**TrAnnoScope** is a comprehensive workflow designed to process PacBio and Illumina reads to generate full-length transcripts and annotate them effectively. 

TrAnnoScope is based on the Snakemake workflow management system. This tool offers a comprehensive suite of features designed to ensure high-quality transcriptome data processing. It performs quality control on Illumina short reads, ensuring that only the best data is used in subsequent analyses. The preprocessing step filters and trims Illumina reads while also handling PacBio long reads, preparing them for further analysis. It effectively removes any contamination from PacBio long reads, ensuring clean data. Error correction is performed by using Illumina reads to correct errors in PacBio reads, enhancing accuracy. The classification feature allows for the clustering and classification of PacBio reads, facilitating better data organization. Annotation of full-length transcripts provides detailed insights into the genetic information. Finally, the quality assessment feature evaluates the overall quality of the transcriptome, ensuring reliable results.

A feature of TrAnnoScope is the ability to run all steps together to obtain annotated transcriptomes from raw reads and to run different steps independently. It provides the user with the opportunity to manage their workflow. The contamination removal and error-correction step can be skipped if you only have PacBio reads. Likewise, if you already have a transcriptome, you can run only the annotation and quality assessment steps. For this purpose, the user defines their data in the config.yaml file and run the related snakefile in the rules/ directory.

## Features
- **Quality Control**: Performs QC on Illumina short reads.
- **Preprocessing**: Filters and trims Illumina reads, processes PacBio long reads
- **Contaminants Removal** Removes contamination from PacBio long reads.
- **Error Correction**: Corrects errors in PacBio reads by Illumina reads.
- **Classification**: Cluster and Classify PacBio reads.
- **Annotation**: Annotates the full-length transcripts.
- **Quality Assessment**: Assesses the quality of the transcriptome.

## Dependencies
TrAnnoScope will use:
- Snakemake 7.25.1
- Python
- FastQC
- MultiQC
- FastQScreen
- Fastp
- isoseq3
- BUSCO
- Bowtie2
- Blast+
- Blobtools2
- Seqkit
- FMLRC
- CD-HIT-EST
- EvidentialGene
- DIAMOND
- Trinotate
- SignalP
- TMHMM2
- NanoPlot
- R
  
All conda dependencies can be installed via precheck.py before starting the analysis. Otherwise, they will be installed automatically by TrAnnoScope.

## Installation 

  1. ## Clone the repository:
  ```bash
  git clone https://github.com/aysevllpkts/TrAnnoScope.git
  ```

  2. ## Create and activate default environment:
  ```bash
  conda env creae -n trannoscope -f TrAnnoScope.yaml
  conda activate trannoscope
  ```

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

To verify the workflow, you can use the provided test data located in the `data/test_data/` directory.

1. **Run the workflow with the test data**:


## License
TrAnnoScope is licensed under the [MIT License](LICENSE).

## Contributions
Contributions are welcome. Please fork the repository and submit a pull request.

## Contact
For any questions or issues, please contact [aysevilpektas@mbg.au.dk](mailto:aysevilpektas@mbg.au.dk).
