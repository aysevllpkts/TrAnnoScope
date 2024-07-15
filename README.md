# TrAnnoScope

**TrAnnoScope** is a comprehensive workflow designed to process PacBio and Illumina reads to generate full-length transcripts and annotate them effectively. 

TrAnnoScope is based on the Snakemake workflow management system. This tool offers a comprehensive suite of features designed to generate the full-length transcriptome from multiple tissues using PacBio and Illumina reads, as well as for functional annotation analysis of the transcriptome. The workflow consists of eight sub-workflows and each part can also be run independently for specific analyses. The initialization part involves the installation of the necessary bioinformatics packages followed by the download and preparation of the relevant reference databases. In the preprocessing step, several quality control tools are used to evaluate and improve the quality of raw Illumina reads. For PacBio reads, isoseq3 is used to obtain high-quality reads. Additionally, optional error correction and contamination-removal steps are defined. FMLRC can be used to increase the quality of long reads obtained by clean Illumina reads, while Blobtools can be used to remove contamination from long reads using coverage information obtained from clean Illumina reads and BUSCO results. During the assembly and reduction step, CDHIT-Est is used to eliminate redundancy between tissue samples, and EvidentialGene is used to remove duplicate reads and obtain protein-coding reads. A quality evaluation of an assembly is conducted based on the completeness, the length distribution, and the ratio of nearly full-length transcripts against databases of well-known biological sequences such as UniProt or closely related species. Trinotate software is used to perform functional annotation of transcripts against various databases as part of the annotation process. Trinotate offers sequence databases, including Pfam, SwissProt, SignalP, TMHMM, EggNog, and Infernal. Additionally, homology searches against the NR and NT databases were added to the workflow as an option. But the users should prepare the necessary files for themselves. As a whole, the workflow includes quality control, reduction of redundancy, obtaining protein-coding transcripts, evaluating assembly, and annotating full-length transcripts and protein sequences.

![image](https://github.com/user-attachments/assets/004eb3db-12a6-4d37-a544-e3aed6ccb8d3)

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

  1. **Clone the repository:**
  ```bash
  git clone https://github.com/aysevllpkts/TrAnnoScope.git
  ```

  2. **Create and activate default environment:**
  ```bash
  conda env creae -n trannoscope -f TrAnnoScope.yaml
  conda activate trannoscope
  ```
  3. **Install necessary packages and databases**
  The initialization part involves the installation of the necessary bioinformatics packages followed by the download and preparation of the relevant reference databases.
  ```bash
  python precheck.py [OPTION STEP]
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
