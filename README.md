# TrAnnoScope

**TrAnnoScope** is a comprehensive workflow designed to process long and short reads to generate full-length transcripts and perform comprehensive functional annotation.

TrAnnoScope is built on the **Snakemake**, offering a robust suite of features for generating full-length transcriptomes from multiple tissues using long and short reads, along with functional annotation of the resulting transcripts. The modular workflow consists of sub-workflows that can be run independently for specific analyses. The initialization step involves installing necessary bioinformatics packages and preparing relevant reference databases.

During preprocessing, various quality control tools are used to assess and improve the quality of raw Illumina reads, while IsoSeq3 is employed to obtain high-quality reads from PacBio data. Optional error correction and contamination removal steps are also defined. FMLRC can be used to enhance the quality of long reads by correcting errors based on clean Illumina reads, while BlobTools2 removes contamination from long reads using coverage information derived from the clean Illumina reads and BUSCO results.

The assembly and reduction step includes CD-HIT-Est to remove redundancy between tissue samples and EvidentialGene to eliminate duplicate reads and extract protein-coding sequences. A quality evaluation of the assembly is performed, assessing completeness, transcript length distribution, and the proportion of nearly full-length transcripts, comparing the results to well-known biological sequence databases such as UniProt or closely related species.

For functional annotation, the workflow employs Trinotate software to annotate transcripts against a range of databases, including Pfam, Rfam (Infernal), SwissProt, SignalP, TMHMM, EggNOG. Additionally, users can opt to include homology searches against the NR and NT databases, with bash scripts provided to help users prepare the necessary files.

In summary, the workflow includes quality control, redundancy reduction, protein-coding transcript extraction, assembly evaluation, and functional annotation of full-length transcripts and protein sequences.

## Workflow
![image](https://github.com/user-attachments/assets/859b75a4-3d46-40fb-bcf0-e9780b49dd7f)



## Features
- **Quality Control**: Perform QC on short reads.
- **Preprocessing**: Filter and trim short reads, preprocess PacBio long reads
- **Contaminants Removal** Remove contamination from long reads (Blobtools2).
- **Error Correction**: Correct errors in long reads by short reads (FMLRC).
- **Classification**: Cluster and Classify long reads.
- **Annotation**: Annotate the full-length transcripts.
- **Quality Assessment**: Assess the quality of the transcriptome.

## Dependencies
- **Tools and Softwares** 
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
  - CD-HIT
  - EvidentialGene
  - DIAMOND
  - Trinotate
  - SQLite
  - HMMER
  - Eggnog-mapper 
  - SignalP
  - TMHMM2
  - NanoPlot
  - R
- **Databases**
  - SwissProt
  - Pfam
  - EGGNOG
  - Rfam
  - NR
  - NT
  - SILVA (LSU, SSU)
  
All conda dependencies can be installed via precheck.py before starting the analysis. Otherwise, they will be installed automatically by TrAnnoScope.

## Installation 

  1. **Clone the repository:**
  ```bash
  git clone https://github.com/aysevllpkts/TrAnnoScope.git
  ```

  2. **Install Miniforge (If not already installed):**
  ```bash
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh
  bash miniforge.sh

  # configuring conda
  conda activate base
  conda config --append channels bioconda
  conda config --set channel_priority disabled # to prevent dependency conflicts
  ```

  3. **Create and activate default environment:**
  ```bash
  conda env create -f TrAnnoScope.yaml
  conda activate trannoscope
  ```

  **Troubleshooting:**
  
  If you still encounter issues while creating or activating the environment using the TrAnnoScope.yaml file, Please follow this step to resolve the problem:
  
  ```bash
  conda env remove -n trannoscope  # Remove existing environment (if any)
  conda env create -n trannoscope -f trannoscope_detailed.yaml  # contains exact versions of all dependencies.
  conda activate trannoscope
  ```

  **System Compatibility:**
  
  Keep in mind that the TrAnnoScope pipeline is developed and optimized for Linux-based systems. If you run the environment on a non-Linux platform, there could be additional compatibility issues. For the best performance and to ensure seamless execution, we recommend using the pipeline in a Linux-based environment.

  
  4. **Modify config.yaml file**

  Edit the config/config.yaml to set file paths and parameters for your analysis before running the installation with precheck.py.

  5. **Install necessary packages and databases**

  The precheck.py script is designed to streamline the setup and installation of conda dependencies required to run the TrAnnoScope pipeline effectively. It ensures that all prerequisites, such as software packages, databases, and configurations, are properly set up before the pipeline execution. The script also allows for modular installation of requirements based on specific pipeline steps.

  **NOTE:** If you want to use SignalP and TmHMM2, first you should install .tar.gz file of them and store it in resources/SignalP and resources/tmhmm2 respectively.
  
  ```bash
  positional arguments:
    STEPS                 Which step do you want to run
      all                 Run full analysis in TrAnnoScope pipeline
      qc_rnaseq           Quality Control for short reads
      preprocessing_rnaseq
                          Filtering and Trimming of short reads
      preprocessing_pacbio
                          Processing PacBio long reads
      remove_contaminants
                          Contamination removal of long reads
      error_correction    Error correction of long reads
      classification      Clustering of long reads
      annotation          Annotation of long reads
      quality_assessment  Quality assessment of transcriptome

  Usage: python precheck.py [STEPS] -c config/config.yaml
  ```
  ### Full Installation
  This command will check and install all necessary tools for the entire pipeline.
  ```bash
  python precheck.py all -c config/config.yaml 
  ```

  ### Step-specific Installation
  You can specify individual steps if you only need to install tools related to a specific part of the pipeline. Replace <STEP> with the desired step name.
  ```bash
  python precheck.py <STEP> -c config/config.yaml
  ```

- **qc_rnaseq**: 
  - No specific actions listed.
  
- **preprocessing_rnaseq**:
  - FastQ Screen Genomes: Download pre-indexed Bowtie2 genomes and configuration files.
  
- **preprocessing_pacbio**:
  - No specific actions listed.
  
- **remove_contaminants**:
  - BUSCO Dataset: Download lineage-specific dataset for transcriptome quality assessment.
  - Download taxonomic information (taxdump).
  
- **error_correction**:
  - No specific actions listed.
  
- **classification**:
  - No specific actions listed.
  
- **annotation**:
  - Install the Trinotate suite and related databases.
  - Create SQLITE.db file and configure the TRINOTATE_DATA_DIR and EGGNOG_DATA_DIR directories.
  - Install SignalP and TmHMM2. (tar.gz files should be located at resources beforehand).
  
- **quality_assessment**:
  - BUSCO Dataset: Download lineage-specific dataset for transcriptome quality assessment.
  
NOTE: Installing FastQ Screen databases and setting up the TRINOTATE_DATA_DIR can be time-consuming due to the large size of databases.

NOTE: When running precheck.py, all conda environments will be installed.


### Additional Configuration for Specific Steps

**For preprocessing_rnaseq step:** 

This step will prompt you to download the FastQScreen genome indexes.

You can modify the fastq_screen.conf file located in resources/fastqscreen to specify which databases to use for removing contaminants and unwanted fragments from the short reads; otherwise, it will default to using all the databases.

```bash
python precheck.py preprocessing_rnaseq -c config/config.yaml
```

Additionally, you can download the SILVA LSU, SSU, and mitochondrial (MT) sequences for your organism of interest.
- SILVA Database:
Download the LSU and SSU reference files from [SILVA Database Exports](https://www.arb-silva.de/no_cache/download/archive/current/Exports)
  -  LSURef: **SILVA_138.1_LSURef_tax_silva.fasta.gz**
  -  SSURef: **SILVA_138.1_SSURef_tax_silva.fasta.gz**

  These are the versions used in the analysis.

- For the MT database, check your organism’s mitochondrial genome on [NCBI](https://www.ncbi.nlm.nih.gov/)

To create Bowtie2 indexes for these databases, use:
```bash
bowtie2-build --threads <THREADS> <INPUT> <OUTPUT>
```
Then, add the path of the Bowtie2 indexes to the resources/fastqscreen/fastq_screen.conf file.

**For the remove_contaminants step:**

This step will prompt you to download the necessary lineage files for BUSCO and taxdump (used by Blobtools2).

```bash
python precheck.py remove_contaminants -c config/config.yaml
```

**For annotation step:** 

This step will prompt you to install Trinotate_v4.2, download the TRINOTATE_DATA_DIR databases, and configure SignalP and TmHMM2.

To use SignalP and TmHMM2, you first need to request the necessary files:
-  Request SignalP 6.0h-fast from:
    -  [SignalP 6.0h-fast Request](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast)
  
    Store the file in: resources/signalp6

-  Request TmHMM2 from:
    - [TMHMM 2.0c Request](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=tmhmm&version=2.0c&packageversion=2.0c&platform=Linux)
  
    Store the file in: resources/tmhmm2

Then run the following command to complete the annotation step:
```bash
python precheck.py annotation -c config/config.yaml
```

**For quality_assessment step:** 

This step will prompt you to download the lineage file for BUSCO. If you have already performed contamination removal, you can answer “NO” to skip this download.

```bash
python precheck.py quality_assessment -c config/config.yaml
```

## Usage
**Running on the local computer (Linux environment):**

Configure settings in config/config.yaml.

    ```bash
    STEP                  Which steps you want to run
      qc_rnaseq           Quality Control for short reads
      preprocessing_rnaseq
                          Filtering and Trimming of short reads
      preprocessing_pacbio
                          Processing PacBio long reads
      remove_contaminants
                          Contamination removal of long reads
      error_correction    Error correction of long reads
      classification      Clustering of long reads
      annotation          Annotation of long reads
      quality_assessment  Quality assessment of the transcriptome for nucleotide and protein sequences
      all                 Run all steps  
    
    Usage: python run_TrAnnoScope.py STEP -c config/config.yaml -t CORES
    Example: python run_TrAnnoScope.py all -c config/config.yaml -t 8
    ```
**Running on SLURM cluster**

To submit a job to the SLURM cluster, configure settings in the config/slurm_config.yaml file.

    ```bash
    Arguments:
      STEP       The step to execute. Available steps are:
                 qc_rnaseq             - Quality Control for short reads
                 preprocessing_rnaseq  - Filtering and Trimming of short reads
                 preprocessing_pacbio  - Processing PacBio long reads
                 remove_contaminants   - Contamination removal of long reads
                 error_correction      - Error correction of long reads
                 classification        - Clustering of long reads
                 annotation            - Annotation of long reads
                 quality_assessment    - Quality assessment of the transcriptome
                 all                   - Run all the steps

    Usage: sbatch slurm_submit.sh STEP [-A]
    Example: sbatch slurm_submit.sh qc_rnaseq -A 
    ```
NOTE: The -A parameter is a boolean (true/false). If used, you must specify a project or account name.

**SLURM Job Configuration**

In the slurm_submit.sh script, you need to specify your SLURM account and email for job notifications. Here’s an example:
    
  ```bash
  #SBATCH --account my_project
  #SBATCH --mail-user=yourname@domain.com
  ```

In the slurm_config.yaml file, ensure that an account_name field is added to specify the project name.

### Test Data

**Verifying the Workflow**

To verify the workflow, you can use the provided test data located in the `data/test_data/` directory.

**1. Configure Test Settings**

First, configure the test settings in `config/test_config.yaml`. You might need to adjust memory and threads. 

To successfully test the `remove_contaminants` step, you need to add the NT database. If the NT database is not added, you can still test specific steps without this requirement.

Once the configuration is complete, you can run the workflow on the test data to ensure everything is set up correctly.

**2. Run the workflow with the test data:**

**Activate environment and prepare necessary files**
```bash
# Activate conda environment
conda activate trannoscope

# Prepare the necessary files
python precheck.py all -c config/test_config.yaml
```

**Run all steps**
```bash
# If you want to run the TrAnnoScope on a local computer for all steps
python run_TrAnnoScope.py all -c config/test_config.yaml

# If you want to run TrAnnoScope slurm cluster for all steps
sbatch slurm_submit.sh all -c config/test_config.yaml [-A]
```

**Run specific steps**

#### Data Paths and Configuration for TrAnnoScope Pipeline

To run the TrAnnoScope pipeline, users need to define specific directories and files for various processing steps in the `config.yaml` file. Below is a description of each data path and the corresponding step:

**Error Correction Step**
- **clean_short_reads**: Directory containing clean short reads (e.g., `<sample>_<fr>.*.fq.gz`).
- **clean_long_reads**: Directory containing clean long reads (e.g., `<sample>.*.fasta`).

> **Note**: Users should define both `clean_short_reads` and `clean_long_reads` for the error correction step.

**Remove Contamination**
- **preprocessed_long_reads**: Directory containing preprocessed long reads (e.g., `<sample>.*.fasta`).
- **clean_short_reads**: Directory containing clean short reads (e.g., `<sample>_<fr>.*.fq.gz`).

> **Note**: For contamination removal, define `preprocessed_long_reads` and `clean_short_reads`.

**Classification**
- **input_for_classification**: Directory containing input files for classification (e.g., `<sample>.*.fasta`).

> **Note**: Users should define `input_for_classification` for the classification step.

**Annotation**
- **input_for_annotation.nucl**: Path to the nucleotide file for annotation (e.g., `transcriptome.fasta`).
- **input_for_annotation.prot**: Path to the protein file for annotation (e.g., `transcriptome_protein.fasta`).

> **Note**: Users should define `input_for_annotation` paths (both nucleotide and protein files) for the annotation step.

**Quality Assessment**
- **input_for_quality.nucl**: Path to the nucleotide file for quality assessment (e.g., `transcriptome.fasta`).
- **input_for_quality.prot**: Path to the protein file for quality assessment (e.g., `transcriptome_protein.fasta`).

> **Note**: Users should define `input_for_quality` paths (both nucleotide and protein files) for quality assessment.

```bash
# If you want to run the TrAnnoScope on a local computer for all steps
python run_TrAnnoScope.py <STEP> -c config/test_config.yaml

# If you want to run TrAnnoScope slurm cluster for all steps
sbatch slurm_submit.sh <STEP> -c config/test_config.yaml [-A]
```

## Citation
The TrAnnoScope paper is currently under review. You can access the preprint version at https://www.preprints.org/manuscript/202411.0489/v1.

## License
TrAnnoScope is licensed under the [MIT License](LICENSE).

## Contact
For any questions or issues, please contact [aysevilpektas@mbg.au.dk](mailto:aysevilpektas@mbg.au.dk).
