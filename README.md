# TrAnnoScope

**TrAnnoScope** is a comprehensive workflow designed to process both long and short reads to generate full-length transcripts and perform comprehensive functional annotation.

TrAnnoScope is built on the Snakemake, offering a robust suite of features for generating full-length transcriptomes from multiple tissues using long and short reads, along with functional annotation of the resulting transcripts. The workflow is modular, consisting of sub-workflows that can be run independently for specific analyses. The initialization step involves installing necessary bioinformatics packages and preparing relevant reference databases.

During preprocessing, various quality control tools are used to assess and improve the quality of raw Illumina reads, while IsoSeq3 is employed to obtain high-quality reads from PacBio data. Optional error correction and contamination removal steps are also defined. FMLRC can be used to enhance the quality of long reads by correcting errors based on clean Illumina reads, while BlobTools2 removes contamination from long reads using coverage information derived from the clean Illumina reads and BUSCO results.

The assembly and reduction step includes CD-HIT-Est to remove redundancy between tissue samples and EvidentialGene to eliminate duplicate reads and extract protein-coding sequences. A quality evaluation of the assembly is performed, assessing completeness, transcript length distribution, and the proportion of nearly full-length transcripts, comparing the results to well-known biological sequence databases such as UniProt or closely related species.

For functional annotation, the workflow employs Trinotate software to annotate transcripts against a range of databases, including Pfam, Rfam (Infernal), SwissProt, SignalP, TMHMM, EggNOG. Additionally, users can opt to include homology searches against the NR and NT databases, with bash scripts provided to help users prepare the necessary files.

In summary, the workflow includes quality control, redundancy reduction, protein-coding transcript extraction, assembly evaluation, and functional annotation of full-length transcripts and protein sequences.

## Workflow
![image](https://github.com/user-attachments/assets/859b75a4-3d46-40fb-bcf0-e9780b49dd7f)



## Features
- **Quality Control**: Performs QC on short reads.
- **Preprocessing**: Filters and trims reads, preprocess PacBio long reads
- **Contaminants Removal** Removes contamination from long reads (Blobtools2).
- **Error Correction**: Corrects errors in long reads by short reads (FMLRC).
- **Classification**: Cluster and Classify long reads.
- **Annotation**: Annotates the full-length transcripts.
- **Quality Assessment**: Assesses the quality of the transcriptome.

## Dependencies
- **Programs** 
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
  
All conda dependencies can be installed via precheck.py before starting the analysis. Otherwise, they will be installed automatically by TrAnnoScope.

## Installation 

  1. **Clone the repository:**
  ```bash
  git clone https://github.com/aysevllpkts/TrAnnoScope.git
  ```

  2. **Create and activate default environment:**
  ```bash
  conda env create -f TrAnnoScope.yaml
  conda activate trannoscope
  ```
  3. **Modify config.yaml file**
    
  Before proceeding the installation, modify the config/config.yaml based on your needs.

  4. **Install necessary packages and databases**
  
  The Python script  called precheck.py is designed to manage the prerequisites and steps involved in running the TrAnnoScope pipeline. It utilizes Snakemake for workflow management and integrates various bioinformatics tools and databases. The script ensures that necessary environments and resources are installed and configured, guiding users through each step of the pipeline setup and execution process. 

  **NOTE:** If you want to use SignalP and TmHMM2, first you should install tar.gz file of them and store it in resources/
  
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
**For preprocessing_rnaseq step:** 

it will ask to download FastQScreen Genome indexes. 

```bash
python precheck.py preprocessing_rnaseq -c config/config.yaml
```

Additionally, you can download SILVA LSU and SSU and also MT sequences of your organism of interest.  
You can download SILVA DB from: https://www.arb-silva.de/no_cache/download/archive/current/Exports \ as LSURef: **SILVA_138.1_LSURef_tax_silva.fasta.gz** and SSURef: **SILVA_138.1_SSURef_tax_silva.fasta.gz** (these are the current one in that time).\
\
For MT DB, you can check your organism in NCBI for MT genome.

To create bowtie2 index for DBs
```bash
bowtie2-build --threads <THREADS> <INPUT> <OUTPUT>
```
Add path of bowtie2-indexes to resources/fastqscreen/fastq_screen.conf

**For the remove_contaminants step:**

it will ask to download interested lineage files for BUSCO and taxdump which is used by Blobtools2.

```bash
python precheck.py remove_contaminants -c config/config.yaml
```

**For annotation step:** 

it will ask to download TRINOTATE_DATA_DIR for databases and register SignalP and TmHMM2.

If you want SignalP and TmHMM2 search, first you have to install request files.

Request to get signalp-6h-fast from: https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast 
The file should be stored at resources/signalp6

Request tmhmm2 from: https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=tmhmm&version=2.0c&packageversion=2.0c&platform=Linux 
The file should be stored at resources/tmhmm2

```bash
python precheck.py annotation -c config/config.yaml
```

**For quality_assessment step:** 

it will ask to download the interested lineage file for BUSCO, if you have already performed, contamination removal, you can say "NO".

```bash
python precheck.py quality_assessment -c config/config.yaml
```

**If you run ALL step:** 

it will ask about all the necessary files to download

```bash
python precheck.py all -c config/config.yaml
```

## Usage
**Running on the local computer:**

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
    
    Usage: python run_FLAnnotTrans.py STEP -c config/config.yaml -t CORES
    Example: python run_FLAnnotTrans.py all -c config/config.yaml -t 2
    ```
**Running on SLURM cluster**

Configure settings in config/slurm_config.yaml.

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

    Usage: sbatch slurm_submit.sh STEP [-A <slurm account name>]
    Example: sbatch slurm_submit.sh qc_rnaseq [-A <slurm account_name>]
    ```


### Test Data

To verify the workflow, you can use the provided test data in the `data/test_data/` directory.

1. **Run the workflow with the test data**:
   
```bash
# Activate conda environment
conda activate trannoscope

# Prepare the necessary files
python precheck.py all -c config/test_config.yaml

# If you want to run the TrAnnoScope on a local computer
python run_TrAnnoScope.py all -c config/test_config.yaml

# If you want to run TrAnnoScope slurm cluster 
sbatch slurm_submit.sh all [-A <slurm account name>]
```


## License
TrAnnoScope is licensed under the [MIT License](LICENSE).

## Contributions
Contributions are welcome. Please fork the repository and submit a pull request.

## Contact
For any questions or issues, please contact [aysevilpektas@mbg.au.dk](mailto:aysevilpektas@mbg.au.dk).
