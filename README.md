# TrAnnoScope

**TrAnnoScope** is a comprehensive workflow designed to process PacBio and Illumina reads to generate full-length transcripts and annotate them effectively. 

TrAnnoScope is based on the Snakemake workflow management system. This tool offers a comprehensive suite of features designed to generate the full-length transcriptome from multiple tissues using PacBio and Illumina reads, as well as for functional annotation analysis of the transcriptome. The workflow consists of eight sub-workflows and each part can also be run independently for specific analyses. The initialization part involves the installation of the necessary bioinformatics packages followed by the download and preparation of the relevant reference databases. In the preprocessing step, several quality control tools are used to evaluate and improve the quality of raw Illumina reads. For PacBio reads, isoseq3 is used to obtain high-quality reads.
Additionally, optional error correction and contamination-removal steps are defined. FMLRC can be used to increase the quality of long reads obtained by clean Illumina reads, while Blobtools remove contamination from long reads using coverage information obtained from clean Illumina reads and BUSCO results. During the assembly and reduction step, CDHIT-Est eliminates redundancy between tissue samples, and EvidentialGene removes duplicate reads and obtains protein-coding reads. A quality evaluation of an assembly is conducted based on the completeness, the length distribution, and the ratio of nearly full-length transcripts against databases of well-known biological sequences such as UniProt or closely related species. Trinotate software is used to perform functional annotation of transcripts against various databases as part of the annotation process. Trinotate offers sequence databases, including Pfam, SwissProt, SignalP, TMHMM, EggNog, and Infernal. Additionally, homology searches against the NR and NT databases were added to the workflow as an option. But the users should prepare the necessary files for themselves. As a whole, the workflow includes quality control, reduction of redundancy, obtaining protein-coding transcripts, evaluating assembly, and annotating full-length transcripts and protein sequences.

## Workflow
<img width="409" alt="image" src="https://github.com/user-attachments/assets/555fa330-bced-4589-9f56-b7b8dd3d60ee">


## Features
- **Quality Control**: Performs QC on Illumina short reads.
- **Preprocessing**: Filters and trims Illumina reads, processes PacBio long reads
- **Contaminants Removal** Removes contamination from PacBio long reads (Blobtools2).
- **Error Correction**: Corrects errors in PacBio reads by Illumina reads (FMLRC).
- **Classification**: Cluster and Classify PacBio reads.
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
  3. **Install necessary packages and databases**
  
  The Python script  called precheck.py is designed to manage the prerequisites and steps involved in running the TrAnnoScope pipeline. It utilizes Snakemake for workflow management and integrates various bioinformatics tools and databases. The script ensures that necessary environments and resources are installed and configured, guiding users through each step of the pipeline setup and execution process. 

  **NOTE:** If you want to use SignalP and TmHMM2, first you should install tar.gz file of them and store it in resources/
  
  ```bash
  positional arguments:
    STEPS                 Which step do you want to run
      all                 Run full analysis in FLTransAnnot pipeline
      qc_rnaseq           Quality Control for Illumina short reads
      preprocessing_rnaseq
                          Filtering and Trimming of Illumina short reads
      preprocessing_pacbio
                          Processing PacBio long reads
      remove_contaminants
                          Contamination removal of PacBio long reads
      error_correction    Error correction of PacBio long reads
      classification      Clustering of PacBio long reads
      annotation          Annotation of PacBio long reads
      quality_assessment  Quality assessment of transcriptome

  Usage: python precheck.py [STEPS] -c config/config.yaml
  ```
**For preprocessing_rnaseq step:** 

it will ask to download FastQScreen Genome indexes. 

```bash
python precheck.py preprocessing_rnaseq -c config/config.yaml
```

Additionally, you can download SILVA LSU and SSU and also MT sequences of your organism of interest.  
You can download SILVA DB from "https://www.arb-silva.de/no_cache/download/archive/current/Exports/SILVA_138.1_LSURef_tax_silva.fasta.gz" and 	SILVA_138.1_SSURef_tax_silva.fasta.gz (these are the current one for now)
For MT DB, Check your organism in NCBI for MT genome

Then create bowtie2 index for them with bowtie2-build --threads <THREADS> <INPUT> <OUTPUT>
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

it will ask to download the interested lineage file for BUSCO, if you already performed, contamination removal, you can say "NO".

```bash
python precheck.py quality_assessment -c config/config.yaml
```

**If you run ALL step:** 

it will ask about all necessary files to download

```bash
python precheck.py all -c config/config.yaml
```

## Usage
**Running on local computer:**

Configure settings in config/config.yaml.

    ```bash
    STEP                  Which step you want to run
      qc_rnaseq           Quality Control for Illumina short reads
      preprocessing_rnaseq
                          Filtering and Trimming of Illumina short reads
      preprocessing_pacbio
                          Processing PacBio long reads
      remove_contaminants
                          Contamination removal of PacBio long reads
      error_correction    Error correction of PacBio long reads
      classification      Clustering of PacBio long reads
      annotation          Annotation of PacBio long reads
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
                 qc_rnaseq             - Quality Control for Illumina short reads
                 preprocessing_rnaseq  - Filtering and Trimming of Illumina short reads
                 preprocessing_pacbio  - Processing PacBio long reads
                 remove_contaminants   - Contamination removal of PacBio long reads
                 error_correction      - Error correction of PacBio long reads
                 classification        - Clustering of PacBio long reads
                 annotation            - Annotation of PacBio long reads
                 quality_assessment    - Quality assessment of the transcriptome
                 all                   - Run all the steps

    Usage: bash slurm_submit.sh STEP [-A <slurm account name>]
    Example: bash slurm_submit.sh qc_rnaseq [-A <slurm account_name>]
    ```


### Test Data

To verify the workflow, you can use the provided test data located in the `data/test_data/` directory.

1. **Run the workflow with the test data**:
   
```bash
# activate conda environment
conda activate trannoscope

# Prepare the necessary files
python precheck.py all -c config/test_config.yaml

# If you want to run the TrAnnoScope on local computer
python run_TrAnnoScope.py all -c config/test_config.yaml

# If you want to run TrAnnoScope slurm cluster
bash slurm_submit.py all [-A <slurm account name>]
```


## License
TrAnnoScope is licensed under the [MIT License](LICENSE).

## Contributions
Contributions are welcome. Please fork the repository and submit a pull request.

## Contact
For any questions or issues, please contact [aysevilpektas@mbg.au.dk](mailto:aysevilpektas@mbg.au.dk).
