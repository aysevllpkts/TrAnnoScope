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
- Programs used by TrAnnoScope
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
- Databases used by TrAnnoScope
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
  conda env creae -n trannoscope -f TrAnnoScope.yaml
  conda activate trannoscope
  ```
  3. **Install necessary packages and databases**
  
  Precheck.py is a helper python script to load the necessary conda environments followed by the download and preparation of the relevant reference databases like FastQscreen genome DBs, BUSCO lineage DB, Taxdump files for Blobtools, Trinotate DBs and registering SignalP and TmHMM2.
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
For preprocessing_rnaseq, it will ask to download FastQScreen Genome indexes

```bash
python precheck.py preprocessing_rnaseq -c config/config.yaml
```

For remove_contaminants, it will ask to download Bnterested lineage files for BUSCO and taxdump which is used by Blobtools2.

```bash
python precheck.py remove_contaminants -c config/config.yaml
```

For annotation step, it will ask to download TRINOTATE_DATA_DIR for databases and register SignalP and TmHMM2

```bash
python precheck.py annotation -c config/config.yaml
```

For quality_assessment step, it will ask to download the interested lineage file for BUSCO, if you already performed, contamination removal, you can say "NO".

```bash
python precheck.py quality_assessment -c config/config.yaml
```

## Usage
**Running on local computer:**

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
    Example: bash slurm_submit.sh qc_rnaseq [-A <account_name>]
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
