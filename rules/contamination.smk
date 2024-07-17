import glob
import os 

# Define config file:
configfile: "config/config.yaml"

# Output directory path 
OUTDIR = "output_" + config["PROJECT"] + "/pacbio"

##------------------------- Define your input file ---------------------------##

# if you already have preprocessed (High Quality) pacbio reads
if config["preprocessed_long_reads"] != "":
    if os.path.isdir(config["preprocessed_long_reads"]):
        print("Your LR input defined in the config file will be used for contamination removal step.")
        print("\nOnly your desired input files in the directory should contain .fasta extension")
        ISOSEQ_DIR = config['preprocessed_long_reads']                                          # path
        fasta_files = glob.glob(os.path.join(config["preprocessed_long_reads"] , "*.fasta"))    # list of input files
        #LR_SAMPLES = get_samples(config['preprocessed_long_reads'], "*.fasta")                  # wildcards
        LR_SAMPLES = [os.path.basename(path).split('.')[0] for path in fasta_files]
        LR_SAMPLES = set(sorted(LR_SAMPLES))
    else:
        print(f"Error: The path for LR provided in the config file does not exist: {config['preprocessed_long_reads']}")
        LR_SAMPLES = None  # Or handle the error accordingly
elif config["preprocessed_long_reads"] == "" and os.path.isdir(OUTDIR + "/clustered"):
    print("LR input files created by ISOSEQ preprocessing step will be used for contamination removal step.")
    ISOSEQ_DIR = OUTDIR + "/clustered"                                                                   # path
    LR_SAMPLES, = glob_wildcards(OUTDIR + "/clustered/{sample}.clustered.fasta")                         # wildcards
    LR_SAMPLES = set(sorted(LR_SAMPLES))
    fasta_files = [os.path.join(OUTDIR, f"clustered/{sample}.clustered.fasta") for sample in LR_SAMPLES] # list of input files

else:
    print("You don't have a input file. Define the path of your input directory correctly or run the previous steps of the workflow.")
    LR_SAMPLES = None 

# rest of the suffix of the input file.
LR_file_suffix = '.'.join(os.path.basename(fasta_files[0]).split('.')[1:]) # Get file name except wildcard part


# if you already have clean short reads
if config["clean_short_reads"] != "":
    if os.path.isdir(config["clean_short_reads"]):
        print("Your SR input defined in the config file will be used for contamination removal step.") 
        print("\nOnly your desired input files in the directory should contain .fq.gz extension")      
        RNASEQ_DIR = config['clean_short_reads']                                        # path
        fq_files = glob.glob(os.path.join(config["clean_short_reads"] , "*.fq.gz"))     # list of files
        SR_SAMPLES = [os.path.basename(path).split('_')[0] for path in fq_files]        # wildcards SAMPLE 
        SR_SAMPLES = set(sorted(SR_SAMPLES))
        SR_file_suffix = fq_files[0].split('.', 1)[1]                                   # file suffix, everything except <sample>
    else:
        print(f"Error: The path provided in the config file does not exist: {config['clean_short_reads']}")
        SR_SAMPLES = None 
elif config["clean_short_reads"] == "" and os.path.isdir(OUTDIR + "/../illumina/trimmed_fastp"):
    print("SR Input files created by RNASEQ preprocessing step will be used for contamination removal step.")
    RNASEQ_DIR = OUTDIR + "/../illumina/trimmed_fastp"                                                      # path
    fq_files = glob.glob(os.path.join(OUTDIR + "/../illumina/trimmed_fastp/", "*.fq.gz"))                   # list of files
    SR_SAMPLES, FRR = glob_wildcards(OUTDIR + "/../illumina/trimmed_fastp/{sample}_{frr}.trimmed.fq.gz")    # wildcards SAMPLE and FR
    SR_SAMPLES = set(sorted(SR_SAMPLES))
    SR_file_suffix = "trimmed.fq.gz"                                                                        # file suffix, everything except SAMPLE and FR
else:
    print("You don't have SR input files. Define the path of your input directory correctly or run the previous steps of the workflow.")
    SR_SAMPLES = None 

# Handling the else case: Ensure that Snakemake will stop execution if the transcriptome is not defined
if SR_SAMPLES is None or LR_SAMPLES is None:
    raise ValueError("Input files are not correctly defined. Please provide the path in the config file or ensure the previous steps are completed.")
elif SR_SAMPLES != LR_SAMPLES:
    print(f"Input LR directory: {ISOSEQ_DIR} and SR directory: {RNASEQ_DIR}")
    raise ValueError("Wildcards for short are long reads do not overlaps, Please ensure the files have same sample wildcards!")
else:
    print(f"Input LR directory: {ISOSEQ_DIR} and SR directory: {RNASEQ_DIR}")
    print(f"Found LR samples: {LR_SAMPLES} and SR samples {SR_SAMPLES}")
    print(f"Input LR files: {fasta_files} and SR files {fq_files}")
    print(f"LR file suffix: {LR_file_suffix} and SR file suffix: {SR_file_suffix}")


##------------------------- Find the number of chunks produced by SEQKIT ---------------------------##

# Function to count sequences in a FASTA file
def count_sequences(fasta_file):
    with open(fasta_file, 'r') as file:
        return sum(1 for line in file if line.startswith('>'))

# Find the FASTA file with the highest number of sequences
max_sequences = max(count_sequences(f) for f in fasta_files)

# Calculate the number of chunks needed
chunk_size = config["blastn"]["split_size"]
chunks = -(-max_sequences // chunk_size)  # Equivalent to math.ceil(max_sequences / chunk_size)
#print(chunks)
# Print the results (for debugging purposes)
#print(f"Number of chunks needed: {chunks}")

# Indices for the fasta split files, in the form of '001', '002'...
indexes = ['{0:03d}'.format(x) for x in range(1, chunks+1)]
#print(indexes)

##-------------------------------------------------------------------------------------------------##

rule all:
    "Target files for the workflow."
    input:
        expand(OUTDIR + "/clustered/newheader/{sample}.clustered.newheader.fasta", sample = LR_SAMPLES),
        expand(OUTDIR + "/remove_contaminants/busco_clustered/" + config["busco"]["lineage"] + "_{sample}/run_" + config["busco"]["lineage"] + "/short_summary.txt", sample = LR_SAMPLES),
        expand(OUTDIR + "/remove_contaminants/coverage/{sample}.clustered.newheader.bam", sample = LR_SAMPLES),
        expand(OUTDIR + "/remove_contaminants/blastn/chunks/{sample}.clustered.newheader.part_{index}.fasta", sample = LR_SAMPLES, index=indexes),
        #expand(OUTDIR + "/remove_contaminants/blastn/chunks/{sample}.part_{chunk}.outfmt6", sample = SAMPLES, chunk = chunks),
        expand(OUTDIR + "/remove_contaminants/blastn/{sample}.merged.outfmt6", sample=LR_SAMPLES),
        #expand(OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}.done", sample = SAMPLES),
        expand(OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}.done", sample = LR_SAMPLES), 
        expand(OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/figures/blobdir_{sample}.{plot}."+config["blobtools"]["plot_extension"], sample = LR_SAMPLES, plot = ["snail", "cumulative", "blob.circle"]),
        expand(OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/table_{sample}.tsv", sample = LR_SAMPLES),
        expand(OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/{sample}.taxa_filtered.fasta", sample = LR_SAMPLES),
        #expand(OUTDIR + "/remove_contaminants/clean_reads/{sample}.done", sample = LR_SAMPLES)
        expand(OUTDIR + "/remove_contaminants/clean_reads/{sample}.clean.fasta", sample = LR_SAMPLES)


ruleorder: modify_header > split_fasta > run_blastn > blastn_merge_outfmt6

rule modify_header: 
    "BUSCO does not accept '/' in the header. Change them with '_' and add sample names as suffix."
    input:
        fasta_files
    output:
        #touch(OUTDIR + "/clustered/newheader/header_modified.done")
        expand(OUTDIR + "/clustered/newheader/{sample}.clustered.newheader.fasta", sample=LR_SAMPLES)
    log:
        "logs/pacbio/modify_header.log"
    params:
        dir = OUTDIR + "/clustered/newheader"
    threads:
        1
    resources:
        mem_mb = 100
    shell:
        """ 
        for file in {input}; do
            sample=$(basename $file .clustered.fasta)
            bash scripts/modify_header.sh $file {params.dir}/$sample.clustered.newheader.fasta $sample {log}
        done
        """

rule busco:  
    """
    Assessing the completeness of a transcriptome or genome assembly by comparing it to a set of BUSCO genes.
    """
    input:
        #touch(OUTDIR + "/clustered/newheader/header_modified.done"),
        OUTDIR + "/clustered/newheader/{sample}.clustered.newheader.fasta"
    output:
        OUTDIR + "/remove_contaminants/busco_clustered/" + config["busco"]["lineage"] + "_{sample}/run_" + config["busco"]["lineage"] + "/short_summary.txt"
    log:
        "logs/pacbio/remove_contaminants/busco_clustered_{sample}.log"
    conda:
        "../envs/busco.yaml"
    params:
        memory = config["busco"]["memory"],
        lineage = config["busco"]["lineage"],
        outdir = OUTDIR + "/remove_contaminants/busco_clustered",
        data_path = config["busco"]["data_path"]
    threads:
        config["busco"]["threads"]
    resources:
        mem_mb = config["busco"]["memory"]
    benchmark:
        "benchmarks/pacbio/remove_contaminants/busco_clustered_{sample}.benchmark.txt"
        #busco -f -i {input} -o {params.lineage}_{wildcards.sample} --out_path {params.outdir} --auto_lineage -m transcriptome -c {threads} --download_path {params.data_path} 2> {log}
        #busco -f -i {input} -o {params.lineage}_{wildcards.sample} --out_path {params.outdir} -l {params.lineage} -m transcriptome -c {threads} --download_path {params.data_path} 2> {log}
    shell:
        """
        if [ -z {params.lineage} ] 2> {log}
        then
            busco -f -i {input} -o {params.lineage}_{wildcards.sample} --out_path {params.outdir} --auto_lineage -m transcriptome -c {threads} 2> {log}
        else
            busco -f -i {input} -o {params.lineage}_{wildcards.sample} --out_path {params.outdir} -l {params.lineage} -m transcriptome -c {threads} --download_path {params.data_path} 2> {log}
        fi
        """
        
rule bowtie2_index:
    """
    Index the reference FASTA file (FL reads from each sample), before mapping short reads to them.
    """
    input:
        OUTDIR + "/clustered/newheader/{sample}.clustered.newheader.fasta"
    output:
        touch(OUTDIR + "/remove_contaminants/coverage/{sample}_index.done")
    log:
        "logs/pacbio/remove_contaminants/bowtie2_index.{sample}.log"
    conda:
        "../envs/bowtie2.yaml"
    params:
        out = OUTDIR + "/remove_contaminants/coverage/{sample}.clustered.newheader.fasta"
    threads:
        config["bowtie2"]["indexing_threads"]
    resources:
        mem_mb = config["bowtie2"]["indexing_memory"]
    benchmark:
        "benchmarks/pacbio/remove_contaminants/coverage/bowtie2_indexing.{sample}.benchmark.txt"
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.out} 2> {log}
        """

if config["END"] == "PE":
    rule bowtie2_mapping_PE:
        """
        Map short PE reads to the reference FASTA file (FL reads from each sample).
        """
        input:
            f = OUTDIR + "/remove_contaminants/coverage/{sample}_index.done",
            R1 = RNASEQ_DIR + "/{sample}_1." + SR_file_suffix,
            R2 = RNASEQ_DIR + "/{sample}_2." + SR_file_suffix,
        output:
            bam = OUTDIR + "/remove_contaminants/coverage/{sample}.clustered.newheader.bam",
            report = OUTDIR + "/remove_contaminants/coverage/{sample}.align_stats.txt"
        log:
            "logs/pacbio/remove_contaminants/coverage/bowtie2_mapping_PE.{sample}.log"
        conda:
            "../envs/bowtie2.yaml"
        params:
            parameters = config["bowtie2"]["mapping_parameters"],
            index = OUTDIR + "/remove_contaminants/coverage/{sample}.clustered.newheader.fasta"
        threads:
            config["bowtie2"]["mapping_threads"]   
        resources:
            mem_mb = config["bowtie2"]["mapping_memory"]
        benchmark:
            "benchmarks/pacbio/remove_contaminants/coverage/bowtie2_mapping_PE.{sample}.benchmark.txt"
        shell:
            """
            bowtie2 -p {threads} -q {params.parameters} -x {params.index} \
                -1 {input.R1} -2 {input.R2} 2> {output.report} | \
                samtools view -Sb -@ 1 -m 1G | samtools sort -@ 1 -m 1G -o {output.bam} 2>{log}
            """

elif config["END"] == "SE":
    rule bowtie2_mapping_SE:
        """
        Map short SE reads to the reference FASTA file (FL reads from each sample).
        """
        input:
            f = OUTDIR + "/remove_contaminants/coverage/{sample}_index.done",
            R1 = RNASEQ_DIR + "/{sample}_1." + SR_file_suffix,
        output:
            bam = OUTDIR + "/remove_contaminants/coverage/{sample}.clustered.newheader.bam",
            report = OUTDIR + "/remove_contaminants/coverage/{sample}.align_stats.txt"
        log:
            "logs/pacbio/remove_contaminants/coverage/bowtie2_mapping_SE.{sample}.log"
        conda:
            "../envs/bowtie2.yaml"
        params:
            parameters = config["bowtie2"]["mapping_parameters"],
            index = OUTDIR + "/remove_contaminants/coverage/{sample}.clustered.newheader.fasta"
        threads:
            config["bowtie2"]["mapping_threads"]   
        resources:
            mem_mb = config["bowtie2"]["mapping_memory"]
        benchmark:
            "benchmarks/pacbio/remove_contaminants/coverage/bowtie2_mapping_SE.{sample}.benchmark.txt"
        shell:
            """
            bowtie2 -p {threads} -q {params.parameters} -x {params.index} \
                -U {input.R1} 2> {output.report} | \
                samtools view -Sb -@ 1 -m 1G | samtools sort -@ 1 -m 1G -o {output.bam} 2>{log}
            """

rule split_fasta:
    """
    Split fasta files into chunks for parallelization.
    """
    input: 
        OUTDIR + "/clustered/newheader/{sample}.clustered.newheader.fasta"
    output: 
        expand(OUTDIR + "/remove_contaminants/blastn/chunks/{{sample}}.clustered.newheader.part_{index}.fasta", index=indexes)
    log: 
        "logs/pacbio/remove_contaminants/split_fasta.{sample}.log"
    params:
        size = chunks,
        dir = OUTDIR + "/remove_contaminants/blastn/chunks"
    conda: 
        "../envs/bowtie2.yaml"
    threads:    
        1
    resources: 
        mem_mb=1000
    shell:
        """
        mkdir -p {params.dir}
        seqkit split --by-part {params.size} {input} -O {params.dir} 2> {log}
        """

rule run_blastn:
    input:
        fasta_chunk = OUTDIR + "/remove_contaminants/blastn/chunks/{sample}.clustered.newheader.part_{index}.fasta"
    output:
        blast_output = OUTDIR + "/remove_contaminants/blastn/chunks/{sample}.part_{index}.outfmt6"
    log:
        "logs/pacbio/remove_contaminants/blastn_{sample}.part_{index}.log"
    params:
        db = config["blastn"]["database"],
        parameters = config["blastn"]["params"],
        outfmt = config["blastn"]["outfmt6"],
        outdir = OUTDIR + "/remove_contaminants/blastn/chunks"
    conda:
        "../envs/blast.yaml"
    threads:
        config["blastn"]["threads"]
    resources:
        mem_mb = config["blastn"]["memory"]
    benchmark:
        "benchmarks/pacbio/remove_contaminants/blastn/parallel_blastn_{sample}.part_{index}.benchmark.txt"
    shell:
        """
        blastn -query {input.fasta_chunk} -db {params.db} -outfmt '{params.outfmt}' {params.parameters} -num_threads {threads} -out {output.blast_output} 2> {log}
        """

rule blastn_merge_outfmt6:
    """
    Merge all blastn results per sample into one file.
    """
    input:
        expand(OUTDIR + "/remove_contaminants/blastn/chunks/{sample}.part_{index}.outfmt6", sample = LR_SAMPLES, index=indexes)
    output:
        OUTDIR + "/remove_contaminants/blastn/{sample}.merged.outfmt6"
    log: 
        "logs/pacbio/blobtools/blastn_merge_outfmt6.{sample}.log"
    threads:    
        1
    resources:  
        mem_mb=1000
    shell:  
        """
        cat {input} | grep {wildcards.sample} > {output} 2>{log}
        """

rule create_Blobdir:
    """
    Create Blobtools directory for each sample.
    """
    input: 
        fasta = OUTDIR + "/clustered/newheader/{sample}.clustered.newheader.fasta", 
        coverage = OUTDIR + "/remove_contaminants/coverage/{sample}.clustered.newheader.bam", 
        blastn = OUTDIR + "/remove_contaminants/blastn/{sample}.merged.outfmt6",
        busco = OUTDIR + "/remove_contaminants/busco_clustered/" + config["busco"]["lineage"] + "_{sample}/run_" + config["busco"]["lineage"] + "/short_summary.txt"
    output: 
        #dir = directory(OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}"),
        check = touch(OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}.done")
    log: 
        "logs/pacbio/blobtools/create_Blobdir.{sample}.log"
    params: 
        taxdir = config["blobtools"]["taxdir"],
        busco_tsv = OUTDIR + "/remove_contaminants/busco_clustered/" + config["busco"]["lineage"] + "_{sample}/run_" + config["busco"]["lineage"] + "/full_table.tsv",
        dir = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}" 
    conda: 
        "../envs/blobtools2.yaml"
    threads:   
        config["blobtools"]["threads"] #2
    resources: 
        mem_mb=config["blobtools"]["memory"] #4000
    benchmark: 
        "benchmarks/pacbio/remove_contaminants/blobtools/create_blobdir_{sample}.txt"
    shell:
        """
        blobtools add \
            --fasta {input.fasta} \
            --cov {input.coverage} \
            --hits {input.blastn} \
            --busco {params.busco_tsv} \
            --taxdump {params.taxdir} \
            --create {params.dir} --threads {threads} 2>{log}
        """

rule create_Blobtool_results:
    """
    Create plots from Blobdir directories.
    """
    input: 
        OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}.done"
        #OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}"
    output: 
        expand(OUTDIR + "/remove_contaminants/blobtools/blobdir_{{sample}}/figures/blobdir_{{sample}}.{plot}."+config["blobtools"]["plot_extension"], plot = ["snail", "cumulative", "blob.circle"]),
        #snail = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/figures/blobdir_{sample}.snail."+config["blobtools"]["plot_extension"],
        #blob = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/figures/blobdir_{sample}.blob.circle."+config["blobtools"]["plot_extension"],
        #cumulative = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/figures/blobdir_{sample}.cumulative."+config["blobtools"]["plot_extension"]
    conda: 
        "../envs/blobtools2.yaml"
    log: 
        "logs/pacbio/remove_contaminants/blobtools_figures_Blobdir.{sample}.log"
    params:
        outdir = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/figures",
        ext = config["blobtools"]["plot_extension"],
        dir = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}"
    threads:   
        config["blobtools"]["threads"] #2
    resources: 
        mem_mb=config["blobtools"]["memory"] #4000
    shell:
        """
        mkdir -p {params.dir}
        blobtools view --plot --format {params.ext} --out {params.outdir} --view snail {params.dir} 2>>{log} 
        blobtools view --plot --format {params.ext} --out {params.outdir} --view blob {params.dir} 2>>{log} 
        blobtools view --plot --format {params.ext} --out {params.outdir} --view cumulative {params.dir} 2>>{log}
        """

rule blobtools_table:
    """
    Create table from Blobdir directory.
    """
    input: 
        #OUTDIR + "clustered/{sample}.clustered.newheader.filtered.fasta"
        #checkpoint = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}.done",
        OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}.done"
    output: # checkpoints _must_ have output.
        #touch(".filtered_fasta.touch")
        OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/table_{sample}.tsv"
    params:
        dir = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}"
    log: 
        "logs/pacbio/remove_contaminants/blobtools_table.{sample}.log"
    conda: 
        "../envs/blobtools2.yaml"
    threads: 
        2
    shell:
        """
        blobtools filter \
            --table {output} \
            --table-fields gc,length,{wildcards.sample}.clustered.newheader_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
            {params.dir} &>> {log}
        """

rule filter_taxa:
    """
    Filter contigs based on taxa and obtain filtered fasta.
    """
    input:
        fasta = OUTDIR + "/clustered/newheader/{sample}.clustered.newheader.fasta",
        blobtools_table = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/table_{sample}.tsv",
    output:
        filtered_ids = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/{sample}_filtered_ids.txt",
        filtered_fasta = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/{sample}.taxa_filtered.fasta"
    log: 
        "logs/pacbio/remove_contaminants/blobtools_filtered_taxa.{sample}.log"
    conda: 
        "../envs/blobtools2.yaml"
    threads: 
        1
    resources:
        mem_mb=1000
    shell:
        """
        python scripts/filter_taxa.py {input.blobtools_table} {output.filtered_ids} &>> {log}
        seqkit grep -f {output.filtered_ids} {input.fasta} > {output.filtered_fasta}
        """

rule remove_extra:
    input:
        filtered_fasta = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}/{sample}.taxa_filtered.fasta",
        blastn = OUTDIR + "/remove_contaminants/blastn/{sample}.merged.outfmt6"
    output: 
        #touch(OUTDIR + "/remove_contaminants/clean_reads/{sample}.done")
        OUTDIR + "/remove_contaminants/clean_reads/{sample}.clean.fasta"
    log:
        "logs/pacbio/remove_contaminants/blobtools/removal_possible_MT_rRNA.{sample}.log"
    params:
        organism = config["organism"],
        qseqid_dir = OUTDIR + "/remove_contaminants/blobtools/blobdir_{sample}",
        mt_removal = config["MT_removal"],
        rRNA_removal = config["rRNA_removal"],
        dir = OUTDIR + "/remove_contaminants/clean_reads"
    threads:
        1
    resources:
        mem_mb=1000
    conda:
        "../envs/blobtools2.yaml"
    shell:
        """
        mkdir -p {params.dir}
        
        if [[ "{params.mt_removal}" == "yes" && "{params.rRNA_removal}" == "yes" ]]; then
            bash scripts/MT_removal.sh {input.blastn} {{params.organism}} {input.filtered_fasta} {params.qseqid_dir}/{wildcards.sample}.mt_nt.qseqid &>> {log}
            bash scripts/rRNA_removal.sh {input.blastn} {input.filtered_fasta} {params.qseqid_dir}/{wildcards.sample}.rRNA_nt.qseqid &>> {log}
            #cat {params.qseqid_dir}/{wildcards.sample}.mt_nt.qseqid {params.qseqid_dir}/{wildcards.sample}.rRNA_nt.qseqid | grep -E "^$" -v | sort | uniq > {params.qseqid_dir}/{wildcards.sample}.mt_rRNA_nt.qseqid
            echo "Running seqkit" &>> {log}
            seqkit grep -n -v -f {params.qseqid_dir}/{wildcards.sample}.mt_nt.qseqid {input.filtered_fasta} | seqkit grep -n -v -f {params.qseqid_dir}/{wildcards.sample}.rRNA_nt.qseqid -o {output} &>> {log}
        elif [[ "{params.mt_removal}" == "yes" && "{params.rRNA_removal}" == "no" ]]; then
            bash scripts/MT_removal.sh {input.blastn} {{params.organism}} {input.filtered_fasta} {params.qseqid_dir}/{wildcards.sample}.mt_nt.qseqid &>> {log}
            seqkit grep -n -v -f {params.qseqid_dir}/{wildcards.sample}.mt_nt.qseqid {input.filtered_fasta} -o {output} 2>> {log}
        elif [[ "{params.mt_removal}" == "no" && "{params.rRNA_removal}" == "yes" ]]; then
            bash scripts/rRNA_removal.sh {input.blastn} {input.filtered_fasta} {params.qseqid_dir}/{wildcards.sample}.rRNA_nt.qseqid &>> {log}
            seqkit grep -n -v -f {params.qseqid_dir}/{wildcards.sample}.rRNA_nt.qseqid {input.filtered_fasta} -o {output} 2>> {log}
        else
            echo "Skipping MT and rRNA removal" &>> {log}
            mkdir -p {params.dir}
            cp {input.filtered_fasta} {output}
        fi
        """

    
          
