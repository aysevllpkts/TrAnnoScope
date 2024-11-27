# Aysevil Pektas
# 25-03-2024
# Trinotate annotation pipeline

########### NOTES ##########
# if you have a problem to create conda environment for trinotate=4.0.2
# USE mamba clean --all, then try again.
# add DBs paths and installation rules if necesssary
# add trinotate_init to register file path for signalp6 and tmhmmv2
# --run "swissprot_blastp swissprot_blastx pfam signalp6 tmhmmv2 infernal EggnogMapper" \
############################

import os
import glob

configfile: "config/config.yaml"

OUTDIR = "output_" + config["PROJECT"] 

##------------------------- Define your transcriptome ---------------------------##

if config["input_for_annotation"]["nucl"] != "" and config["input_for_annotation"]["prot"] != "":
    if os.path.isfile(config["input_for_annotation"]["nucl"]) and config["input_for_annotation"]["prot"]: 
        print("Nucleotide and protein sequences of your transcriptome from the config file will be used for annotation step.")
        transcriptome_nucl = config["input_for_annotation"]["nucl"]
        transcriptome_prot = config["input_for_annotation"]["prot"]
    else:
        raise ValueError(f"Error: The path provided in the config file does not exist: {config['input_for_annotation']}")
        #transcriptome_= None  # Or handle the error accordingly
    #print("Your transcriptome from config file will be used for annotation.")
    #transcriptome = config["fasta_for_annotation"]
elif config["input_for_annotation"]["nucl"] == "" and os.path.isfile(OUTDIR + "/pacbio/classification/evigene/okayset/merged_clean_corrected.okay.mrna"):
    print("Transcriptome file created by classification step will be used for annotation.")
    transcriptome_nucl = os.path.join(OUTDIR + "/pacbio/classification/evigene/okayset/merged_clean_corrected.okay.mrna")
    transcriptome_prot = os.path.join(OUTDIR + "/pacbio/classification/evigene/okayset/merged_clean_corrected.okay.aa")
else:
    raise ValueError("Transcriptome is not defined. Please provide the path in the config file or ensure the previous steps are completed.")
    transcriptome_nucl = None  # Or set it to an error value
    transcriptome_nucl = None

##------------------------- Find the number of chunks produced by SEQKIT ---------------------------##

# Function to count sequences in a FASTA file
def count_sequences(file):
    with open(transcriptome_nucl, 'r') as file:
        return sum(1 for line in file if line.startswith('>'))

# Find the FASTA file with the highest number of sequences
max_sequences = count_sequences(transcriptome_nucl)

# Calculate the number of chunks needed
chunk_size = config["trinotate"]["split_size"]
chunks = -(-max_sequences // chunk_size)  # Equivalent to math.ceil(max_sequences / chunk_size)
#print(chunks)
# Print the results (for debugging purposes)
#print(f"Number of chunks needed: {chunks}")

# Indices for the fasta split files, in the form of '001', '002'...
indexes = ['{0:03d}'.format(x) for x in range(1, chunks+1)]
#print(indexes)

##-------------------------------------------------------------------------------------------------##

# For optional database runs:
conditional_outputs = list()
nr_plot = list()
eggnog_plots = list()

if "swissprot_blastx" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/swissprot_blastx/swissprot_blastx_trinotate_load.done")
if "swissprot_blastp" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/swissprot_blastp/swissprot_blastp_trinotate_load.done")
if "pfam" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/pfam/pfam_trinotate_load.done")
if "signalp6" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/signalp6/signalp6_trinotate_load.done") 
if "tmhmmv2" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/tmhmm2/tmhmm2_trinotate_load.done")
if "eggnog_mapper" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/eggnog_mapper/eggnog_mapper_load.done")
    eggnog_plots.append(OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_KEGG_dist." + config["format"])
    eggnog_plots.append(OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_KEGG_slim_dist." + config["format"])
    eggnog_plots.append(OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_COG_dist." + config["format"])
if "infernal" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/infernal/infernal_load.done")
if "nr_blastx" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/nr_blastx/nr_blastx_load.done")
if "nr_blastp" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/nr_blastp/nr_blastp_load.done")
if "nt" in config["database_run"]:
    conditional_outputs.append(OUTDIR + "/annotation/nt_blastn/nt_trinotate_load.done")
    nr_plot.append(OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_top10_species_dist_NR_blastx." + config["format"])


#print(conditional_outputs)

localrules: all, copy_transcriptome, modify_prot_header, create_gene_trans_map, copy_trinotate_sqlite, split_pep, split_nucl, get_summary

rule all:
    input:
        OUTDIR + "/transcriptome/transcriptome_nt.fasta",
        #OUTDIR + "/transcriptome/transcriptome_prot.fasta",
        OUTDIR + "/transcriptome/transcriptome_prot_modified.fasta",
        OUTDIR + "/transcriptome/gene_trans_map.txt",
        OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        OUTDIR + "/annotation/trinotate_init.done",
        #expand(OUTDIR + "/annotation/signalp6/chunks/sigP6outdir_{index}", index = indexes),
        #conditional_outputs,
        OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report.xls",
        OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report.summary",
        OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report_trans.goslims",
        expand(OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_{category}_GOSlims.{ext}", category = ["ALL", "Cellular_Component", "Molecular_Function", "Biological_Process"], ext = config["format"]),
        eggnog_plots,
        nr_plot
 
        
checkpoint copy_transcriptome:
    input:
        transcriptome_nt = transcriptome_nucl,
        transcriptome_prot = transcriptome_prot
    output: 
        copied_nt = OUTDIR + "/transcriptome/transcriptome_nt.fasta",
        copied_prot = OUTDIR + "/transcriptome/transcriptome_prot.fasta" 
    log: 
        "logs/transcriptome/transcriptome_copy.log"
    threads: 1
    resources: mem_mb = 1000
    shell:
        """
        cp {input.transcriptome_nt} {output.copied_nt}
        cp {input.transcriptome_prot} {output.copied_prot}
        echo 'transcriptome_nt and transcriptome_prot were copied to transcriptome folder' &> {log}
        """

if config["modify_protein_header"] == "yes":
    rule modify_prot_header:
        """
        Trinotate only except transdecoder like header format. 
        Modify the evigene format to transdecoder-like format.
        """
        input:  
            OUTDIR + "/transcriptome/transcriptome_prot.fasta"
        output:
            OUTDIR + "/transcriptome/transcriptome_prot_modified.fasta"
        log:
            "logs/transcriptome/modify_prot_header.log"
        threads: 1
        resources:  mem_mb=1000 
        shell: "python scripts/reformat_header.py {input} {output} &> {log}" 

else:
    rule copy_prot_file:
        """
        Copy file to match name!
        """
        input:  
            OUTDIR + "/transcriptome/transcriptome_prot.fasta"
        output:
            OUTDIR + "/transcriptome/transcriptome_prot_modified.fasta"
        log:
            "logs/transcriptome/modify_prot_header.log"
        threads: 1
        resources:  mem_mb=1000 
        shell: "cp {input} {output}" 

checkpoint create_gene_trans_map:
    input:
        transcriptome = OUTDIR + "/transcriptome/transcriptome_nt.fasta"
    output:
        gene_trans_map = OUTDIR + "/transcriptome/gene_trans_map.txt"
    log:
        "logs/transcriptome/create_gene_trans_map.log"
    threads: 1
    resources:  mem_mb = 1000
    shell:
        """
        if [ ! -f {output.gene_trans_map} ]; then
            bash scripts/create_gene_trans_map.sh {input.transcriptome} {output.gene_trans_map} &> {log}
            echo "gene_trans_map was created" &>> {log}
        fi
        """

checkpoint copy_trinotate_sqlite:
    input: 
        "resources/myTrinotate.sqlite"
    output:
        copied = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        check = touch(OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite_copied.done")
    log: 
        "logs/annotation/copy_trinotate_sqlite.log"
    threads: 1
    resources: mem_mb = 2000 
    shell:
        """
        cp {input} {output.copied} 2> {log}
        """

rule split_pep:
    """
    Split protein fasta file into chunks for parallelization.
    """
    input: 
        OUTDIR + "/transcriptome/transcriptome_prot_modified.fasta",
    output: 
        expand(OUTDIR + "/transcriptome/split/prot/transcriptome_prot_modified.part_{index}.fasta", index=indexes)
    log: 
        "logs/transcriptome/split_transcriptome_prot.log"
    conda:
        "../envs/bowtie2.yaml"
    params:
        size = chunks,
        outdir = OUTDIR + "/transcriptome/split/prot"
    threads:    
        1
    resources: 
        mem_mb=1000
    shell:
        """
        seqkit split --by-part {params.size} {input} -O {params.outdir} 2> {log}
        """

rule split_nucl:
    """
    Split transcriptome nucl fasta file into chunks for parallelization.
    """
    input: 
        OUTDIR + "/transcriptome/transcriptome_nt.fasta", 
    output: 
        expand(OUTDIR + "/transcriptome/split/nucl/transcriptome_nt.part_{index}.fasta", index=indexes)
    log: 
        "logs/transcriptome/split_transcriptome_nucl.log"
    conda:
        "../envs/bowtie2.yaml"
    params:
        size = chunks,
        outdir = OUTDIR + "/transcriptome/split/nucl"
    threads:    
        1
    resources: 
        mem_mb=1000
    shell:
        """
        seqkit split --by-part {params.size} {input} -O {params.outdir} 2> {log}
        """

rule init_trinotate:
    "Initial import of transcriptome and protein data."
    input:
        transcriptome = OUTDIR + "/transcriptome/transcriptome_nt.fasta",
        protein = OUTDIR + "/transcriptome/transcriptome_prot_modified.fasta",
        gene_trans_map = OUTDIR + "/transcriptome/gene_trans_map.txt",
        check = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite_copied.done"
    output:
        touch(OUTDIR + "/annotation/trinotate_init.done")
    log:
        "logs/annotation/trinotate_init.log"
    conda: 
        "../envs/trinotate_v4.0.2.yaml"
    params:
        db = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        Trinotate = config["trinotate"]["path"]
    threads:
        1
    resources:
        mem_mb = 2000
    shell:
        """
        # Check if the file _init.ok exists in the directory - it is preventing to initiate the sqlite file
        if [ -f "__init.ok" ]; then
            # If the file exists, delete it
            rm "__init.ok"
            echo "__init.ok has been deleted. Ready for new TRINOTATE SQLITE initialization"
        else
             echo "__init.ok does not exist in the current directory."
        fi
        {params.Trinotate}/Trinotate --db {params.db} --init --gene_trans_map {input.gene_trans_map} --transcript_fasta {input.transcriptome} --transdecoder_pep {input.protein}  &> {log}
        echo "Trinotate was initialized" &>> {log}
        """

if "swissprot_blastx" in config["database_run"]:
    rule parallel_swissprot_blastx:
        "Run swissprot blastx in parallel"
        input:
            OUTDIR + "/transcriptome/split/nucl/transcriptome_nt.part_{index}.fasta"
        output:
            OUTDIR + "/annotation/swissprot_blastx/chunks/transcriptome_nt_swissprot_diamond.part_{index}.outfmt6"
        log:
            "logs/annotation/paralell_swissprot_blastx.part_{index}.log"
        params:
            db = config["trinotate"]["trinotate_data_dir"] + "/uniprot_sprot.dmnd",
            outfmt = config["trinotate"]["outfmt"],
            evalue=config["trinotate"]["e-value"]
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads:
            config["parallel_run"]["threads"]
        resources:
            mem_mb = config["parallel_run"]["memory"]
        benchmark:
            "benchmarks/annotation/swissprot_blastx_parallel.part_{index}.benchmark.txt"
        shell:
            """
            diamond blastx --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --max-target-seqs 1 --threads {threads} --evalue {params.evalue} &> {log}
            """

    rule merge_swissprot_blastx:
        """
        Merge all swissprot blastx results.
        """
        input:
            expand(OUTDIR + "/annotation/swissprot_blastx/chunks/transcriptome_nt_swissprot_diamond.part_{index}.outfmt6", index=indexes)
        output:
            OUTDIR + "/annotation/swissprot_blastx/transcriptome_swissprot_blastx.outfmt6"
        log: 
            "logs/annotation/swissprot_blastx_merge.log"
        threads:    
            1
        resources:  
            mem_mb=1000
        shell:  
            """
            cat {input} > {output} 2>{log}
            """

    rule load_swissprot_blastx:
        """
        Load OUTFMT6 to SQLITE DB.
        """
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/swissprot_blastx/transcriptome_swissprot_blastx.outfmt6"
        output:
            touch(OUTDIR + "/annotation/swissprot_blastx/swissprot_blastx_trinotate_load.done")
        log:
            "logs/annotation/swissprot_blastx_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads:
            2
        resources:  mem_mb=2000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_swissprot_blastx {input.hits} 2> {log}
            """

if "swissprot_blastp" in config["database_run"]:
    rule parallel_swissprot_blastp:
        """
        Run swissprot blastp.
        """
        input: 
            OUTDIR + "/transcriptome/split/prot/transcriptome_prot_modified.part_{index}.fasta"
        output: 
            OUTDIR + "/annotation/swissprot_blastp/chunks/transcriptome_prot_swissprot_diamond.part_{index}.outfmt6"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        log: 
            "logs/annotation/swissprot_blastp_parallel.part_{index}.log"
        params:
            db = config["trinotate"]["trinotate_data_dir"] + "/uniprot_sprot.dmnd",
            outfmt = config["trinotate"]["outfmt"],
            evalue=config["trinotate"]["e-value"]
        threads:
            config["parallel_run"]["threads"]
        resources:
            mem_mb = config["parallel_run"]["memory"]
        shell:
            """
            diamond blastp --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --max-target-seqs 1 --threads {threads} --evalue {params.evalue} &> {log}
            """

    rule merge_swissprot_blastp:
        """
        Merge all swissprot blastp results.
        """
        input:
            expand(OUTDIR + "/annotation/swissprot_blastp/chunks/transcriptome_prot_swissprot_diamond.part_{index}.outfmt6", index=indexes)
        output:
            OUTDIR + "/annotation/swissprot_blastp/transcriptome_swissprot_blastp.outfmt6"
        log: 
            "logs/annotation/swissprot_blastp_merge.log"
        threads:    
            1
        resources:  
            mem_mb=1000
        shell:  
            """
            cat {input} > {output} 2>{log}
            """

    rule load_swissprot_blastp:
        """
        Load OUTFMT6 to SQLITE DB.
        """
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/swissprot_blastp/transcriptome_swissprot_blastp.outfmt6"
        output:
            touch(OUTDIR + "/annotation/swissprot_blastp/swissprot_blastp_trinotate_load.done")
        log:
            "logs/annotation/swissprot_blastp_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads:
            2
        resources:  
            mem_mb=2000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_swissprot_blastp {input.hits} 2> {log}
            """

if "pfam" in config["database_run"]:
    rule parallel_pfam:
        """
        Run pfam in parallel.
        """
        input:
            OUTDIR + "/transcriptome/split/prot/transcriptome_prot_modified.part_{index}.fasta",
        output: 
            OUTDIR + "/annotation/pfam/chunks/TrinoattePFAM_{index}.out"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        log: 
            "logs/annotation/pfam_parallel_{index}.log"
        params:
            db = config["trinotate"]["trinotate_data_dir"] + "/Pfam-A.hmm",
        threads:
            config["parallel_run"]["threads"]
        resources:
            mem_mb = config["parallel_run"]["memory"]
        shell:
            """
            hmmsearch --cpu {threads} --noali --domtblout {output} {params.db} {input} > {log}
            """

    rule merge_pfam:
        """
        Merge PFAM outputs. 
        """
        input:
            expand(OUTDIR + "/annotation/pfam/chunks/TrinoattePFAM_{index}.out", index = indexes)
        output:
            OUTDIR + "/annotation/pfam/TrinotatePFAM.out"
        log:
            "logs/annotation/pfam_merge.log"
        threads:
            1
        resources:
            mem_mb = 1000
        shell:
            """
            cat {input} > {output} 2>{log}
            """

    rule load_pfam:
        """
        Load PFAM output to SQLITE db. 
        """
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/pfam/TrinotatePFAM.out"
        output:
            touch(OUTDIR + "/annotation/pfam/pfam_trinotate_load.done")
        log:
            "logs/annotation/pfam_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads: 1
        resources:
            mem_mb = 1000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_pfam {input.hits} &> {log}
            """

if "signalp6" in config["database_run"]:
    rule parallel_signalp6:
        """
        Run signalp6 in parallel
        """
        input:
            OUTDIR + "/transcriptome/split/prot/transcriptome_prot_modified.part_{index}.fasta",
        output:
            OUTDIR + "/annotation/signalp6/chunks/sigP6outdir_{index}/output.gff3"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        log: 
            "logs/annotation/signalp6_parallel_{index}.log"
        params:
            dir = OUTDIR + "/annotation/signalp6/chunks/sigP6outdir_{index}",
        threads: config["parallel_run"]["threads"]
        resources: mem_mb = config["parallel_run"]["memory"]
        shell:
            """
            signalp6 --fastafile {input} --output_dir {params.dir} --format none --organism euk --mode fast &> {log}
            """
    rule merge_signalp6:
        """
        Merge GFF3 file from each signalp6 run 
        """
        input:
            expand(OUTDIR + "/annotation/signalp6/chunks/sigP6outdir_{index}/output.gff3", index = indexes)
        output:
            OUTDIR + "/annotation/signalp6/sigP6out.gff3"
        log:
            "logs/annotation/signalp6_merge.log"
        threads: 1
        resources: mem_mb = 1000
        shell:
            """
            cat {input} > {output} 2>{log}
            """
    rule load_signalp6:
        """
        Load TMHMM2 output to SQLITE db.
        """
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/signalp6/sigP6out.gff3"
        output:
            touch(OUTDIR + "/annotation/signalp6/signalp6_trinotate_load.done")
        log:
            "logs/annotation/signalp6_trinotate_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads:   2
        resources:  mem_mb = 1000 
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_signalp {input.hits}
            """


if "tmhmmv2" in config["database_run"]:
    rule parallel_tmhmm2:
        """
        Run tmhmmv2.
        """
        input: 
            OUTDIR + "/transcriptome/split/prot/transcriptome_prot_modified.part_{index}.fasta"
        output:
            OUTDIR + "/annotation/tmhmm2/chunks/tmhmm2_{index}.out",
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        log: 
            "logs/annotation/tmhmm2_parallel_{index}.log"
        threads: config["parallel_run"]["threads"]
        resources: mem_mb=config["parallel_run"]["memory"]
        shell:
            """
            tmhmm --short {input} > {output}
            """

    rule merge_tmhmm2:
        """
        Merge TMHMM2 outputs file from each signalp6 run 
        """
        input:
            expand(OUTDIR + "/annotation/tmhmm2/chunks/tmhmm2_{index}.out", index = indexes)
        output:
            OUTDIR + "/annotation/tmhmm2/tmhmm2.out"
        log:
            "logs/annotation/tmhmm2_merge.log"
        threads: 1
        resources: mem_mb = 1000
        shell:
            """
            cat {input} > {output} 2>{log}
            """
    rule load_tmhmm2:
        """
        """
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/tmhmm2/tmhmm2.out"
        output:
            touch(OUTDIR + "/annotation/tmhmm2/tmhmm2_trinotate_load.done")
        log:
            "logs/annotation/tmhmm2_trinotate_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads:   2
        resources:  mem_mb = 1000 
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_tmhmmv2 {input.hits}
            """

if "eggnog_mapper" in config["database_run"]:
    rule parallel_eggnog_mapper:
        """
        Run eggnog mapper in parallel for each chunk.
        """
        input: 
            OUTDIR + "/transcriptome/split/prot/transcriptome_prot_modified.part_{index}.fasta"
        output: 
            touch(OUTDIR + "/annotation/eggnog_mapper/chunks/eggnog_mapper_{index}.done")
        log: 
            "logs/annotation/parallel_eggnog_{index}.log"
        params: 
            db = config["trinotate"]["trinotate_data_dir"] + "/EGGNOG_DATA_DIR/",
            dir = OUTDIR + "/annotation/eggnog_mapper/chunks",
            prefix = "transcriptome_prot_modified.part_{index}"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        threads: config["parallel_run"]["threads"]
        resources: mem_mb=config["parallel_run"]["memory"]
        benchmark: 
            "benchmarks/annotation/eggnog_mapper_parallel_{index}.txt"
        shell:
            """
            mkdir -p {params.dir} &> {log}
            emapper.py -i {input} --cpu {threads} --data_dir {params.db} --output {params.prefix} --output_dir {params.dir} &>> {log}
            """

    rule merge_eggnog_mapper:
        """
        Merge all eggnog mapper results per sample into one file.
        """
        input:
            expand(OUTDIR + "/annotation/eggnog_mapper/chunks/eggnog_mapper_{index}.done", index = indexes)
        output:
            annotation = OUTDIR + "/annotation/eggnog_mapper/eggnog_mapper.emapper.annotations",
            hits = OUTDIR + "/annotation/eggnog_mapper/eggnog_mapper.emapper.hits",
            seeds = OUTDIR + "/annotation/eggnog_mapper/eggnog_mapper.emapper.seed_orthologs",
        log: 
            "logs/annotation/eggnog_mapper_merge.log"
        params:
            dir = OUTDIR + "/annotation/eggnog_mapper/chunks"
        threads:    
            1
        resources:  
            mem_mb=1000
        shell: 
            """
            echo "Merge annotation files" 2> {log}
            cat {params.dir}/transcriptome_prot_modified.part_*.emapper.annotations > {output.annotation} 2>> {log}
            echo "Merge hits files" 2>> {log}
            cat {params.dir}/transcriptome_prot_modified.part_*.emapper.hits > {output.hits} 2>> {log}
            echo "Merge seed files" 2>> {log}
            cat {params.dir}/transcriptome_prot_modified.part_*.emapper.seed_orthologs > {output.seeds} 2>> {log}
            """
        
    rule load_eggnog_mapper:
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/eggnog_mapper/eggnog_mapper.emapper.annotations"
        output:
            touch(OUTDIR + "/annotation/eggnog_mapper/eggnog_mapper_load.done")
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        log:
            "logs/annotation/eggnog_mapper_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        threads:
            2
        resources:
            mem_mb=2000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_EggnogMapper {input.hits} &>{log}
            """
        

if "infernal" in config["database_run"]:
    rule parallel_infernal:
        """
        Run infernal in parallel.
        """
        input: 
            OUTDIR + "/transcriptome/split/nucl/transcriptome_nt.part_{index}.fasta"
        output: 
            log = OUTDIR + "/annotation/infernal/chunks/infernal_{index}.log",
            out = OUTDIR + "/annotation/infernal/chunks/infernal_{index}.out"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        params:
            dir = OUTDIR + "/annotation/infernal/chunks",
            db_clanin = config["trinotate"]["trinotate_data_dir"] + "/Rfam.clanin",
            db_cm = config["trinotate"]["trinotate_data_dir"] + "/Rfam.cm",
        threads: config["parallel_run"]["threads"]
        resources: mem_mb=config["parallel_run"]["memory"]
        shell:
            """
            mkdir -p {params.dir}
            cmscan -Z 5 --cut_ga --rfam --nohmmonly --tblout {output.out} --fmt 2 --cpu {threads} --clanin {params.db_clanin} {params.db_cm} {input} > {output.log} 
            """
    
    rule merge_infernal:
        """
        Merge all infernal results per sample into one file.
        """
        input:
            log_files = expand(OUTDIR + "/annotation/infernal/chunks/infernal_{index}.log", index = indexes),
            out_files = expand(OUTDIR + "/annotation/infernal/chunks/infernal_{index}.out", index = indexes)
        output:
            log_file = OUTDIR + "/annotation/infernal/infernal.log",
            out_file = OUTDIR + "/annotation/infernal/infernal.out"
        log: 
            "logs/annotation/infernal_merge.log"
        threads: 1
        resources: mem_mb=1000
        shell: 
            """
            echo "Merge log files" 2> {log}
            cat {input.log_files} > {output.log_file} 2>> {log}
            echo "Merge output files" 2>> {log}
            cat {input.out_files} > {output.out_file} 2>> {log}
            """

    rule load_infernal:
        """
        Load infernal output to SQLITE db. 
        """
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/infernal/infernal.out"
        output:
            touch(OUTDIR + "/annotation/infernal/infernal_load.done")
        log:
            "logs/annotation/infernal_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        threads: 2
        resources:  mem_mb = 1000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_infernal {input.hits} &>{log}
            """

if "nr_blastx" in config["database_run"]:
    rule parallel_nr_blastx:
        """
        Run NR blastx in parallel for each chunk.
        """
        input: 
            OUTDIR + "/transcriptome/split/nucl/transcriptome_nt.part_{index}.fasta"
        output: 
            OUTDIR + "/annotation/nr_blastx/chunks/transcriptome_nt_nr_diamond.part_{index}.outfmt6"
        log: 
            "logs/annotation/nr_blastx_parallel_{index}.log"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        params:
            db = config["trinotate"]["nr_path"],
            outfmt = config["trinotate"]["outfmt"],
            evalue = config["trinotate"]["e-value"]
        threads: config["parallel_run"]["threads"]
        resources: mem_mb=config["parallel_run"]["memory"]
        benchmark: 
            "benchmarks/annotation/nr_blastx_parallel_{index}.txt"
        shell:
            """
            diamond blastx --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --max-target-seqs 1 --threads {threads} --evalue {params.evalue} &> {log}
            """

    rule merge_nr_blastx:
        """
        Merge all NR blastx hits results.
        """
        input: 
            expand(OUTDIR + "/annotation/nr_blastx/chunks/transcriptome_nt_nr_diamond.part_{index}.outfmt6", index = indexes)
        output: 
            OUTDIR + "/annotation/nr_blastx/transcriptome_blastx_nr_diamond.outfmt6",
        log: 
            "logs/annotation/nr_blastx_merge.log"
        threads: 1
        resources: mem_mb=1000
        shell:  
            """
            cat {input} >> {output} 2>> {log}
            """

    rule load_nr_blastx:
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/nr_blastx/transcriptome_blastx_nr_diamond.outfmt6"
        output:
            touch(OUTDIR + "/annotation/nr_blastx/nr_blastx_load.done")
        conda: 
            "../envs/trinotate_v4.0.2.yaml"        
        log:
            "logs/annotation/nr_blastx_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
            db_name = "NR"
        threads: 1
        resources: mem_mb = 2000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_custom_blast {input.hits} --blast_type blastx --custom_db_name {params.db_name} &> {log}
            """
if "nr_blastp" in config["database_run"]:
    rule parallel_nr_blastp:
        """
        Run NR blastp in parallel for each chunk.
        """
        input:
            OUTDIR + "/transcriptome/split/prot/transcriptome_prot_modified.part_{index}.fasta"
        output: 
            OUTDIR + "/annotation/nr_blastp/chunks/transcriptome_prot_modified_nr_diamond.part_{index}.outfmt6"     
        log: 
            "logs/annotation/nr_blastp_parallel_{index}.log"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"  
        params:
            db = config["trinotate"]["nr_path"],
            outfmt = config["trinotate"]["outfmt"],
            evalue = config["trinotate"]["e-value"]
        threads: config["parallel_run"]["threads"]
        resources: mem_mb=config["parallel_run"]["memory"]
        benchmark: 
            "benchmarks/annotation/nr_blastp_parallel_{index}.txt"
        shell:
            """
            diamond blastp --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --max-target-seqs 1 --threads {threads} --evalue {params.evalue} &> {log}
            """

    rule merge_nr_blastp:
        """
        Merge all NR blastx hits results.
        """
        input: 
            expand(OUTDIR + "/annotation/nr_blastp/chunks/transcriptome_prot_modified_nr_diamond.part_{index}.outfmt6", index = indexes)
        output: 
            OUTDIR + "/annotation/nr_blastp/transcriptome_blastp_nr_diamond.outfmt6",
        log: 
            "logs/annotation/nr_blastp_merge.log"
        threads: 1
        resources: mem_mb=1000
        shell:  
            """
            cat {input} >> {output} 2>> {log}
            """

    rule load_nr_blastp:
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/nr_blastp/transcriptome_blastp_nr_diamond.outfmt6"
        output:
            touch(OUTDIR + "/annotation/nr_blastp/nr_blastp_load.done")
        conda: 
            "../envs/trinotate_v4.0.2.yaml"        
        log:
            "logs/annotation/nr_blastp_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
            db_name = "NR"
        threads: 1
        resources: mem_mb = 2000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_custom_blast {input.hits} --blast_type blastp --custom_db_name {params.db_name} &> {log}
            """

if "nt" in config["database_run"]:
    rule parallel_nt:
        input: 
            OUTDIR + "/transcriptome/split/nucl/transcriptome_nt.part_{index}.fasta"
        output: 
            OUTDIR + "/annotation/nt_blastn/chunks/transcriptome_blastn_nt.part_{index}.outfmt6"
        conda: 
            "../envs/trinotate_v4.0.2.yaml"
        log: 
            "logs/annotation/nt_blastn_parallel.{index}.log"
        params:
            db = config["trinotate"]["nt_path"],
            outfmt = config["trinotate"]["outfmt"],
            evalue = config["trinotate"]["e-value"] 
        threads: 
            config["parallel_run"]["threads"]
        resources: 
            mem_mb = config["parallel_run"]["memory"]
        benchmark: 
            "benchmarks/annotation/nt_blastn_parallel_{index}.txt"
        shell:
            """
            blastn -db {params.db} -query {input} -outfmt "{params.outfmt}" -out {output} -max_hsps 1 -max_target_seqs 1 -num_threads {threads} -evalue {params.evalue} &> {log}
            """

    rule merge_nt:
        """
        Merge all NT hits results.
        """
        input: 
            expand(OUTDIR + "/annotation/nt_blastn/chunks/transcriptome_blastn_nt.part_{index}.outfmt6", index = indexes)
        output: 
            OUTDIR + "/annotation/nt_blastn/transcriptome_blastn_nt.outfmt6",
        log: 
            "logs/annotation/nt_merge.log"
        threads:  1
        resources:  mem_mb=1000
        shell:  
            """
            cat {input} >> {output} 2>{log}
            """

    rule load_nt_hits:
        input:
            check = OUTDIR + "/annotation/trinotate_init.done",
            hits = OUTDIR + "/annotation/nt_blastn/transcriptome_blastn_nt.outfmt6"
        output:
            touch(OUTDIR + "/annotation/nt_blastn/nt_trinotate_load.done")
        conda: 
            "../envs/trinotate_v4.0.2.yaml"        
        log:
            "logs/annotation/nt_trinotate_load.log"
        params:
            path = config["trinotate"]["path"],
            sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite",
            db_name = "NT"
        threads: 1
        resources:  mem_mb=2000
        shell:
            """
            {params.path}/Trinotate --db {params.sqlite} --LOAD_custom_blast {input.hits} --blast_type blastx --custom_db_name {params.db_name} &> {log}
            """

rule get_report:
    input:
        conditional_outputs
    output:
        OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report.xls"
    log:
        "logs/annotation/" + config["PROJECT"] + "_trinotate_report.log"
    conda:
        "../envs/trinotate_v4.0.2.yaml"
    params:
        path = config["trinotate"]["path"],
        sqlite = OUTDIR + "/annotation/" + config["PROJECT"] + "_myTrinotate.sqlite"
    threads: 2
    resources:  mem_mb = 2000
    shell:
        """
        {params.path}/Trinotate --db {params.sqlite} --report > {output} 2>> {log}
        """
        
rule get_summary:
    input:
        annotation = OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report.xls"
    output:
        summary_table = OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report.summary",
        go_slims = OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report_trans.goslims",
        GO_plots = expand(OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_{category}_GOSlims.{ext}", category = ["ALL", "Cellular_Component", "Molecular_Function", "Biological_Process"], ext = config["format"])
    log:
        "logs/trinotate_report.log"
    conda:
        "../envs/trinotate_v4.0.2.yaml"
    params:
        path = config["trinotate"]["path"],
        prefix = config["PROJECT"], 
        fmt = config["format"], # plot file extention
        outdir = OUTDIR + "/annotation"
    threads:    1
    resources:  mem_mb = 2000
    shell:
        """
        # Summary table
        {params.path}/util/count_table_fields.pl {input.annotation} > {output.summary_table} 2> {log}
        # extract go terms - include ancestral terms (-I), gene (-G), transcripts (--trans, -T)
        resources/extract_GO_assignments_from_Trinotate_xls_updated.pl --Trinotate_xls {input.annotation} --trans > {params.outdir}/{params.prefix}_trans.goterms 2>> {log}
        # Get GOSlims 
        {params.path}/util/gene_ontology/Trinotate_GO_to_SLIM.pl  {params.outdir}/{params.prefix}_trans.goterms > {output.go_slims} 2>> {log}
        # Create GO plots
        Rscript scripts/GOSlim_plots.R {output.go_slims} {params.prefix} {params.fmt} {params.outdir}/plots/ 2&>> {log}
        """    

if "eggnog_mapper" in config["database_run"]:
    rule get_eggnog_plots:
        input:
            annotation = OUTDIR + "/annotation/" + config["PROJECT"] + "_trinotate_report.xls"
        output:
            KEGG_plots = expand(OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_KEGG_{name}dist." + config["format"], name = ["", "slim_"]),
            COG_plot = OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_COG_dist." + config["format"]
        log:
            "logs/trinotate_report.log"
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        params:
            prefix = config["PROJECT"], 
            kegg_file = "resources/kegg_df.csv",
            cog_file = "resources/cog_fun.txt",
            fmt = config["format"], # plot file extention
            outdir = OUTDIR + "/annotation"
        threads:    1
        resources:  mem_mb = 2000
        shell:
            """
            # KEGG
            Rscript scripts/kegg_plot.R {input.annotation} {params.prefix} {params.kegg_file} {params.fmt} {params.outdir}/plots/ 2&>> {log}
            # COG
            Rscript scripts/cog_plot.R {input.annotation} {params.cog_file} {params.prefix}  {params.fmt} {params.outdir}/plots/ 2&>> {log}
            """ 

if "nr_blastx" in config["database_run"]:
    rule species_plot:
        input:
            outfmt = OUTDIR + "/annotation/nr_blastx/transcriptome_blastx_nr_diamond.outfmt6"
        output:
            species_plot = OUTDIR + "/annotation/plots/" + config["PROJECT"] + "_top10_species_dist_NR_blastx." + config["format"]
        log:
            "logs/species_plot.log"
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        params:
            prefix = config["PROJECT"], 
            fmt = config["format"], # plot file extention
            outdir = OUTDIR + "/annotation"
        threads:    1
        resources:  mem_mb = 2000
        shell:
            """
            # Species - NR BLASTX
            Rscript scripts/species_plot.R {input.outfmt} {params.prefix} {params.fmt} {params.outdir}/plots/ 2&>>{log}
            """  
