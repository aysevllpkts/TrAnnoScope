# Aysevil Pektas
# 11-04-2024
# quality assessment of the contigs

##################### NOTES #####################
# add options for only prot and only nucl   ++  # 
# add busco plot                            ++  #
# add blastp for FL-representation          ++  #
# add length distribution for protein       ++  #
# add protein completeness plot                 #
#################################################


import os
import glob


configfile: "config/config.yaml"

OUTDIR = "output_" + config["PROJECT"] 

##------------------------- Define your input file ---------------------------##

NUCL = config["input_for_quality"]["nucl"]
PROT = config["input_for_quality"]["prot"]

if NUCL != "":
    if os.path.isfile(NUCL):
        print("Your nucleotide input defined in the config file will be used for quality assessment.")
    elif  NUCL == "None": 
        print("Quality assessment for nucleotide sequences of transcriptome will be skipped.")
    else:
        raise ValueError(f"Error: The input file for nucleotide sequences provided in the config file does not exist: {NUCL}")
else:
    print("Input for nucleotide sequences created by annotation step will be used for quality assessment.")
    NUCL = OUTDIR + "/transcriptome/transcriptome_nt.fasta"

if PROT != "":
    if os.path.isfile(PROT): 
        print("Your protein input defined in the config file will be used for quality assessment.")
    elif PROT == "None":
        print("Quality assessment for protein sequences of transcriptome will be skipped.")
    else:
        raise ValueError(f"Error: The input file for protein sequences provided in the config file does not exist: {PROT}")
else:
    print("Input for protein sequences created by annotatıon step will be used for quality assessment.")
    PROT = OUTDIR + "/transcriptome/transcriptome_prot_modified.fasta"

# Check if the TRINOTATE_DATA_DIR directory exists for SwissProt DB
if os.path.isdir("resources/TRINOTATE_DATA_DIR"):
    swissprot_diamond = "resources/TRINOTATE_DATA_DIR/uniprot_sprot.dmnd"
    swissprot_fasta = "resources/TRINOTATE_DATA_DIR/uniprot_sprot.pep"
else:
    swissprot_diamond = config['swissprot']["diamond"]
    swissprot_fasta = config['swissprot']["fasta"] 

## --------------------------------------------------------------------------- ##

rule all:
    input:
        OUTDIR + "/quality_assessment/stats/transcriptome_nt_NanoStats.txt" if NUCL != "None" else [],
        OUTDIR + "/quality_assessment/busco/"+config["busco"]["lineage"]+"_transcriptome/short_summary.specific."+config["busco"]["lineage"]+"."+config["busco"]["lineage"]+"_transcriptome.txt" if NUCL != "None" else [],
        OUTDIR + "/quality_assessment/busco/"+config["busco"]["lineage"]+"_transcriptome/busco_figure.png" if NUCL != "None" else [],
        OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastx.outfmt6.coverage.hist" if NUCL != "None" else [],
        OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastx.outfmt6.coverage.hist" if NUCL != "None" and config['custom_db']["diamond"] != "" else [],
        OUTDIR + "/quality_assessment/stats/transcriptome_prot_summary.txt" if PROT != "None" else [],
        OUTDIR + "/quality_assessment/stats/transcriptome_prot_dist_plot.png" if PROT != "None" else [],
        OUTDIR + "/quality_assessment/stats/protein_completeness.done" if PROT != "None" else [],
        OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastp.outfmt6.coverage.hist" if PROT != "None" else [],
        OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastp.outfmt6.coverage.hist" if PROT != "None" and config['custom_db']["diamond"] != "" else []


## --- RULES --- ##

if NUCL != "None":
    rule nanoplot:
        "Summary tool for long sequencing data."
        input:
            NUCL
        output:
            OUTDIR + "/quality_assessment/stats/transcriptome_nt_NanoStats.txt"
        log:
            "logs/quality_assessment/nanoplot.log"
        conda:
            "../envs/quality_assessment.yaml"
        params:
            prefix = "transcriptome_nt_",
            out_dir = OUTDIR + "/quality_assessment/stats"
        threads:
            config["nanoplot"]["threads"]
        resources:
            mem_mb = config["nanoplot"]["memory"]
        shell:
            """
            mkdir -p {params.out_dir} &> {log}
            NanoPlot -o {params.out_dir} -p {params.prefix} -t {threads} --N50 --fasta {input} &>> {log}
            """

    rule busco:
        """
        Assessing the completeness of a transcriptome by comparing it to a set of BUSCO genes.
        """
        input: 
            NUCL
        output:
            OUTDIR + "/quality_assessment/busco/"+config["busco"]["lineage"]+"_transcriptome/short_summary.specific."+config["busco"]["lineage"]+"."+config["busco"]["lineage"]+"_transcriptome.txt"
        log:
            "logs/quality_assessment/busco.log"
        conda:
            "../envs/quality_assessment.yaml"
        params:
            lineage = config["busco"]["lineage"],
            outdir = OUTDIR + "/quality_assessment/busco",
            data_path = config["busco"]["data_path"]
        threads:
            config["busco"]["threads"]
        resources:
            mem_mb = config["busco"]["memory"]
        shell:
            """
            mkdir -p {params.outdir}
            if [ -z {params.lineage} ] &> {log}
            then
                busco -f -i {input} -o {params.lineage}_transcriptome --out_path {params.outdir} --auto-lineage -m trans -c {threads} &>> {log}
            else
                busco -f -i {input} -o {params.lineage}_transcriptome --out_path {params.outdir} -l {params.lineage} -m trans -c {threads} --download_path {params.data_path} &>> {log}
            fi     
            """

    rule generate_busco_plot:
        """
        Generate BUSCO plot.
        """
        input:
            OUTDIR + "/quality_assessment/busco/"+config["busco"]["lineage"]+"_transcriptome/short_summary.specific."+config["busco"]["lineage"]+"."+config["busco"]["lineage"]+"_transcriptome.txt"
        output:
            OUTDIR + "/quality_assessment/busco/"+config["busco"]["lineage"]+"_transcriptome/busco_figure.png"
        log:
            "logs/quality_assessment/busco_plot.log" 
        params:
            dir = OUTDIR + "/quality_assessment/busco/"+config["busco"]["lineage"]+"_transcriptome"
        conda:
            "../envs/quality_assessment.yaml"
        threads:
            1
        resources: 
            mem_mb = 1000
        shell:
            """
            generate_plot.py -wd {params.dir}
            """

    rule blastx_uniprot:
        "Blast transcriptome against Uniprot for highly similar sequences."
        input:
            NUCL
        output: 
            OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastx.outfmt6"
        log:
            "logs/quality_assessment/fl_blastx.log"
        conda:
            "../envs/quality_assessment.yaml"
        params:
            outfmt = config["diamond"]["outfmt"],
            db = swissprot_diamond,
            extra = config["diamond"]["extra_parameters"],
            dir = OUTDIR + "/quality_assessment/fl_representation"
        threads:
            config["diamond"]["threads"]
        resources:
            mem_mb = config["diamond"]["memory"]
        shell:
            """
            mkdir -p {params.dir}
            diamond blastx --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --threads {threads} {params.extra} &> {log}       
            """

    rule fl_summary_uniprot_blastx:
        "Coverage summary of the blastx results."
        input:
            OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastx.outfmt6",
        output:
            OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastx.outfmt6.coverage.hist"
        log:
            "logs/quality_assessment/fl_summary_uniprot.log"
        threads:
            1
        resources:
            mem_mb = 1000
        shell:
            """
            bash scripts/coverage_table.sh {input} >> {output} 2> {log}
            """

    if config["custom_db"]["diamond"] != "":
        rule blastx_custom:
            "Blast transcriptome against custom database - specific to your data - for highly similar sequences."
            input:
                NUCL
            output: 
                OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastx.outfmt6"
            log:
                "logs/quality_assessment/fl_blastx_custom.log"
            conda:
                "../envs/quality_assessment.yaml"
            params:
                outfmt = config["diamond"]["outfmt"],
                db = config["custom_db"]["diamond"],
                extra = config["diamond"]["extra_parameters"]
            threads:
                config["diamond"]["threads"]
            resources:
                mem_mb = config["diamond"]["memory"]
            shell:
                """
                diamond blastx --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --threads {threads} {params.extra} &> {log}       
                """

        rule fl_summary_custom_blastx:
            "Coverage summary of the blastx results."
            input:
                outfmt = OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastx.outfmt6",
            output:
                OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastx.outfmt6.coverage.hist"
            log:
                "logs/quality_assessment/fl_summary_customDB.log"
            threads:
                1
            resources:
                mem_mb = 100
            shell:
                """
                bash scripts/coverage_table.sh {input} >> {output} 2> {log}
                """
else:
   print("Steps for nucleotide sequence have been skipped.")

if PROT != "None":
    rule summary_stats_prot:
        """
        Create Summary Stats TXT file from protein FASTA file. 
        """
        input:
            PROT
        output:
            OUTDIR + "/quality_assessment/stats/transcriptome_prot_summary.txt"
        log:
            "logs/quality_assessment/prot_summary_stats.log"
        conda:
            "../envs/quality_assessment.yaml"
        threads:
            2
        resources:  mem_mb = 1000
        shell:
            "python3 scripts/summary_stats_prot.py {input} {output} {log}"

    rule plot_dist_prot:
        """
        Create a length distribution histogram file from protein FASTA file. 
        """
        input:
            PROT
        output:
            OUTDIR + "/quality_assessment/stats/transcriptome_prot_dist_plot.png"
        log:
            "logs/quality_assessment/prot_dist_plot.log"
        conda:
            "../envs/quality_assessment.yaml"
        threads:
            2
        resources:  mem_mb = 1000
        shell:
            """
            Rscript scripts/prot_dist_hist.R {input} {output} {log}
            """  

    rule protein_completeness:
        """
        Generate protein completeness plot from protein FASTA file.
        """
        input:
            PROT
        output:
            OUTDIR + "/quality_assessment/stats/protein_completeness.done"
        log:
           "logs/quality_assessment/prot_completeness.log" 
        params:
            dir = OUTDIR + "/quality_assessment/stats"
        conda:
           "../envs/quality_assessment.yaml" 
        threads:
            2
        resources: mem_mb = 2000
        shell:
            """
            Rscript scripts/protein_completeness.R {input} {params.dir} {log}
            touch {output}
            """

    rule blastp_uniprot:
        "Protein blast transcriptome against Uniprot for highly similar sequences."
        input:
            PROT
        output: 
            OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastp.outfmt6"
        log:
            "logs/quality_assessment/fl_blastp.log"
        conda:
            "../envs/quality_assessment.yaml"
        params:
            outfmt = config["diamond"]["outfmt"],
            db = swissprot_diamond,
            extra = config["diamond"]["extra_parameters"],
            dir = OUTDIR + "/quality_assessment/fl_representation"
        threads:
            config["diamond"]["threads"]
        resources:
            mem_mb = config["diamond"]["memory"]
        shell:
            """
            mkdir -p {params.dir}
            diamond blastp --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --threads {threads} {params.extra} &> {log}       
            """

    rule fl_summary_uniprot_blastp:
        "Coverage summary of the blastx results."
        input:
            OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastp.outfmt6",
        output:
            OUTDIR + "/quality_assessment/fl_representation/transcriptome_uniprot_blastp.outfmt6.coverage.hist"
        log:
            "logs/quality_assessment/fl_summary_uniprot.log"
        threads:
            1
        resources:
            mem_mb = 1000
        shell:
            """
            bash scripts/coverage_table.sh {input} >> {output} 2> {log}
            """

    if config["custom_db"]["diamond"] != "":
        rule blastp_custom:
            "Protein blast transcriptome against custom database - specific to your data - for highly similar sequences."
            input:
                PROT
            output: 
                OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastp.outfmt6"
            log:
                "logs/quality_assessment/fl_blastx_custom.log"
            conda:
                "../envs/quality_assessment.yaml"
            params:
                outfmt = config["diamond"]["outfmt"],
                db = config["custom_db"]["diamond"],
                extra = config["diamond"]["extra_parameters"]
            threads:
                config["diamond"]["threads"]
            resources:
                mem_mb = config["diamond"]["memory"]
            shell:
                """
                diamond blastp --db {params.db} --query {input} --outfmt {params.outfmt} --out {output} --threads {threads} {params.extra} &> {log}       
                """

        rule fl_summary_custom_blastp:
            "Coverage summary of the blastx results."
            input:
                outfmt = OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastp.outfmt6",
            output:
                OUTDIR + "/quality_assessment/fl_representation/transcriptome_customDB_blastp.outfmt6.coverage.hist"
            log:
                "logs/quality_assessment/fl_summary_customDB.log"
            threads:
                1
            resources:
                mem_mb = 100
            shell:
                """
                bash scripts/coverage_table.sh {input} >> {output} 2> {log}
                """
else:
   print("Steps for protein sequence have been skipped.") 
