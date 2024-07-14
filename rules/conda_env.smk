
rule all:
    input:
        "logs/precheck/blast.done",
        "logs/precheck/blobtools2.done",
        "logs/precheck/bowtie2.done",
        "logs/precheck/busco.done",
        "logs/precheck/cd_hit_est.done",
        "logs/precheck/evigene.done",
        "logs/precheck/fmlrc.done",
        "logs/precheck/pacbio.done",
        "logs/precheck/qc.done",
        "logs/precheck/quality_assessment.done",
        "logs/precheck/trinotate.done"

        


rule blast:
    output: "logs/precheck/blast.done"
    conda:
        "../envs/blast.yaml"
    shell:
        " "
        
rule blobtools2:
    output: "logs/precheck/blobtools2.done"
    conda:
        "../envs/blobtools2.yaml"
    shell:
        " "

rule bowtie2:
    output: "logs/precheck/bowtie2.done"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        " "

rule busco:
    output: "logs/precheck/busco.done"
    conda:
        "../envs/busco.yaml"
    shell:
        " "

rule cd_hit_est:
    output: "logs/precheck/cd_hit_est.done"
    conda:
        "../envs/cd-hit-est.yaml"
    shell:
        " "

rule evigene:
    output: "logs/precheck/evigene.done"
    conda:
        "../envs/evigene.yaml"
    shell:
        " "

rule fmlrc:
    output: "logs/precheck/fmlrc.done"
    conda:
        "../envs/fmlrc.yaml"
    shell:
        " "

rule pacbio:
    output: "logs/precheck/pacbio.done"
    conda:
        "../envs/pacbio.yaml"
    shell:
        " "

rule qc:
    output: "logs/precheck/qc.done"
    conda:
        "../envs/qc.yaml"
    shell:
        " "

rule quality_assessment:
    output: "logs/precheck/quality_assessment.done"
    conda:  
        "../envs/quality_assessment.yaml"
    shell: " "   

rule trinotate:
    output: "logs/precheck/trinotate.done"
    conda:
        "../envs/trinotate_v4.0.2.yaml"
    shell:
        " "

