import os
import glob

configfile: "config/config.yaml"

rule all:
    input: 
        "resources/TRINOTATE_DATA_DIR/trinotate_init.done",
        "resources/signalp6/signalp6_register.done" if "signalp6" in config["database_run"] else [],
        "resources/tmhmm2/tmhmm2_register.done" if "tmhmmv2" in config["database_run"] else []

rule create_sqlite_db:
    "Initial creation of the Trinotate sqlite database and downloading of the required data sets"
    input:
        "resources/Trinotate-Trinotate-v4.0.2/Trinotate"
    output:
        #directory("data/TRINOTATE_DATA_DIR")
        touch("resources/TRINOTATE_DATA_DIR/trinotate_init.done")
    log:
        "logs/annotation/trinotate_create_sqlite_db.log"
    conda: 
        "../envs/trinotate_v4.0.2.yaml"
    params:
        data_dir = "resources/TRINOTATE_DATA_DIR",
        db = "resources/myTrinotate.sqlite"
    threads: 4
    resources:  mem_mb = 8000
    shell: 
        """
        {input} --db {params.db} --create --trinotate_data_dir {params.data_dir} 2> {log}
        """

if "signalp6" in config["database_run"]:
    rule register_signalp6:
        """
        Initializing signalp6. You should have been install the TAR.GZ file before!
        """
        input:
            config["signalp6_register"]
        output:
            touch("resources/signalp6/signalp6_register.done")
        log:
            "logs/annotation/register_signalp6.log"
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads:
            1
        resources:
            mem_mb=1000
        shell:
            """
            signalp6-register {input}
            echo "Signalp6 was registered" &>> {log}
            """

if "tmhmmv2" in config["database_run"]:
    rule register_tmhmm2:
        """
        Initializing tmhmmv2. You should have been download TAR.GZ file before!
        """
        input:
            config["tmhmmv2_register"] 
        output:
            touch("resources/tmhmm2/tmhmm2_register.done")
        log:
            "logs/annotation/register.tmhmm2.log"
        conda:
            "../envs/trinotate_v4.0.2.yaml"
        threads: 1
        resources:  mem_mb=1000
        shell:
            """
            tmhmm2-register {input}
            echo "TMHMMV2 was registered" &>> {log}
            touch resources/tmhmm2/tmhmm2_register.done
            """
