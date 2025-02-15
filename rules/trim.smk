
rule trim:
    input:
        read = get_fq_list
    output:
        trimmed_read = (
            "{clean_out}/trim_galore/{sample}_trimmed.fq.gz"
            if config["dt"] == "SE"
            else [
                "{clean_out}/trim_galore/{sample}_1_val_1.fq.gz",
                "{clean_out}/trim_galore/{sample}_2_val_2.fq.gz"
            ]
        )
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["trim"]["params"],
        clean_out = directories["clean_out"],
        trim_out = f"{directories['clean_out']}/trim_galore"
    log:
        log = "{clean_out}/trim_galore/logs/{sample}.log"
    shell:
        """
        if [ "{config[dt]}" == "SE" ]; then
           trim_galore {params.option} {input.read} -o {params.trim_out} > {log.log} 2>&1
        else
           trim_galore {params.option} --paired {input.read[0]} {input.read[1]} -o {params.trim_out} > {log.log} 2>&1
        fi
        """
