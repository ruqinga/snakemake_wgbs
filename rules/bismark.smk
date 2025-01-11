rule bismark:
    input:
        trimmed_read = get_trimmed_list
    output:
        bam = (
            "{bis_out}/{sample}_trimmed_bismark_bt2_se.bam"
            if config["dt"] == "SE"
            else [
                "{bis_out}/{sample}_1_val_1_bismark_bt2_pe.bam"
            ]
        )
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["bismark"]["params"],
        genome = config["bismark"]["index"],
        bis_out = directories["bis_out"]
    log:
        log="{bis_out}/logs/{sample}.log"
    shell:
        """
        if [ "{config[dt]}" == "SE" ]; then
            bismark {params.option} --genome {params.genome} {input.trimmed_read} -o {params.bis_out} > {log.log} 2>&1
        else
            bismark {params.option} --genome {params.genome} -1 {input.trimmed_read[0]} -2 {input.trimmed_read[1]} -o {params.bis_out} > {log.log} 2>&1
        fi
        """
