rule trim_pe:
    input:
        read = get_fq_list
    output:
        trimmed_read = [
                "Results/02_cleandata/trim_galore/{sample}_1_val_1.fq.gz",
                "Results/02_cleandata/trim_galore/{sample}_2_val_2.fq.gz"
        ]
    conda:
        config["conda_env"]
    params:
        option = config["trim"]["params"],
        trim_out= "Results/02_cleandata/trim_galore"
    log:
        log = "Results/02_cleandata/trim_galore/logs/{sample}.log"
    shell:
        """
           trim_galore {params.option} --paired {input.read[0]} {input.read[1]} -o {params.trim_out} > {log.log} 2>&1
        """

rule trim_se:
    input:
        read = get_fq_list
    output:
        trimmed_read = "Results/02_cleandata/trim_galore/{sample}_trimmed.fq.gz"
    conda:
        config["conda_env"]
    params:
        option = config["trim"]["params"],
        trim_out= "Results/02_cleandata/trim_galore"
    log:
        log = "Results/02_cleandata/trim_galore/logs/{sample}.log"
    shell:
        """
           trim_galore {params.option} {input.read} -o {params.trim_out} > {log.log} 2>&1
        """