rule bismark_pe:
    input:
        trimmed_read = get_trimmed_list
    output:
        bam = "Results/04_bismark/{sample}_1_val_1_bismark_bt2_pe.bam"
    conda:
        config["conda_env"]
    params:
        genome = config["bismark"]["index"],
        strategy= config["bismark_strategy"],
        bis_out = "Results/04_bismark"
    log:
        log="Results/04_bismark/logs/{sample}.log"
    shell:
        """
        bismark {params.strategy} --parallel 8 --genome {params.genome} -1 {input.trimmed_read[0]} -2 {input.trimmed_read[1]} -o {params.bis_out} > {log.log} 2>&1
        """


rule bismark_se:
    input:
        trimmed_read = get_trimmed_list
    output:
        bam = "Results/04_bismark/{sample}_trimmed_bismark_bt2.bam"
    conda:
        config["conda_env"]
    params:
        genome = config["bismark"]["index"],
        strategy= config["bismark_strategy"],
        bis_out = "Results/04_bismark"
    log:
        log="Results/04_bismark/logs/{sample}.log"
    shell:
        """
        bismark {params.strategy} --parallel 8 --genome {params.genome} {input.trimmed_read} -o {params.bis_out} > {log.log} 2>&1
        """