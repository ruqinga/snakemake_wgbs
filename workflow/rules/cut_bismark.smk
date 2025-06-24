# zhenglab使用的wgbs的kit里会额外加上一段序列，需要在trim之后再去除
rule cut_pe:
    input:
        trimmed_read = get_trimmed_list
    output:
        cutted_read = [
                "Results/02_cleandata/cut/{sample}_1_val_1.fq.gz",
                "Results/02_cleandata/cut/{sample}_2_val_2.fq.gz"
            ]
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        options_pe = config["cutadapt"]["pe"],
        clean_out = "Results/02_cleandata/cut"
    log:
        log="Results/02_cleandata/logs/{sample}.log"
    shell:
        """
        cutadapt {params.options_pe} -o {output.cutted_read[0]} -p {output.cutted_read[1]} {input.trimmed_read[0]} {input.trimmed_read[1]} > {log.log} 2>&1
        """

rule cut_se:
    input:
        trimmed_read = get_trimmed_list
    output:
        cutted_read = "Results/02_cleandata/cut/{sample}_trimmed.fq.gz"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        options_se = config["cutadapt"]["se"],
        clean_out = "Results/02_cleandata/cut"
    log:
        log="Results/02_cleandata/logs/{sample}.log"
    shell:
        """
        cutadapt {params.options_se} -o {output.cutted_read} {input.trimmed_read} > {log.log} 2>&1
        """

rule bismark_pe:
    input:
        cutted_read = get_cutted_list
    output:
        bam = "Results/04_bismark/{sample}_1_val_1_bismark_bt2_pe.bam"
    conda:
        config["conda_env"]
    params:
        option = config["bismark"]["params"],
        genome = config["bismark"]["index"],
        strategy= config["bismark_strategy"],
        bis_out = "Results/04_bismark"
    log:
        log="Results/04_bismark/logs/{sample}.log"
    shell:
        """
        bismark {params.option} {params.strategy} --genome {params.genome} -1 {input.cutted_read[0]} -2 {input.cutted_read[1]} -o {params.bis_out} > {log.log} 2>&1
        """


rule bismark_se:
    input:
        cutted_read = get_cutted_list
    output:
        bam = "Results/04_bismark/{sample}_trimmed_bismark_bt2_se.bam"
    conda:
        config["conda_env"]
    params:
        option = config["bismark"]["params"],
        genome = config["bismark"]["index"],
        strategy= config["bismark_strategy"],
        bis_out = "Results/04_bismark"
    log:
        log="Results/04_bismark/logs/{sample}.log"
    shell:
        """
        bismark {params.option} {params.strategy} --genome {params.genome} {input.cutted_read} -o {params.bis_out} > {log.log} 2>&1
        """

