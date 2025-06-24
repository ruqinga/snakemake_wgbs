rule extract_methylation_pe:
    input:
        deduplicated_bam = get_dedu_out
    output:
        cytosine_result = "Results/06_extract/{sample}/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    conda:
        config["conda_env"]
    params:
        genome_folder = config["bismark"]["index"],
        option= config["bis_extractor"]["params"],
        output_folder = "Results/06_extract/{sample}"
    log:
        log="Results/06_extract/logs/{sample}.log"
    shell:
        """
        bismark_methylation_extractor --paired-end {params.option} \
            --genome_folder {params.genome_folder} \
            {input.deduplicated_bam} \
            -o {params.output_folder} > {log.log} 2>&1
        """

rule extract_methylation_se:
    input:
        deduplicated_bam = get_dedu_out
    output:
        cytosine_result = "Results/06_extract/{sample}/{sample}_trimmed_bismark_bt2_se.deduplicated.bismark.cov.gz"
    conda:
        config["conda_env"]
    params:
        genome_folder = config["bismark"]["index"],
        option= config["bis_extractor"]["params"],
        output_folder = "Results/06_extract/{sample}"
    log:
        log="Results/06_extract/logs/{sample}.log"
    shell:
        """
        bismark_methylation_extractor {params.option} \
            --genome_folder {params.genome_folder} \
            {input.deduplicated_bam} \
            -o {params.output_folder} > {log.log} 2>&1
        """