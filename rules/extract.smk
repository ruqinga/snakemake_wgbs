rule extract_methylation:
    input:
        deduplicated_bam = get_dedu_out
    output:
        cytosine_report = (
            "{extract_out}/{sample}/{sample}_trimmed_bismark_bt2_se.deduplicated.bedGraph.gz"
            if config["dt"] == "SE"
            else
            "{extract_out}/{sample}/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
        )
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        genome_folder = config["bismark"]["index"],
        extract_out = directories["extract_out"],
        output_folder = "{extract_out}/{sample}",
        option = config["bis_extractor"]["params"]
    log:
        log="{extract_out}/logs/{sample}.log"
    shell:
        """
        if [ "{config[dt]}" == "SE" ]; then
            bismark_methylation_extractor {params.option} \
                --genome_folder {params.genome_folder} \
                {input.deduplicated_bam} \
                -o {params.output_folder} > {log.log} 2>&1
        else
            bismark_methylation_extractor --paired-end {params.option} \
                --genome_folder {params.genome_folder} \
                {input.deduplicated_bam} \
                -o {params.output_folder} > {log.log} 2>&1
        fi
        """