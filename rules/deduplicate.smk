rule deduplicate_bismark:
    input:
        bam = get_bismark_out
    output:
        deduplicated_bam = (
            "{dedu_out}/{sample}_trimmed_bismark_bt2_se.deduplicated.bam"
            if config["dt"] == "SE"
            else
            "{dedu_out}/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bam"
        )
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        dedu_out = directories["dedu_out"]
    log:
        log="{bis_out}/logs/{sample}.log"
    shell:
        """
        deduplicate_bismark --bam {input.bam} --output_dir {params.dedu_out} > {log.log} 2>&1
        """
