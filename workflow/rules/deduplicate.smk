rule deduplicate_bismark_pe:
    input:
        bam = get_bismark_out
    output:
        deduplicated_bam = temp("Results/05_dedu/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bam"),
        deduplicated_sorted_bam = "Results/05_dedu/{sample}_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        dedu_out = "Results/05_dedu"
    log:
        log="Results/05_dedu/logs/{sample}.log"
    shell:
        """
        deduplicate_bismark --bam {input.bam} --output_dir {params.dedu_out} > {log.log} 2>&1
        # sort
        samtools sort -@ 20 -o {output.deduplicated_sorted_bam} {output.deduplicated_bam}
        """

rule deduplicate_bismark:
    input:
        bam = get_bismark_out
    output:
        deduplicated_bam = temp("Results/05_dedu/{sample}_trimmed_bismark_bt2_se.deduplicated.bam"),
        deduplicated_sorted_bam = "Results/05_dedu/{sample}_trimmed_bismark_bt2_se.deduplicated.sorted.bam"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        dedu_out = "Results/05_dedu"
    log:
        log="Results/05_dedu/logs/{sample}.log"
    shell:
        """
        deduplicate_bismark --bam {input.bam} --output_dir {params.dedu_out} > {log.log} 2>&1
        # sort
        samtools sort -@ 20 -o {output.deduplicated_sorted_bam} {output.deduplicated_bam}
        """