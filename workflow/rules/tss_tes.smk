rule tss_tes_repeats:
    input:
        sorted_bed="Results/06_extract/{sample}/sorted.bed"
    output:
        sorted_v2_bedGraph=temp("Results/06_extract/{sample}/sorted_v2.bedGraph"),
        intersect=temp("Results/08_tss_tes/repeats/{sample}_intersect.bed"),
        norm_dis=temp("Results/08_tss_tes/repeats/{sample}_level_dis.txt"),
        norm_mean="Results/08_tss_tes/repeats/{sample}_level_dis_mean.txt"
    conda:
        config["conda_env"]
    group: "Additional_analysis"
    params:
        index_bed=config["tss_tes"]["repeats"],
        bw_threshold=config["bw_threshold"],
        script="./workflow/scripts/extract_tss_tes.sh",
        mean_py="./workflow/scripts/mean.py"
    log:
        log="Results/08_tss_tes/logs/{sample}.log"
    shell:
        """
        bash {params.script} \
            {input.sorted_bed} \
            {output.sorted_v2_bedGraph} \
            {params.index_bed} \
            {output.intersect} \
            {params.bw_threshold} \
            {output.norm_dis} \
            {output.norm_mean} \
            {params.mean_py} \
            {log.log}
        """
