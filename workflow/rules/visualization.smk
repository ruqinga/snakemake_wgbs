rule sort_bed:
    input:
        cytosine_result="Results/06_extract/{sample}/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    output:
        sorted_bed = "Results/06_extract/{sample}/sorted.bed"
    conda:
        config["conda_env"]
    group: "processing_group"
    shell:
        """
        # sort 共6列 chr start end level methy unmethy
        zcat {input.cytosine_result} | LC_COLLATE=C sort -k1,1 -k2,2n > {output.sorted_bed} 
        """


rule bw:
    input:
        sorted_bed = "Results/06_extract/{sample}/sorted.bed"
    output:
        filtered_bed = "Results/07_visualization/bed/{sample}_t{t}_forbw.bed",
        filtered_bw = "Results/07_visualization/bw/{sample}_t{t}.bw"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        genome_fasta_fai = config["bismark"]["genome_fasta_fai"],
        extract_threshold = config["bw_threshold"],
        scripts_dir = "./workflow/scripts/filter_bed.sh"
    log:
        log="Results/07_visualization/logs/{sample}_{t}.log"
    shell:
        """
        # 筛选total count
        bash {params.scripts_dir} {params.genome_fasta_fai} {input.sorted_bed} {output.filtered_bed} {params.extract_threshold} > {log.log} 2>&1
     
        # 转为bw
        /home_data/home/slst/leixy2023/software/bedGraphToBigWig {output.filtered_bed} {params.genome_fasta_fai} {output.filtered_bw} > {log.log} 2>&1
        """

