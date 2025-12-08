rule sort_bed:
    input:
        cytosine_result=get_cytosine_result
    output:
        sorted_bed = "Results/06_extract/{sample}/sorted.bed"
    conda:
        config["conda_env"]
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
        bedGraphToBigWig {output.filtered_bed} {params.genome_fasta_fai} {output.filtered_bw} > {log.log} 2>&1
        """

rule bin_bw:
    input:
        sorted_bed="Results/06_extract/{sample}/sorted.bed"
    output:
        filtered_bed="Results/07_visualization/bed/{sample}_bin_100_t{t}.bed",
        filtered_bw="Results/07_visualization/bw/{sample}_bin_100_t{t}.bw"
    conda:
        config["conda_env"]
    params:
        genome_fasta_fai=config["bismark"]["genome_fasta_fai"],
        extract_threshold=config["bin_threshold"],
        scripts_dir="./workflow/scripts/binning.py",
        out_dir="Results/07_visualization"
    log:
        log="Results/07_visualization/logs/{sample}_bin_100_t{t}.log"
    shell:
        """
        python {params.scripts_dir} -i {input.sorted_bed} -o {params.out_dir} --chrom_fai {params.genome_fasta_fai} --count_threshold {params.extract_threshold} --sorted > {log.log} 2>&1
        """