# visualization 结束之后运行summary
rule summary:
    input:
        extract_finish = expand("Results/06_extract/{sample}/sorted.bed", sample = samples)
    output:
        summary = "Results/summary.csv"
    params:
        trim_log_dir="Results/02_cleandata/trim_galore/logs/", # resource of raw reads
        bismark_report_dir="Results/04_bismark/", # resource of clean and mapped(aligned) reads
        dedu_report_dir="Results/05_dedu/", # resource of dedu reads
        extract_methy_logs_dir="Results/06_extract/logs/", # resource of dedu_genome_methy_ratio
        bed_dir="Results/06_extract/"
    conda:
        config["conda_env"]
    group: "global_process"
    script:
        "../scripts/summary.py"
