# 解析input，生成output
class SampleProcessor:
    def __init__(self, config):
        self.reads = config.get("reads", [])
        self.threshold = config.get("bw_threshold",3)
        self.sample_names = []
        self.sample_info = {}
        self.process_reads()

    # 生成sample和对应sample_info
    def process_reads(self):
        for read in self.reads:
            sample_name = read["read1"].split("/")[-1].replace("_1.fastq.gz", "").replace(".fastq.gz", "")
            if sample_name in self.sample_names:
                print(f"Error: Duplicate sample name '{sample_name}' found.")
                exit(1)
            self.sample_names.append(sample_name)

            # 判断是否包含read2
            if "read2" in read and read["read2"]:
                self.sample_info[sample_name] = "PE"
            else:
                self.sample_info[sample_name] = "SE"

    # 根据样品类型（SE/PE）确定trim之后得到的后缀
    def get_trim_ext(self, sample):
        return "1_val_1" if self.sample_info[sample] == "PE" else "trimmed"

    # 所有要生成的文件
    def generate_targets(self, sample):
        dt = 'pe' if self.sample_info[sample] == 'PE' else 'se'
        t = self.threshold
        return [
            f"Results/02_cleandata/trim_galore/{sample}_{self.get_trim_ext(sample)}.fq.gz",
            f"Results/03_qc/rawdata/multiqc_report.html",
            f"Results/03_qc/cleandata/multiqc_report.html",
            f"Results/04_bismark/{sample}_{self.get_trim_ext(sample)}_bismark_bt2_{dt}.bam",
            f"Results/05_dedu/{sample}_{self.get_trim_ext(sample)}_bismark_bt2_{dt}.deduplicated.sorted.bam",
            f"Results/06_extract/{sample}/{sample}_{self.get_trim_ext(sample)}_bismark_bt2_{dt}.deduplicated.bismark.cov.gz",
            f"Results/07_visualization/bw/{sample}_t{t}.bw",
            f"Results/08_tss_tes/repeats/{sample}_level_dis_mean.txt",
            f"Results/summary.csv"
        ]

    # 获取所有目标路径
    def get_all_targets(self):
        all_paths = []
        for sample in self.sample_names:
            paths = self.generate_targets(sample)
            all_paths.extend(paths)
        return all_paths


# Create input file list based on configuration
def get_fq_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"{config['fq_dir']}/{wildcards.sample}.fastq.gz"
    elif sample_info[wildcards.sample] == "PE":
        return [
            f"{config['fq_dir']}/{wildcards.sample}_1.fastq.gz",
            f"{config['fq_dir']}/{wildcards.sample}_2.fastq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_trimmed_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/02_cleandata/trim_galore/{wildcards.sample}_trimmed.fq.gz"
    elif sample_info[wildcards.sample] == "PE":
        return [
            f"Results/02_cleandata/trim_galore/{wildcards.sample}_1_val_1.fq.gz",
            f"Results/02_cleandata/trim_galore/{wildcards.sample}_2_val_2.fq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_cutted_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/02_cleandata/cut/{wildcards.sample}_trimmed.fq.gz"
    elif sample_info[wildcards.sample] == "PE":
        return [
            f"Results/02_cleandata/cut/{wildcards.sample}_1_val_1.fq.gz",
            f"Results/02_cleandata/cut/{wildcards.sample}_2_val_2.fq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_bismark_out(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/04_bismark/{wildcards.sample}_trimmed_bismark_bt2_se.bam"
    elif sample_info[wildcards.sample] == "PE":
        return f"Results/04_bismark/{wildcards.sample}_1_val_1_bismark_bt2_pe.bam"
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_dedu_out(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/05_dedu/{wildcards.sample}_trimmed_bismark_bt2_se.deduplicated.bam"
    elif sample_info[wildcards.sample] == "PE":
        return f"Results/05_dedu/{wildcards.sample}_1_val_1_bismark_bt2_pe.deduplicated.bam"
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")
