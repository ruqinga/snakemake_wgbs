def get_directories(config):
    """
    根据配置文件中的 work_dir 和目录结构返回所有目录的路径。
    """
    work_dir = config["work_dir"]
    directories = {
        key: value.format(work_dir=work_dir) if "{work_dir}" in value else value
        for key, value in config["directories"].items()
    }
    # 打印调试信息
    print(f"Resolved directories: {directories}")
    return directories


def get_sample_name(filepath, dt):
    """
    从文件路径中提取样本名，基于数据类型（单端或双端）。
    """
    try:
        if dt == "PE":
            return filepath.split("/")[-1].replace("_1.fastq.gz", "").replace("_2.fastq.gz", "")
        elif dt == "SE":
            return filepath.split("/")[-1].replace(".fastq.gz", "")
        else:
            raise ValueError(f"Unsupported 'dt': {dt}")
    except Exception as e:
        raise ValueError(f"Error in get_sample_name for {filepath}: {e}")


def get_sample_list(config):
    """
    从配置文件中提取样本列表。
    """
    dt = config.get("dt")
    reads = config.get("reads", [])
    samples = {get_sample_name(read["read1"], dt) for read in reads}
    samples_list = list(samples)
    # 打印调试信息
    print("Samples:", samples_list)
    return samples_list


def get_all(directories, samples):
    """
    根据配置生成所有需要的目标文件路径列表。
    """
    dt = config.get("dt")
    all_targets = []

    if dt == "PE":
        all_targets += expand(
            "{clean_out}/trim_galore/{sample}_1_val_1.fq.gz",
            clean_out=directories["clean_out"],
            sample=samples
        )
        all_targets += expand(
            "{bis_out}/{sample}_1_val_1_bismark_bt2_pe.bam",
            bis_out=directories["bis_out"],
            sample=samples
        )
        all_targets += expand(
            "{dedu_out}/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bam",
            dedu_out=directories["dedu_out"],
            sample=samples
        )
        all_targets += expand(
            "{extract_out}/{sample}/{sample}_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz",
            extract_out=directories["extract_out"],
            sample=samples
        )

    elif dt == "SE":
        all_targets += expand(
            "{bis_out}/{sample}_trimmed_bismark_bt2_se.bam",
            bis_out=directories["bis_out"],
            sample=samples
        )
        all_targets += expand(
            "{dedu_out}/{sample}_trimmed_bismark_bt2_se.deduplicated.bam",
            dedu_out=directories["dedu_out"],
            sample=samples
        )
        all_targets += expand(
            "{extract_out}/{sample}/{sample}_trimmed_bismark_bt2_se.deduplicated.bedGraph.gz",
            extract_out=directories["extract_out"],
            sample=samples
        )


    else:
        raise ValueError(f"Unsupported 'dt': {dt}")

    # 打印调试信息
    print(f"Generated targets for dt={dt}: {all_targets}")
    return all_targets


# Create input file list based on configuration
def get_fq_list(wildcards):
    if config["dt"] == "SE":
        return f"{config['fq_dir']}/{wildcards.sample}.fastq.gz"
    elif config["dt"] == "PE":
        return [
            f"{config['fq_dir']}/{wildcards.sample}_1.fastq.gz",
            f"{config['fq_dir']}/{wildcards.sample}_2.fastq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_trimmed_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['trim_out']}/{wildcards.sample}_trimmed.fq.gz"
    elif config["dt"] == "PE":
        return [
            f"{directories['trim_out']}/{wildcards.sample}_1_val_1.fq.gz",
            f"{directories['trim_out']}/{wildcards.sample}_2_val_2.fq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_cutted_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['clean_out']}/cut/{wildcards.sample}_trimmed.fq.gz"
    elif config["dt"] == "PE":
        return [
            f"{directories['clean_out']}/cut/{wildcards.sample}_1_val_1.fq.gz",
            f"{directories['clean_out']}/cut/{wildcards.sample}_2_val_2.fq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_bismark_out(wildcards):
    if config["dt"] == "SE":
        return f"{directories['bis_out']}/{wildcards.sample}_trimmed_bismark_bt2_se.bam"
    elif config["dt"] == "PE":
        return f"{directories['bis_out']}/{wildcards.sample}_1_val_1_bismark_bt2_pe.bam"
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_dedu_out(wildcards):
    if config["dt"] == "SE":
        return f"{directories['dedu_out']}/{wildcards.sample}_trimmed_bismark_bt2_se.deduplicated.bam"
    elif config["dt"] == "PE":
        return f"{directories['dedu_out']}/{wildcards.sample}_1_val_1_bismark_bt2_pe.deduplicated.bam"
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")
