import re
import logging
import pandas as pd
from pathlib import Path
import sys

# trim_logs_dir = Path("Results/02_cleandata/trim_galore/logs/")  # resource of raw reads
# bismark_report_dir = Path("Results/04_bismark/")  # resource of clean and mapped(aligned) reads
# dedu_report_dir = Path("Results/05_dedu/")  # resource of dedu reads
# extract_methy_logs_dir = Path("Results/06_extract/logs/")  # resource of dedu_genome_methy_ratio
# bed_dir = Path("Results/06_extract/")
# outputfile = Path("Results/summary.csv")

# 获取 Snakemake 参数
trim_logs_dir = Path(snakemake.params.trim_log_dir)
bismark_report_dir = Path(snakemake.params.bismark_report_dir)
dedu_report_dir = Path(snakemake.params.dedu_report_dir)
extract_methy_logs_dir = Path(snakemake.params.extract_methy_logs_dir)
bed_dir = Path(snakemake.params.bed_dir)
outputfile = Path(snakemake.output.summary)

# 设置日志格式
logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

def extract_from_trim_log(log_file):
    with open(log_file, 'r') as f:
        content = f.read()

    input_reads = re.search(r"Total number of sequences analysed:\s*(\d+)", content)
    if not input_reads:
        input_reads = re.search(r"(\d+) sequences processed in total", content)

    if not input_reads:
        logging.warning(f"{log_file} 没有找到 input_reads")
        return 0

    return int(input_reads.group(1))


def extract_bis_report(report_file):
    with open(report_file, 'r') as f:
        content = f.read()

    clean = re.search(r"Sequence pairs analysed in total:\s*(\d+)", content)
    aligned = re.search(r"Number of paired-end alignments with a unique best hit:\s*(\d+)", content)
    methy = re.search(r"C methylated in CpG context:\s*([\d.]+)", content)

    methy = float(methy.group(1)) * 0.01 # report里显示的是百分比

    if not clean or not aligned or not methy:
        logging.error(f"{report_file} 缺失bis配对相关字段，单端未实现")
        sys.exit(1)

    return int(clean.group(1)), int(aligned.group(1)), methy


def extract_dedu_report(report_file):
    with open(report_file, 'r') as f:
        content = f.read()

    dedu = re.search(r"Total count of deduplicated leftover sequences:\s*(\d+)", content)
    if not dedu:
        logging.error(f"{report_file} 缺失dedu配对相关字段，单端未实现")
        sys.exit(1)

    return int(dedu.group(1))


def extract_methy_extract_reports(log_file):
    with open(log_file, 'r') as f:
        content = f.read()

    methy = re.search(r"C methylated in CpG context:\s*([\d.]+)", content)
    if not methy:
        logging.error(f"{log_file} 缺失methy配对相关字段，单端未实现")
        sys.exit(1)

    methy = float(methy.group(1)) * 0.01  # report里显示的是百分比

    return methy


def methy_ratio_of_lambda(file_path):
    total_methy = 0
    total_count = 0

    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if columns[0] == 'chrL':
                try:
                    fifth_value = float(columns[4])
                    sixth_value = float(columns[5])
                    total_methy += fifth_value
                    total_count += (fifth_value + sixth_value)
                except ValueError:
                    continue

    return total_methy / total_count if total_count else 0


def main():
    raw_reads_dict = {}
    clean_reads_dict = {}
    mapped_reads_dict = {}
    mapping_rate_dict = {}
    methy_ratio_dict = {}
    dedu_reads_dict = {}
    dedu_ratio_dict = {}
    methy_ratio_dedu_dict = {}
    convertion_ratio = {}

    # 1. trim log
    for log_file in trim_logs_dir.glob("*.log"):
        sample = log_file.stem
        raw_reads_dict[sample] = extract_from_trim_log(log_file)

    # 2. bismark report
    for report_file in bismark_report_dir.glob("*_1_val_1_bismark_bt2_PE_report.txt"):
        sample = report_file.stem.replace("_1_val_1_bismark_bt2_PE_report", "")
        #print(sample)

        clean, aligned, methy = extract_bis_report(report_file)

        clean_reads_dict[sample] = clean
        mapped_reads_dict[sample] = aligned
        methy_ratio_dict[sample] = methy
        mapping_rate_dict[sample] = aligned / clean if clean else 1

    # 3. deduplication report
    for report_file in dedu_report_dir.glob("*_1_val_1_bismark_bt2_pe.deduplication_report.txt"):
        sample = report_file.stem.replace("_1_val_1_bismark_bt2_pe.deduplication_report", "")
        dedu_reads = extract_dedu_report(report_file)
        dedu_reads_dict[sample] = dedu_reads
        dedu_ratio_dict[sample] = dedu_reads / mapped_reads_dict.get(sample, 1)

    # 4. extract methylation logs
    for log_file in extract_methy_logs_dir.glob("*.log"):
        sample = log_file.stem
        methy_ratio_dedu_dict[sample] = extract_methy_extract_reports(log_file)

    # 5. lambda DNA 甲基化比率
    for bed_file in bed_dir.glob("*/sorted.bed"):
        sample = bed_file.parent.name
        print(sample)
        convertion_ratio[sample] = 1 - methy_ratio_of_lambda(bed_file)

    # 6. 构建 DataFrame
    df = pd.DataFrame({
        "raw_reads": raw_reads_dict,
        "clean_reads": clean_reads_dict,
        "mapped_reads": mapped_reads_dict,
        "mapping_rate": mapping_rate_dict,
        "methy_ratio": methy_ratio_dict,
        "dedu_reads": dedu_reads_dict,
        "dedu_ratio": dedu_ratio_dict,
        "methy_ratio_dedu": methy_ratio_dedu_dict,
        "convertion_ratio": convertion_ratio,
    }).reset_index().rename(columns={"index": "filename"})

    # 格式化
    df["raw_reads"] = df["raw_reads"].fillna(0).astype(int)
    df["clean_reads"] = df["clean_reads"].fillna(0).astype(int)
    df["mapped_reads"] = df["mapped_reads"].fillna(0).astype(int)
    df["dedu_reads"] = df["dedu_reads"].fillna(0).astype(int)
    for col in ["mapping_rate", "dedu_ratio", "methy_ratio", "methy_ratio_dedu", "convertion_ratio"]:
        if col in df.columns:
            df[col] = df[col].astype(float).round(3)

    # 输出结果
    df.to_csv(outputfile, sep='\t', index=False, encoding='utf-8')
    print("Wrote summary file to {}".format(outputfile))


if __name__ == "__main__":
    main()
