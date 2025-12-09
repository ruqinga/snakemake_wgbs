#!/usr/bin/env python3
"""
Minimal patched version of the original data_summary_extractor.py.

This patch fixes:
- incorrect use of '|' between glob strings
- broken .replace() with malformed string literals
- brittle regex handling for methylation percentage
- avoids sys.exit inside helper functions (return None on missing fields)
- consistent output filename/separator (.tsv with '\t')
- safer sample name extraction using regex
- minor logging improvements

This is a minimal patch — structure and overall logic kept as in your original script.
"""
from __future__ import annotations

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


def extract_from_trim_log(log_file: Path) -> int:
    with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    input_reads = re.search(r"Total number of sequences analysed:\s*([\d,]+)", content)
    if not input_reads:
        input_reads = re.search(r"(\d+)\s+sequences processed in total", content, re.IGNORECASE)

    if not input_reads:
        logging.warning(f"{log_file} 没有找到 input_reads")
        return 0

    return int(input_reads.group(1).replace(",", ""))


def extract_bis_report(report_file: Path):
    """
    返回 (clean_reads:int or None, aligned:int or None, methy:float in 0-1 or None).
    不在此处 sys.exit，调用者根据返回值决定如何处理。
    """
    with open(report_file, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    # Clean reads: try PE first then SE patterns
    clean = re.search(r"Sequence pairs analysed in total:\s*([\d,]+)", content)
    if not clean:
        clean = re.search(r"Sequences analysed in total:\s*([\d,]+)", content)

    # Aligned reads: try paired-end then single-end patterns
    aligned = re.search(
        r"Number of paired-end alignments with a unique best hit:\s*([\d,]+)",
        content,
    )
    if not aligned:
        aligned = re.search(
            r"Number of alignments with a unique best hit from the different alignments:\s*([\d,]+)", content
        )

    # Methylation: prefer percentage forms
    methy = re.search(r"C methylated in CpG context:\s*([\d.]+)%", content)
    if not methy:
        # sometimes the report uses "Total methylated C's in CpG context: 123 (3.45%)"
        methy = re.search(
            r"Total methylated C's in CpG context:.*\(([\d.]+)%\)", content
        )
    if not methy:
        # fallback: numeric without percent
        methy = re.search(r"C methylated in CpG context:\s*([\d.]+)", content)

    # parse values
    clean_val = int(clean.group(1).replace(",", "")) if clean else None
    aligned_val = int(aligned.group(1).replace(",", "")) if aligned else None
    methy_val = None
    if methy:
        s = methy.group(1)
        try:
            f = float(s)
            # decide if it's percentage: prefer percent-form matched earlier,
            # but if fallback matched no-percent and value > 1 assume percent
            if "%" in methy.group(0) or f > 1.0:
                methy_val = f * 0.01
            else:
                methy_val = f
        except ValueError:
            methy_val = None

    # If critical fields are missing, warn and return Nones
    if clean_val is None or aligned_val is None:
        logging.warning(f"{report_file} clean{clean_val},aligned{aligned_val},methy{methy_val}")
        logging.warning(f"{report_file} 缺失bis配对相关字段（clean/aligned），跳过该文件")
        return None, None, None

    return clean_val, aligned_val, methy_val


def extract_dedu_report(report_file: Path) -> int:
    with open(report_file, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    dedu = re.search(
        r"Total count of deduplicated leftover sequences:\s*([\d,]+)", content, re.IGNORECASE
    )

    if not dedu:
        logging.warning(f"{report_file} 没有找到 dedu")
        return 0

    return int(dedu.group(1).replace(",", ""))


def extract_methy_extract_reports(log_file: Path):
    with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    methy = re.search(r"C methylated in CpG context:\s*([\d.]+)%", content)
    if not methy:
        methy = re.search(r"C methylated in CpG context:\s*([\d.]+)", content)

    if not methy:
        logging.warning(f"{log_file} 缺失methy配对相关字段，返回 None")
        return None

    try:
        val = float(methy.group(1))
    except ValueError:
        logging.warning(f"{log_file} methy 无法转换为数值，返回 None")
        return None

    # 如果匹配到带%或数值>1则视为百分比
    if "%" in methy.group(0) or val > 1.0:
        return val * 0.01
    return val


def methy_ratio_of_lambda(file_path: Path):
    total_methy = 0.0
    total_count = 0.0

    with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            columns = line.strip().split()
            if not columns:
                continue
            if columns[0] == "chrL":
                if len(columns) < 6:
                    continue
                try:
                    fifth_value = float(columns[4].rstrip("%"))
                    sixth_value = float(columns[5].rstrip("%"))
                    total_methy += fifth_value
                    total_count += (fifth_value + sixth_value)
                except ValueError:
                    continue

    return total_methy / total_count if total_count else None


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
    # iterate over the two patterns instead of using '|' on strings
    patterns = [
        "*_1_val_1_bismark_bt2_PE_report.txt",
        "*_trimmed_bismark_bt2_SE_report.txt",
        "*_bismark_bt2_PE_report.txt",
        "*_bismark_bt2_SE_report.txt",
    ]
    for pat in patterns:
        for report_file in bismark_report_dir.glob(pat):
            # derive sample name by removing known suffixes using regex
            fname = report_file.name
            sample = re.sub(
                r"(_1_val_1_bismark_bt2_PE_report|_trimmed_bismark_bt2_SE_report|_bismark_bt2_PE_report|_bismark_bt2_SE_report)\.txt$",
                "",
                fname,
                flags=re.IGNORECASE,
            )
            # parse report
            clean, aligned, methy = extract_bis_report(report_file)

            if clean is None or aligned is None:
                # already logged inside extract_bis_report
                continue

            clean_reads_dict[sample] = clean
            mapped_reads_dict[sample] = aligned
            methy_ratio_dict[sample] = methy if methy is not None else 0.0
            mapping_rate_dict[sample] = aligned / clean if clean else 0.0

    # 3. deduplication report
    dedu_patterns = [
        "*_1_val_1_bismark_bt2_pe.deduplication_report.txt",
        "*_trimmed_bismark_bt2.deduplication_report.txt",
        "*_bismark_bt2_pe.deduplication_report.txt",
        "*_bismark_bt2.deduplication_report.txt",
    ]
    for pat in dedu_patterns:
        for report_file in dedu_report_dir.glob(pat):
            fname = report_file.name
            sample = re.sub(
                r"(_1_val_1_bismark_bt2_pe\.deduplication_report|_trimmed_bismark_bt2\.deduplication_report|_bismark_bt2_pe\.deduplication_report|_bismark_bt2\.deduplication_report)\.txt?$",
                "",
                fname,
                flags=re.IGNORECASE,
            )
            dedu_reads = extract_dedu_report(report_file)
            dedu_reads_dict[sample] = dedu_reads
            mapped = mapped_reads_dict.get(sample)
            try:
                dedu_ratio_dict[sample] = dedu_reads / mapped if (mapped and mapped > 0) else 0.0
            except Exception:
                dedu_ratio_dict[sample] = 0.0

    # 4. extract methylation logs
    for log_file in extract_methy_logs_dir.glob("*.log"):
        sample = log_file.stem
        methy_ratio_dedu_dict[sample] = extract_methy_extract_reports(log_file)

    # 5. lambda DNA 甲基化比率
    for bed_file in bed_dir.glob("*/sorted.bed"):
        sample = bed_file.parent.name
        val = methy_ratio_of_lambda(bed_file)
        convertion_ratio[sample] = 1 - val if val is not None else 0.0

    # 6. 构建 DataFrame
    df = pd.DataFrame(
        {
            "raw_reads": raw_reads_dict,
            "clean_reads": clean_reads_dict,
            "mapped_reads": mapped_reads_dict,
            "mapping_rate": mapping_rate_dict,
            "methy_ratio": methy_ratio_dict,
            "dedu_reads": dedu_reads_dict,
            "dedu_ratio": dedu_ratio_dict,
            "methy_ratio_dedu": methy_ratio_dedu_dict,
            "convertion_ratio": convertion_ratio,
        }
    ).reset_index().rename(columns={"index": "filename"})

    # 格式化
    df["raw_reads"] = df["raw_reads"].fillna(0).astype(int)
    df["clean_reads"] = df["clean_reads"].fillna(0).astype(int)
    df["mapped_reads"] = df["mapped_reads"].fillna(0).astype(int)
    df["dedu_reads"] = df["dedu_reads"].fillna(0).astype(int)
    for col in ["mapping_rate", "dedu_ratio", "methy_ratio", "methy_ratio_dedu", "convertion_ratio"]:
        if col in df.columns:
            df[col] = df[col].astype(float).round(4)

    # 输出结果
    outputfile.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outputfile, sep="\t", index=False, encoding="utf-8")
    print("Wrote summary file to {}".format(outputfile))


if __name__ == "__main__":
    main()