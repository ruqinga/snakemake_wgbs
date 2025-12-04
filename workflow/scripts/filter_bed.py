# -*- coding: utf-8 -*-
import argparse, logging, pandas as pd, subprocess, sys, os

def setup_logging(log_level=logging.INFO):
    logging.basicConfig(level=log_level, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

def parse_args():
    parser = argparse.ArgumentParser(description="Filter BED file and convert to BigWig.")
    parser.add_argument("--fai_file", type=str, default="/public/slst/home/qushy/toolkit/reference_genome/index/WGBS-index/mm10-lambdaDNA.fa.fai",
                        help="Path to .fai genome index file.") #使用mm10-lambdaDNA.fa.fai是为了计算lambdaDNA的转换率
    parser.add_argument("-i","--input", type=str, required=True, help="Input sorted BED file (chr, start, end, ratio, methy, unmethy).")
    parser.add_argument("--output_bed", type=str, help="Output filtered bedGraph file (default: basename.filtered.bed).")
    parser.add_argument("--output_bw", type=str, help="Output BigWig file (default: basename.filtered.bw).")
    parser.add_argument("-t","--extract_threshold", type=float, required=True, help="Threshold for filtering based on total count.")
    parser.add_argument("--bedGraphToBigWig", type=str, default="bedGraphToBigWig", help="Path to bedGraphToBigWig executable.")
    return parser.parse_args()

def load_fai(fai_file):
    chrom_sizes = {}
    with open(fai_file, "r") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom, size = parts[0], int(parts[1])
                chrom_sizes[chrom] = size
    logging.info(f"Loaded {len(chrom_sizes)} chromosomes from {fai_file}")
    return chrom_sizes

def filter_bed(sorted_bed, filtered_bed, chrom_sizes, extract_threshold):
    df = pd.read_csv(sorted_bed, sep='\t', header=None, usecols=[0, 1, 4, 5],
                     names=['chr', 'start', 'methy', 'unmethy'],
                     dtype={'chr': str, 'start': int, 'methy': float, 'unmethy': float})
    df['total'] = df['methy'] + df['unmethy']
    df['ratio'] = df['methy'] / df['total']
    df = df[(df['total'] > extract_threshold) & (df['total'] > 0)]
    df = df[df.apply(lambda row: row['start'] < chrom_sizes.get(row['chr'], 0), axis=1)]
    df['end'] = df['start'] + 1
    df_out = df[['chr', 'start', 'end', 'ratio']]
    df_out.to_csv(filtered_bed, sep='\t', header=False, index=False, float_format='%.6f')
    logging.info(f"Filtered BED written to {filtered_bed}, {len(df_out)} records.")

def run_bedGraphToBigWig(filtered_bed, fai_file, output_bw, bedGraphToBigWig="bedGraphToBigWig"):
    cmd = [bedGraphToBigWig, filtered_bed, fai_file, output_bw]
    logging.info(f"Running: {' '.join(cmd)}")
    try:
        subprocess.check_call(cmd)
        logging.info(f"BigWig file generated: {output_bw}")
    except subprocess.CalledProcessError as e:
        logging.error(f"bedGraphToBigWig failed: {e}")
        sys.exit(1)

def main():
    args = parse_args()
    setup_logging()

    # 自动生成默认输出文件名
    base = os.path.splitext(os.path.basename(args.input))[0]
    output_bed = args.output_bed if args.output_bed else f"{base}_t{args.extract_threshold}.bed"
    output_bw = args.output_bw if args.output_bw else f"{base}_t{args.extract_threshold}.bw"

    chrom_sizes = load_fai(args.fai_file)
    filter_bed(args.input, output_bed, chrom_sizes, args.extract_threshold)
    run_bedGraphToBigWig(output_bed, args.fai_file, output_bw, args.bedGraphToBigWig)

if __name__ == "__main__":
    main()
