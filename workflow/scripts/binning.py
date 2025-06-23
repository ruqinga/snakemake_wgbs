import argparse
import os
import subprocess
from pathlib import Path

def sort_bed_file(input_file, output_file):
    subprocess.run(f"sort -k1,1 -k2,2n {input_file} > {output_file}", shell=True, check=True)

def bin_bed_file(input_file, output_file, bin_size):
    awk_script = f"""
    awk -v bin_size={bin_size} '
    BEGIN {{ OFS="\\t" }}
    {{
        chromosome = $1;
        bin_group = int($2 / bin_size);
        key = chromosome ":" bin_group;
        total_meth[key] += $5;
        total_unmeth[key] += $6;
        count[key]++;
    }}
    END {{
        for (key in count) {{
            if ((total_meth[key] + total_unmeth[key]) == 0) continue;
            split(key, parts, ":");
            chromosome = parts[1];
            bin_group = parts[2];
            start_pos = bin_group * bin_size;
            end_pos = start_pos + bin_size;
            total_count = total_meth[key] + total_unmeth[key]
            ratio_mean = total_meth[key] / (total_count);
            print chromosome, start_pos, end_pos, ratio_mean, total_count;
        }}
    }}' {input_file} > {output_file}
    """
    subprocess.run(awk_script, shell=True, check=True)

def filter_binned_file(input_file, output_file, count_threshold):
    awk_cmd = f"awk -v total_count_thread={count_threshold} '$5 > total_count_thread {{print $1, $2, $3, $4}}' {input_file} > {output_file}"
    subprocess.run(awk_cmd, shell=True, check=True)

def bed_to_bigwig(filtered_bed, chrom_fai, output_bw):
    tmp_bed = f"{filtered_bed}.tmp"
    subprocess.run(f"sed -i 's/ \\+/\\t/g' {filtered_bed}", shell=True, check=True)
    subprocess.run(f"bedtools slop -i {filtered_bed} -g {chrom_fai} -b 0 > {tmp_bed}", shell=True, check=True)
    subprocess.run(f"bedGraphToBigWig {tmp_bed} {chrom_fai} {output_bw}", shell=True, check=True)
    # Optional cleanup:
    # os.remove(filtered_bed)
    os.remove(tmp_bed)

def process_file(file_path, output_dir, chrom_fai, bin_size, count_threshold, sorted_input):
    base = Path(file_path).parent.name
    print(f"Processing {base}")
    file_to_process = file_path

    if not sorted_input:
        sorted_tmp = os.path.join(output_dir, f"bed/{base}_sorted_tmp.bed")
        sort_bed_file(file_path, sorted_tmp)
        file_to_process = sorted_tmp

    binned_tmp = os.path.join(output_dir, f"bed/{base}_binning_{bin_size}_tmp.bed")
    bin_bed_file(file_to_process, binned_tmp, bin_size)

    sorted_binned = os.path.join(output_dir, f"bed/{base}_bin_{bin_size}.bed")
    sort_bed_file(binned_tmp, sorted_binned)

    filtered_bed = os.path.join(output_dir, f"bed/{base}_bin_{bin_size}_t_{count_threshold}_tmp.bed")
    filter_binned_file(sorted_binned, filtered_bed, count_threshold)

    output_bw = os.path.join(output_dir, f"bw/{base}_bin_{bin_size}_t_{count_threshold}.bw")
    bed_to_bigwig(filtered_bed, chrom_fai, output_bw)

    print(f"✓ Finished processing: {file_path}")
    print(f"→ Output BigWig: {output_bw}\n")

def main():
    parser = argparse.ArgumentParser(description="BED binning and BigWig conversion pipeline.")
    parser.add_argument('-i','--input_dir', required=True, help='Input directory with BED files.')
    parser.add_argument('-o','--output_dir', required=True, help='Directory to save processed files.')
    parser.add_argument('--chrom_fai', required=True, help='Path to .fa.fai file for genome.')
    parser.add_argument('--bin_size', type=int, default=100, help='Size of bins.')
    parser.add_argument('--count_threshold', type=int, default=10, help='Minimum total count to keep bin.')
    parser.add_argument('--sorted', action='store_true', help='Set if BED files are already sorted.')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    bed_files = list(Path(args.input_dir).rglob("*.bed"))
    if not bed_files:
        print(f"No .bed files found in {args.input_dir}")
        return

    for file_path in bed_files:
        process_file(
            str(file_path),
            args.output_dir,
            args.chrom_fai,
            args.bin_size,
            args.count_threshold,
            args.sorted
        )

if __name__ == '__main__':
    main()
