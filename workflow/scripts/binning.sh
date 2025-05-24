#PBS -N binning.pbs
#PBS -l nodes=1:ppn=10
#PBS -S /bin/bash
#PBS -j oe
#PBS -q slst_pub

source ~/.bashrc
conda activate base-omics

chrom_fai="/public/slst/home/leixy2023/database/mm10/mm10.fa.fai"

# define param
bin_size=100
total_count_thread=10
sorted=TRUE

input_dir="/home_data/home/slst/leixy2023/data/project/DNMT3C/BS/2025_rebuild_cell_line/250515_DKO/snakemake_wgbs_3/Results/07_visualization"
output_dir="/home_data/home/slst/leixy2023/data/project/DNMT3C/BS/2025_rebuild_cell_line/250515_DKO/snakemake_wgbs_3/Results/07_visualization/bed"

cd $output_dir

# produce each .bed
find $input_dir -name '*.bed' | while read -r file; do
    echo -e "Processing file: $file\n"
    base=$(basename "${file}" .bed)

    # sort file if not already sorted
    if [ "$sorted" = FALSE ]; then
        tmp_sort_file="${base}_sorted_tmp.bed"
        sort -k1,1 -k2,2n "$file" > "$tmp_sort_file"
        file="$tmp_sort_file"
    fi

    # binning
    tmp_file="${base}_binning_tmp.bed"
    awk -v bin_size=$bin_size '
    BEGIN { OFS="\t" }
    {
        chromosome = $1;
        bin_group = int($2 / bin_size);
        key = chromosome ":" bin_group;
        total_count[key] += $5;
        ratio_sum[key] += $4;
        count[key]++;
    }
    END {
        for (key in ratio_sum) {
            split(key, parts, ":");
            chromosome = parts[1];
            bin_group = parts[2];
            start_pos = bin_group * bin_size;
            end_pos = start_pos + bin_size;
            ratio_mean = ratio_sum[key] / count[key];
            print chromosome, start_pos, end_pos, ratio_mean, total_count[key];
        }
    }' "$file" > "$tmp_file"

    # sort the binned output
    sorted_output_file="${base}_bin_${bin_size}.bed"
    sort -k1,1 -k2,2n "$tmp_file" > "$sorted_output_file"
    echo -e "Sorted output saved to: $sorted_output_file\n"

    # filter and convert to BigWig format
    filtered_bed="${base}_bin_${bin_size}_t_${total_count_thread}_tmp.bed"
    awk -v total_count_thread=$total_count_thread '$5 > total_count_thread {print $1, $2, $3, $4}' "$sorted_output_file" > "$filtered_bed"
    # 截断超过部分
    sed -i 's/ \+/\t/g' "$filtered_bed"
    bedtools slop -i "$filtered_bed" -g "${chrom_fai}" -b 0 > "tmp_${filtered_bed}"

    bedGraphToBigWig "tmp_$filtered_bed" "${chrom_fai}" "${base}_bin_${bin_size}_t_${total_count_thread}.bw"
    #rm "$filtered_bed" "tmp_${filtered_bed}"

    echo -e "bw saved to: ${base}_bin_${bin_size}_t_${total_count_thread}.bw\n"
done
