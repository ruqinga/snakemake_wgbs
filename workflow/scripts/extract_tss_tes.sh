#!/usr/bin/env bash
set -euo pipefail

# 读取参数
sorted_bed=$1
sorted_v2_bedGraph=$2
index_bed=$3
intersect_output=$4
bw_threshold=$5
relative_distance_table=$6
norm_mean_output=$7
mean_py=$8
log_file=$9

# 最后一列改为 total count
awk '{print $1,$2,$3,$4,$5+$6}' "$sorted_bed" > "$sorted_v2_bedGraph"
sed -i 's/ \+/\t/g' "$sorted_v2_bedGraph"

# intersect
bedtools intersect -a "$index_bed" -b "$sorted_v2_bedGraph" -wb > "$intersect_output"
echo "$(date) intersect finish, total lines: $(wc -l < "$intersect_output")" >> "$log_file"
head "$intersect_output" >> "$log_file"

# 计算相对于全长的相对长度
awk -v threshold="$bw_threshold" 'BEGIN { OFS = "\t"; }
{
    if (NF >= 13 && $13 >= threshold) {
        if ($8 == "+") {
            tss_edge = $5;
        } else if ($8 == "-") {
            tss_edge = $6;
        } else {
            next;
        }
        mehtylation_pos = $10;
        distance = tss_edge - mehtylation_pos;
        if (distance < 0) distance = -distance;
        normdis = distance / $4;
        print $7, $4, normdis, $12;
    }
}' "$intersect_output" > "$relative_distance_table"
echo "$(date) relative distance calculation complete, lines: $(wc -l < "$relative_distance_table")" >> "$log_file"

# 求平均
python "$mean_py" "$relative_distance_table" "$norm_mean_output" -g 1,3 -v 4
echo "$(date) mean calculation done" >> "$log_file"
