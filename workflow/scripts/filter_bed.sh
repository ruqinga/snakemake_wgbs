#!/bin/bash

# 参数
fai_file=$1
sorted_bed=$2
filtered_bed=$3
extract_threshold=$4


# 筛选total count
# 使用awk过滤，确保结束坐标不超出染色体的大小
awk -F'\t' -v OFS='\t' -v fai_file="$fai_file" -v extract_threshold="$extract_threshold" '
  # 读取.fai文件并存入数组
  NR==FNR {sizes[$1]=$2; next}

  # 如果条件满足，且结束坐标不超出染色体长度
  $5 + $6 > extract_threshold && $5 + $6 > 0 && $2 < sizes[$1] {
    print $1, $2, $2+1, ($5 / ($5 + $6))
  }
' "$fai_file" "$sorted_bed" > "$filtered_bed"