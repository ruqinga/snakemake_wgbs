#!/bin/bash
# 用于处理单端和双端数据并调用 Snakemake 运行流程
# 用法: bash run.sh <fq_dir> [--cut] [--pbat] [-y]
# nohup用法：nohup bash run.sh ../rawdata/WT_ln [--cut] [--pbat] -y > run.log 2>&1 &

# 默认值
Snakefile="./workflow/Snakefile"  # 默认选择 Snakefile
bismark_strategy=""

# 解析命令行选项
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --cut)
            Snakefile="./workflow/Snakefile_cut"  # 如果指定了 --cut，则选择 Snakefile
            shift
            ;;
        --pbat)
            bismark_strategy="--pbat"  # 如果指定了 --pbat，则设置 bismark_strategy
            shift
            ;;
        -y)
            confirm_run="y"  # -y 标志，表示确认
            shift
            ;;
        *)
            fq_dir="$1"  # 处理第一个参数 fq_dir
            shift
            ;;
    esac
done

# 确保 fq_dir 参数存在
if [ -z "$fq_dir" ]; then
    echo "错误: 未指定输入目录!"
    exit 1
fi

echo "输入数据目录: $fq_dir"
# 判断是否使用cut
if [ "$Snakefile" == "./workflow/Snakefile" ]; then
    echo "不使用 cut"
else
    echo "运行 cut_bismark 处理 zhenglab data"
fi

# 判断 bismark_strategy 的值来输出不同信息
if [ -z "$bismark_strategy" ]; then
    echo "处理 wgbs data"
else
    echo "处理 pbat data"
fi
echo "---------------"


# 激活 Snakemake 的 Conda 环境
source ~/.bashrc
conda activate snakemake_env

# 调用脚本生成单端和双端数据列表

# 初始化数组
json_array_pe=()
json_array_se=()
json_array=()

# 遍历所有 .fastq.gz 文件
for file in $(find "$fq_dir" -maxdepth 1 -name "*.fastq.gz" | sort); do

    file=$(realpath "$file")  # 处理路径中的特殊字符，比如空格

    # 检查是否是带 _1.fastq.gz 或 _2.fastq.gz 的 PE 文件
    if [[ "$file" =~ _1.fastq.gz$ ]]; then
        # 获取对应的 _2 文件
        file2="${file/_1.fastq.gz/_2.fastq.gz}"

        if [ -f "$file2" ]; then
            # 如果 _2 文件存在，则添加到 PE 数组
            json_entry_pe="{\"read1\": \"$file\", \"read2\": \"$file2\"}"
            json_array_pe+=("$json_entry_pe")
            json_array+=("$json_entry_pe")
        fi
    elif [[ ! "$file" =~ _1.fastq.gz$ && ! "$file" =~ _2.fastq.gz$ ]]; then
        # 如果不是 _1 或 _2，视为 SE 文件
        json_entry_se="{\"read1\": \"$file\"}"
        json_array_se+=("$json_entry_se")
        json_array+=("$json_entry_se")
    fi
done

# 将数组转换为 JSON 字符串
json_output_pe=$(IFS=,; echo "[${json_array_pe[*]}]")
json_output_se=$(IFS=,; echo "[${json_array_se[*]}]")
json_output=$(IFS=,; echo "[${json_array[*]}]")

# 保存 json_output 到文件
json_output_file="$fq_dir/reads_json.json"
echo "$json_output" > "$json_output_file"
echo "save input reads information to $json_output_file"

# 检查 JSON 数据是否为空并输出相应的提示信息
if [[ -z "$json_output_se" || "$json_output_se" == "[]" ]]; then
    echo "没有单端数据"
    json_output_se=""
else
    echo "即将处理如下单端数据:"
    echo "$json_output_se" | jq .
fi

if [[ -z "$json_output_pe" || "$json_output_pe" == "[]" ]]; then
    echo "没有双端数据"
    json_output_pe=""
else
    echo "即将处理如下双端数据:"
    echo "$json_output_pe" | jq .
fi


# 运行 Snakemake 工作流（预览模式）
echo "运行 Snakemake （仅预览）..."
snakemake \
    -np \
    --use-conda \
    --snakefile "$Snakefile" \
    --config fq_dir="$fq_dir" reads="$json_output"

# 提示是否确认实际执行任务
# 检查是否传递了 -y 参数
if [[ -z "$confirm_run" ]]; then
    read -p "是否确认执行任务（实际提交作业）？(y/n): " confirm_run
fi

# 如果 confirm_run 不是 "y" 或 "Y"，则取消任务
if [[ "$confirm_run" != "y" && "$confirm_run" != "Y" ]]; then
    echo "任务已取消！"
    exit 0
fi

# 实际运行 Snakemake 工作流
echo "运行 Snakemake ..."
snakemake \
    --snakefile "$Snakefile" \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "python workflow/scripts/submit_job.py --config config/cluster_config.yaml --seqtype "bs" --sample {wildcards} --rule {rule}" \
    --latency-wait 60 \
    --jobs 5 \
    --use-conda \
    --groups processing_group=20 Additional_analysis=10 \
    --config fq_dir="$fq_dir" reads="$json_output" bismark_strategy="$bismark_strategy"

bismark_pe
echo "任务已完成！"
