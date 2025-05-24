#PBS -N extract_Dnmt_bai.pbs
#PBS -l nodes=1:ppn=20
#PBS -S /bin/bash
#PBS -j oe
#PBS -q pub_fat

source ~/.bashrc
conda activate base-omics

# 定义目标区域
Dnmt1="chr9:20907209-20959888"
Dnmt3A="chr12:3806007-3917443"
Dnmt3B="chr2:153649450-153687730"
Dnmt3C="chr2:153696652-153729907"

# 定义输入文件夹路径
input_dir="/home_data/home/slst/leixy2023/data/project/DNMT3C/BS/2025_rebuild_cell_line/250515_DKO/snakemake_wgbs_3/Results/05_dedu"

# 进入输入文件夹
cd "$input_dir"

# 遍历文件夹中的所有 .bam 文件
for bam_file in *.bam; do
    # 获取文件名前缀
    prefix=$(basename "$bam_file" .bam)
    echo $prefix

    # 对每个文件进行排序和索引
    samtools index -@ 20 "${prefix}.bam"

    # 提取所有区域并合并到一个文件中
    samtools view -b "${prefix}.bam" "$Dnmt1" "$Dnmt3A" "$Dnmt3B" "$Dnmt3C" > "${prefix}_Dnmt.bam"

    # 为合并后的文件生成索引
    samtools index "${prefix}_Dnmt.bam"
done

echo "处理完成！"