#PBS -N visualization.pbs
#PBS -l nodes=1:ppn=10
#PBS -S /bin/bash
#PBS -j oe
#PBS -q slst_fat
#PBS -l walltime=60:00:00 


source ~/.bashrc
conda activate base-omics

work_dir="/home_data/home/slst/leixy2023/data/project/DNMT3C/BS/2025_rebuild_cell_line/251124_QKO-3COE-PGC-test/snakemake_wgbs"

out_dir="${work_dir}/Results/07_visualization_manual"

mkdir -p "${out_dir}"/bed "${out_dir}"/bw

bin_cutoff=10
base_cutoff=5

# 单独指定
# for sample in 3CTKO-7-rep1_raw 3CTKO-7-rep2_raw QKO-27-rep1_raw; do

# 批量处理
for dir in "${work_dir}"/Results/06_extract/*/; do
  sample=$(basename "$dir")

  echo "Processing sample: $sample"
  # sort bed
  #zcat ${work_dir}/Results/06_extract/${sample}/${sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | sort -k1,1 -k2,2n > ${work_dir}/Results/06_extract/${sample}/sorted.bed

  # 筛选total count
  #bash ${work_dir}/workflow/scripts/filter_bed.sh /public/slst/home/qushy/toolkit/reference_genome/index/WGBS-index/mm10-lambdaDNA.fa.fai ${work_dir}/Results/06_extract/${sample}/sorted.bed ${out_dir}/bed/${sample}_t${base_cutoff}_forbw.bed ${base_cutoff}
  # 转为bw
  #bedGraphToBigWig ${out_dir}/bed/${sample}_t${base_cutoff}_forbw.bed /public/slst/home/qushy/toolkit/reference_genome/index/WGBS-index/mm10-lambdaDNA.fa.fai ${out_dir}/bw/${sample}_t${base_cutoff}.bw

  # bin resolution
  python ${work_dir}/workflow/scripts/binning.py \
      -i ${work_dir}/Results/06_extract/${sample}/sorted.bed \
      -o ${out_dir} \
      --chrom_fai /public/slst/home/qushy/toolkit/reference_genome/index/WGBS-index/mm10-lambdaDNA.fa.fai \
      --bin_size 100 \
      --count_threshold ${bin_cutoff} \
      --sorted
done



