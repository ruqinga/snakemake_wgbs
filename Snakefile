configfile: "config.yaml"

include: "rules/common.smk"

# 获取全局参数
print(f"Config: {config}")

directories = get_directories(config)
samples = get_sample_list(config)

rule all:
    input:
        get_all(directories, samples)


#load rules
include: "rules/trim.smk"
#include: "rules/cut_bismark.smk" # 处理我们的数据
include: "rules/bismark.smk" # 处理公开数据
include: "rules/deduplicate.smk" # 去重，无法修改去重阈值
include: "rules/extract.smk" # 提取甲基化信息

