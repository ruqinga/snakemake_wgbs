library(dplyr)
library(ggplot2)
library(argparse)

# 输入参数解析
parser <- ArgumentParser(description = "Process input folder and output file name.")
parser$add_argument("--input_folder", type = "character", default = "./", help = "Input folder path")
parser$add_argument("--output_file", type = "character", default = "output.pdf", help = "Output file name")
parser$add_argument("--filter_group", type = "character", default = "L1", help = "Filter group (L1 or LTR)")

args <- parser$parse_args()

##### load data #####
source("./scr/load_data.R")
# methy
data <- load_data(path = args$input_folder, pattern = "*_t5_forbw_level_dis_mean.txt",
                  col_name = c(class = "V1", norm_dis = "V2", level = "V3"),
                  output = "df")
print("load data finish")
##### tidy #####
source("./scr/binning.R")
# 筛选数据
LTR <- c( "IAPEz-int","IAPEy-int","ETnERV3-int", #ERVK  "IAP-d-int",
             "MURVY-int", "MMERGLN-int", "MMERGLN_LTR", #ERV1
             "ORR1A1-int","ORR1A0","MERVL-int" # ERVL "ORR1A0-int",
             )
L1 <- c("L1Md_A","L1Md_T","L1Md_F2","L1Md_F3","L1Md_Gf","Lx")

# 如果class在LTR里，则增加rep_group列为LTR，如果在L1中，则增加rep_group列为L1
data_target <- data %>%
  filter(class %in% c(LTR, L1)) %>%
  mutate(rep_group = ifelse(class %in% LTR,"LTR", "L1"))

# 计算bin
df_bin <- my_bin(data_target, "level", c("group", "bin","class","rep_group"))
write.table(df_bin, "df_bin.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##### plot #####
col <- c( "#7E6148FF", "#168aad","#0e6ba8","#B09C85FF", "#3C5488FF", "#ff6f00")
glevel <- c("WT-R1","DKO-1","DKO-29", "J1_WT","J1_DKO", "J1_TKO")
draw <- df_bin %>%
  filter(group %in% glevel) %>%
  mutate(class = factor(class, levels = c(L1, LTR)), group = factor(group, levels = glevel)) %>%
  filter(rep_group == args$filter_group)

library(ggplot2)
my_theme <- theme(
      ### 背景 ###
      panel.background = element_blank(),  # 去掉图形背景
      plot.background = element_blank(),    # 去掉整个图的背景
      panel.grid = element_blank(),  # Remove grid lines

      ### 坐标轴 ###
      axis.line = element_line(color = 'black'),    # 保留坐标轴
      axis.ticks = element_line(color = 'black'),   # 保留坐标轴刻度
      axis.ticks.length = unit(0.05, "cm"),         # 设置坐标轴刻度长度

      ### 文字 ###
      text = element_text(size = 16),
      axis.text = element_text(color = 'black', size = 12, face = 'bold'),
      # axis.text.x = element_text(size = 12,face = 'bold'),  # angle = 45, hjust = 1 旋转x轴标签，使其可读性更强
      # axis.text.y = element_text(size = 12), # 修改y轴刻度标签的字体大小
      plot.title = element_text(size = 14), # 修改标题的字体大小
      # axis.title.x = element_text(size = 14),       # 修改x轴标签的字体大小
      # axis.title.y = element_text(size = 14),       # 修改y轴标签的字体大小
      legend.title = element_text(size = 14),   # 修改图例标题的字体大小
      legend.text = element_text( size = 14)# 修改图例文本的字体大小
      )

p1 <- ggplot(draw, aes(x = bin, y = level*100, color = group)) +
      geom_point(shape = 15, size = 2) +  # 添加散点图层
      geom_line(linewidth = 1)+
      # geom_smooth(method = "loess", se = FALSE) +  # 添加平滑线，使用loess方法
      scale_x_continuous(
        breaks = c(0, 0.5, 1.0, 1.5,2.0),  # 设置横轴的刻度
        labels = c("", "5'UTR", "","3'UTR", "")  # 自定义刻度标签（可选）
      ) +
      scale_y_continuous(limits = c(0, 100)) +
      facet_wrap(~ class, ncol = 3, scales = "free_y") +
      labs(x = " ", y = "CpG methylation (%)",title = " ")+
      scale_color_manual(values = col) + my_theme

ggsave(args$output_file, width = 10, height = 6,dpi = 600)