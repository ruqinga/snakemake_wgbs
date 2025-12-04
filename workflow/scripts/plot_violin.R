library(dplyr)
library(ggplot2)

##### load data #####
source("./scr/load_data.R")
# methy
data <- load_data(path = "./bed", pattern = "*_bin_100_t10.bed.slop_tmp",
                  col_name = c(level = "V4"),
                  output = "df")
print("load data finish")

##### tidy #####
df_rename <- data %>%
  mutate(
    cell_type = case_when(
      grepl("TKO", group) ~ "TKO",
      grepl("DKO", group) ~ "DKO",
      grepl("WT", group) ~ "WT",
      grepl("D1KO", group) ~ "D1KO",
      TRUE ~ NA_character_  # 如果没有匹配到任何内容，返回 NA
    ),
    source = case_when(
      grepl("J1", group) ~ "public",
      TRUE ~ "inhouse"  # 如果没有匹配到 J1，默认设置为 inhouse
    )
  )
write.table(df_rename, "df_rename.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##### plot #####

draw <- df_rename
col <- c("#BB0021FF","#3B4992FF")
point_width <- 1/(length(unique(draw$cell_type))^(4/3))

p <- ggplot(draw, aes(x = group, y = level, fill = source)) +
  geom_violin(position = "identity", linewidth = 0.5, width = 0.75)+
  geom_boxplot(width = 0.02, outlier.color = NA, fill = "black") +
  stat_summary(fun = median, geom = "crossbar", width = point_width, color = "white", linewidth = 0.6)+
      scale_fill_manual(values = col) +
  facet_grid(~ cell_type, scales = "free") +
  #scale_fill_manual(values = col) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 10, color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("plot_violin_v2.png", width = 10, height = 6, dpi = 600)

