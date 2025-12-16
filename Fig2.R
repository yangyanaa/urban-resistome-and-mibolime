#Fig2 a-b
##arg HGT 转移的event
rm(list=ls())
folder <- "/Users/yangy/Documents/E/city/metagenome/HGT/diamond"

# 找到所有 .tsv 文件（全路径）
files <- list.files(folder, pattern = "\\.f6$", full.names = TRUE)

# 结果表
result <- data.frame(File = character(), Lines_without_header = integer(), stringsAsFactors = FALSE)

for (f in files) {
  # 读取所有行数（如果文件很大会占内存，这种方法在多数常见文件下是可行的）
  total_lines <- length(readLines(f))
  
  # 去掉 _SARG_diamond.f6 后缀
  sample_name <- sub("_SARG_diamond\\.f6$", "", basename(f))
  
  # 保存结果
  result <- rbind(result, data.frame(File = sample_name,
                                     Lines_without_header = total_lines,
                                     stringsAsFactors = FALSE))
}




colnames(result)[1]<-"sample.id"
colnames(result)[2]<-"count"


#write.csv(result, "/Users/yangy/Documents/E/city/metagenome/HGT/results/hgt_arg_gene_count.csv")


##计算HGT的数量
library(agricolae)
library(ggplot2)
result<-read.csv("/Users/yangy/Documents/E/city/metagenome/HGT/results/hgt_arg_gene_count.csv",sep=",",header=T)
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group <- group0[order(group0$sample.id), ]
result <- result[order(result$sample.id), ]

result$group <- group$group[match(result$sample.id, group$sample.id)]
group<-result

kruskal.test(group$count,group$group)
library(broom)
tidy(kruskal.test(group$count,group$group))

ST <- kruskal(group$count,group$group, p.adj = "fdr")
tidy(kruskal.test(group$count,group$group))

statis <- as.data.frame(rbind(ST$statistics[,1:2]))
rp.value <- data.frame(
  label=sprintf("P = %.3f (χ² = %.3g)", 
                statis$p.chisq,statis$Chisq),
  x = "Tamarix chinensis",
  y=0.55)
statis1 <- as.data.frame(rbind(ST$groups))
statis2 <- as.data.frame(rbind(ST$means))
p1 <- ggplot(group, aes(y=count, x = factor(group))) +
  stat_boxplot(geom = "errorbar", width=0.3) +
  geom_boxplot(width = 0.5)
p1
ymaxdata <- layer_data(p1) 
statis2$ymax <- ymaxdata$ymax
statis1$id <- rownames(statis1)
statis2$id <- rownames(statis2)
statis2 <- merge(statis2, statis1, by = "id")

multcomp.data <- data.frame(
  abc = factor(statis2$groups),
  group     = statis2$id,
  y     = as.numeric(statis2$ymax))

library(gghalves)


chi2_value=14.3
p_value=0.00249
format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 科学计数法
  } else {
    round(p, 3)  # 保留三位小数
  }
}
p_to_plotmath <- function(p) {
  # 科学计数法，如 "1.86e-03"
  sci <- format(p, scientific = TRUE)
  
  # 转成 plotmath，例如 "1.86%*%10^-3"
  sci <- gsub("e-", "%*%10^-", sci)
  sci <- gsub("e\\+", "%*%10^", sci)
  
  return(sci)
}
p <- ggplot(group, aes(x = group, y = count, fill = group, color = group)) +
  # 设置颜色
  scale_fill_manual(values = c(
    "Forest" = "#6cb509",
    "Cropland" = "cyan",
    "Industry" = "#5e9ceb",
    "Urban" = "#ff7f27"
  )) +
  scale_color_manual(values = c(
    "Forest" = "#6cb509",
    "Cropland" = "cyan",
    "Industry" = "#5e9ceb",
    "Urban" = "#ff7f27"
  )) +
  scale_x_discrete(limits = c("Forest", "Cropland", "Industry", "Urban")) +
  
  # 小提琴图
  geom_violin(trim = FALSE, alpha = 0.6, color = NA, width = 0.8) +
  
  # 箱线图（白色箱体）
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.8) +
  
  # 散点
  geom_jitter(size = 2.5, alpha = 0.6, width = 0.15, aes(color = group)) +
  
  # 显著性字母
  geom_text(
    data = multcomp.data,
    aes(x = group, y = y + 0.01, label = abc),
    color = "black", size = 6, hjust = 0.5
  ) +
  
  # 自动格式化的卡方与P值标注
  annotate(
    "text", x = 2.5, y = 25,
    label = paste0(
      "chi^2==", chi2_value,
      "~~P==", p_to_plotmath(p_value)
    ),
    parse = TRUE,
    size = 5, color = "black"
  )+
  
  # y轴标签
  ylab("HGT events(ARG)") +
  
  # 主题
  theme_light() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title.y = element_text(size = 16, color = "black", face = "bold",
                                vjust = 1.9, hjust = 0.5, angle = 90),
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "grey80", size = 0.1),
    panel.grid.minor = element_line(color = "grey90", size = 0.1),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 15, colour = "black", face = "bold", hjust = 0.5)
  )

p


#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/hgt_counts.pdf",p,height=5,width=6,device = cairo_pdf)



##计算arg转移基因的richness
rm(list=ls())
folder <- "E:/city/metagenome/HGT/diamond"

# 找到所有 .f6 文件
files <- list.files(folder, pattern = "\\.f6$", full.names = TRUE)

# 结果表
result <- data.frame(sample.id = character(), richness = integer(), stringsAsFactors = FALSE)

for (f in files) {
  # 检查文件是否为空
  if (file.info(f)$size == 0) {
    gene_richness <- 0
  } else {
    df <- tryCatch(
      read.table(f, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(df) || nrow(df) == 0) {
      gene_richness <- 0
    } else {
      gene_richness <- length(unique(df[[2]]))
    }
  }
  
  sample_name <- sub("_SARG_diamond\\.f6$", "", basename(f))
  
  result <- rbind(result, data.frame(sample.id = sample_name,
                                     richness = gene_richness,
                                     stringsAsFactors = FALSE))
}

result
colnames(result)[1]<-"sample.id"
colnames(result)[2]<-"richness"


#write.csv(result, "E:/city/metagenome/HGT/results/hgt_arg_gene_richness.csv")

result<-read.csv("/Users/yangy/Documents/E/city/metagenome/HGT/results/hgt_arg_gene_richness.csv",sep=",",header=T)
library(agricolae)
library(ggplot2)
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group <- group0[order(group0$sample.id), ]
result <- result[order(result$sample.id), ]

result$group <- group$group[match(result$sample.id, group$sample.id)]
group<-result

kruskal.test(group$richness,group$group)
library(broom)
tidy(kruskal.test(group$richness,group$group))

ST <- kruskal(group$richness,group$group, p.adj = "fdr")
tidy(kruskal.test(group$richness,group$group))

statis <- as.data.frame(rbind(ST$statistics[,1:2]))
rp.value <- data.frame(
  label=sprintf("P = %.3f (χ² = %.3g)", 
                statis$p.chisq,statis$Chisq),
  x = "Tamarix chinensis",
  y=0.55)
statis1 <- as.data.frame(rbind(ST$groups))
statis2 <- as.data.frame(rbind(ST$means))
p1 <- ggplot(group, aes(y=richness, x = factor(group))) +
  stat_boxplot(geom = "errorbar", width=0.3) +
  geom_boxplot(width = 0.5)
p1
ymaxdata <- layer_data(p1) 
statis2$ymax <- ymaxdata$ymax
statis1$id <- rownames(statis1)
statis2$id <- rownames(statis2)
statis2 <- merge(statis2, statis1, by = "id")

multcomp.data <- data.frame(
  abc = factor(statis2$groups),
  group     = statis2$id,
  y     = as.numeric(statis2$ymax))

library(gghalves)


chi2_value=15
p_value=0.00186
format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 科学计数法
  } else {
    round(p, 3)  # 保留三位小数
  }
}

p_to_plotmath <- function(p) {
  # 科学计数法，如 "1.86e-03"
  sci <- format(p, scientific = TRUE)
  
  # 转成 plotmath，例如 "1.86%*%10^-3"
  sci <- gsub("e-", "%*%10^-", sci)
  sci <- gsub("e\\+", "%*%10^", sci)
  
  return(sci)
}
annotate(
  "text", x = 2.5, y = 25,
  label = paste0(
    "chi^2==", chi2_value,
    "~~P==", p_to_plotmath(p_value)
  ),
  parse = TRUE,
  size = 5, color = "black"
)

p <- ggplot(group, aes(x = group, y = richness, fill = group, color = group)) +
  # 设置颜色
  scale_fill_manual(values = c(
    "Forest" = "#6cb509",
    "Cropland" = "cyan",
    "Industry" = "#5e9ceb",
    "Urban" = "#ff7f27"
  )) +
  scale_color_manual(values = c(
    "Forest" = "#6cb509",
    "Cropland" = "cyan",
    "Industry" = "#5e9ceb",
    "Urban" = "#ff7f27"
  )) +
  scale_x_discrete(limits = c("Forest", "Cropland", "Industry", "Urban")) +
  
  # 小提琴图
  geom_violin(trim = FALSE, alpha = 0.6, color = NA, width = 0.8) +
  
  # 箱线图（白色箱体）
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.8) +
  
  # 散点
  geom_jitter(size = 2.5, alpha = 0.6, width = 0.15, aes(color = group)) +
  
  # 显著性字母
  geom_text(
    data = multcomp.data,
    aes(x = group, y = y + 0.01, label = abc),
    color = "black", size = 6, hjust = 0.5
  ) +
  
  # 自动格式化的卡方与P值标注
  annotate(
    "text", x = 2.5, y = 25,
    label = paste0(
      "chi^2==", chi2_value,
      "~~P==", p_to_plotmath(p_value)
    ),
    parse = TRUE,
    size = 5, color = "black"
  )+
  
  # y轴标签
  ylab("The richness of\n horizontal gene transferred ARGs") +
  
  # 主题
  theme_light() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title.y = element_text(size = 16, color = "black", face = "bold",
                                vjust = 1.9, hjust = 0.5, angle = 90),
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "grey80", size = 0.1),
    panel.grid.minor = element_line(color = "grey90", size = 0.1),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 15, colour = "black", face = "bold", hjust = 0.5)
  )
p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/hgt_richnesss.pdf",p,height = 5,width = 6,device = cairo_pdf)

#c-f弦图
#HGT taxa
rm(list = ls())
diamond_dir <- "/Users/yangy/Documents/E/city/metagenome/HGT/diamond"
lgt_dir     <- "/Users/yangy/Documents/E/city/metagenome/HGT/hgt1"
library(circlize)

# 列出所有 .f6 文件
f6_files <- list.files(diamond_dir, pattern = "\\.f6$", full.names = TRUE)

# 最终结果
all_results <- data.frame()

for (f in f6_files) {
  
  # 提取样品名，例如 FZA2
  sample_id <- sub("_SARG.*", "", basename(f))
  message("Processing sample: ", sample_id)
  
  ### 1. 判断 f6 文件是否为空 ###
  if (file.info(f)$size == 0) {
    message("  Empty f6 file, writing zeros...")
    
    summary_sample <- data.frame(
      from = NA,
      to = NA,
      count = 0,
      sample = sample_id
    )
    
    all_results <- rbind(all_results, summary_sample)
    next
  }
  
  ### 2. f6 不是空文件再继续 ###
  df <- read.table(f, stringsAsFactors = FALSE)
  df$contig_clean <- sub("^([^_]+_[^_]+).*", "\\1", df$V1)
  
  # 找对应的 lgt 文件
  lgt_file <- file.path(lgt_dir, paste0(sample_id, ".lgt.tsv"))
  
  if (!file.exists(lgt_file)) {
    warning("No LGT file for: ", sample_id)
    
    summary_sample <- data.frame(
      from = NA,
      to = NA,
      count = 0,
      sample = sample_id
    )
    
    all_results <- rbind(all_results, summary_sample)
    next
  }
  
  # 读取 LGT 表
  lgt <- read.delim(lgt_file, stringsAsFactors = FALSE)
  colnames(lgt)[1] <- "contig"
  
  # 合并
  merged <- merge(df["contig_clean"], lgt,
                  by.x = "contig_clean", 
                  by.y = "contig")
  
  if (nrow(merged) == 0) {
    message("  No matched contigs, writing zeros...")
    
    summary_sample <- data.frame(
      from = NA,
      to = NA,
      count = 0,
      sample = sample_id
    )
    
    all_results <- rbind(all_results, summary_sample)
    next
  }
  
  ### 统计 ###
  colA <- colnames(lgt)[13]
  colB <- colnames(lgt)[14]
  
  merged$count <- 1
  
  summary_sample <- aggregate(count ~ merged[[colA]] + merged[[colB]], 
                              data = merged, sum)
  
  colnames(summary_sample) <- c("from", "to", "count")
  summary_sample$sample <- sample_id
  
  all_results <- rbind(all_results, summary_sample)
}

# 查看最终统计表
all_results
colnames(all_results)[4]<-"sample.id"
group<-read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T,row.names=1)

all_results2 <- merge(all_results, group[, c("sample.id", "group")],
                      by = "sample.id", 
                      all.x = TRUE)
all_results2 <- na.omit(all_results2)
all_results2$from_phylum <- sub(".*p__([^|]+).*", "\\1", all_results2$from)
all_results2$to_phylum <- sub(".*p__([^|]+).*", "\\1", all_results2$to)

phyla <- c(
  "Actinobacteria",
  "Proteobacteria",
  "Planctomycetes",
  "Acidobacteria",
  "Bacteroidetes",
  "Chloroflexi",
  "Firmicutes",
  "Gemmatimonadetes",
  "Spirochaetes",
  "Verrucomicrobia",
  "Cyanobacteria",
  "Deinococcus_Thermus"
)



all_results2_forest <- subset(all_results2, group == "Forest")
sum(all_results2_forest$from_phylum == "Proteobacteria" &
      all_results2_forest$to_phylum   == "Proteobacteria")
sum(all_results2_forest$from_phylum == "Actinobacteria" &
      all_results2_forest$to_phylum   == "Actinobacteria")
library(dplyr)

undirected_summary_forest <- all_results2_forest %>%
  mutate(
    # 生成无向 pair（例如 A-B / B-A 都变成 A|B）
    pair = ifelse(from_phylum < to_phylum,
                  paste(from_phylum, to_phylum, sep = "|"),
                  paste(to_phylum, from_phylum, sep = "|"))
  ) %>%
  group_by(pair) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  tidyr::separate(pair, into = c("phylum1", "phylum2"), sep = "\\|")

undirected_summary_forest$phylum1 <- paste0(undirected_summary_forest$phylum1,"1")
undirected_summary_forest$phylum2 <- paste0( undirected_summary_forest$phylum2,"2")

colnames(undirected_summary_forest)[1]<-"from"
colnames(undirected_summary_forest)[2]<-"to"
colnames(undirected_summary_forest)[3]<-"value"
chordDiagram(undirected_summary_forest)
undirected_summary_forest

grid.col <- c(
  "Actinobacteria1" = "#81b21f",
  "Actinobacteria2" = "#81b21f",
  "Proteobacteria2" = "#FF7F00",
  "Proteobacteria1" = "#FF7F00",
  "Firmicutes1" = "#f59b7b",
  "Firmicutes2" = "#f59b7b",
  "Bacteroidetes1" = "#fcc41e",
  "Bacteroidetes2" = "#fcc41e",
  "Acidobacteria1" = "#8d7cba",
  "Acidobacteria2" = "#8d7cba",
  "Chloroflexi1"="#bebfb4",
  "Chloroflexi2"="#bebfb4",
  "Planctomycetes1"="#65A2D2",
  "Planctomycetes2"="#65A2D2",
  "Gemmatimonadetes1"="pink",
  "Gemmatimonadetes2"="pink")
#chordDiagram(undirected_summary_forest, grid.col = grid.col, column.col = grid.col)
#pdf("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/HGT_foret弦图1.pdf", width = 5, height = 4)  # 设置 PDF 文件名与大小
chordDiagram(undirected_summary_forest, 
             grid.col = grid.col, 
             column.col = grid.col,
             transparency = 0.8)
#dev.off()   # 关闭绘图设备并保存 PDF


#cropland

all_results2_cropland <- subset(all_results2, group == "Cropland")
sum(all_results2_cropland$from_phylum == "Proteobacteria" &
      all_results2_cropland$to_phylum   == "Proteobacteria")
sum(all_results2_cropland$from_phylum == "Actinobacteria" &
      all_results2_cropland$to_phylum   == "Actinobacteria")



library(dplyr)

undirected_summary_cropland <- all_results2_cropland %>%
  mutate(
    # 生成无向 pair（例如 A-B / B-A 都变成 A|B）
    pair = ifelse(from_phylum < to_phylum,
                  paste(from_phylum, to_phylum, sep = "|"),
                  paste(to_phylum, from_phylum, sep = "|"))
  ) %>%
  group_by(pair) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  tidyr::separate(pair, into = c("phylum1", "phylum2"), sep = "\\|")

undirected_summary_cropland$phylum1 <- paste0(undirected_summary_cropland$phylum1,"1")
undirected_summary_cropland$phylum2 <- paste0( undirected_summary_cropland$phylum2,"2")

colnames(undirected_summary_cropland)[1]<-"from"
colnames(undirected_summary_cropland)[2]<-"to"
colnames(undirected_summary_cropland)[3]<-"value"

grid.col <- c(
  "Actinobacteria1" = "#81b21f",
  "Actinobacteria2" = "#81b21f",
  "Proteobacteria2" = "#FF7F00",
  "Proteobacteria1" = "#FF7F00",
  "Firmicutes1" = "#f59b7b",
  "Firmicutes2" = "#f59b7b",
  "Bacteroidetes1" = "#fcc41e",
  "Bacteroidetes2" = "#fcc41e",
  "Acidobacteria1" = "#8d7cba",
  "Acidobacteria2" = "#8d7cba",
  "Chloroflexi1"="#bebfb4",
  "Chloroflexi2"="#bebfb4",
  "Planctomycetes1"="#65A2D2",
  "Planctomycetes2"="#65A2D2",
  "Gemmatimonadetes1"="pink",
  "Gemmatimonadetes2"="pink",
  "Spirochaetes1"="#974c01",
  "Spirochaetes2"="#974c01")
#chordDiagram(undirected_summary_cropland, grid.col = grid.col, column.col = grid.col)
#pdf("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/HGT_cropland弦图1.pdf", width = 5, height = 4)  # 设置 PDF 文件名与大小
chordDiagram(undirected_summary_cropland, 
             grid.col = grid.col, 
             column.col = grid.col,
             transparency = 0.8)
#dev.off()   # 关闭绘图设备并保存 PDF
all_results2_industry <- subset(all_results2, group == "Industry")

sum(all_results2_industry$from_phylum == "Proteobacteria" &
      all_results2_industry$to_phylum   == "Proteobacteria")
sum(all_results2_industry$from_phylum == "Actinobacteria" &
      all_results2_industry$to_phylum   == "Actinobacteria")



library(dplyr)

undirected_summary_industry <- all_results2_industry %>%
  mutate(
    # 生成无向 pair（例如 A-B / B-A 都变成 A|B）
    pair = ifelse(from_phylum < to_phylum,
                  paste(from_phylum, to_phylum, sep = "|"),
                  paste(to_phylum, from_phylum, sep = "|"))
  ) %>%
  group_by(pair) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  tidyr::separate(pair, into = c("phylum1", "phylum2"), sep = "\\|")

undirected_summary_industry$phylum1 <- paste0(undirected_summary_industry$phylum1,"1")
undirected_summary_industry$phylum2 <- paste0( undirected_summary_industry$phylum2,"2")

colnames(undirected_summary_industry)[1]<-"from"
colnames(undirected_summary_industry)[2]<-"to"
colnames(undirected_summary_industry)[3]<-"value"
#chordDiagram(undirected_summary_industry)
#undirected_summary_industry



grid.col <- c(
  "Actinobacteria1" = "#81b21f",
  "Actinobacteria2" = "#81b21f",
  "Proteobacteria2" = "#FF7F00",
  "Proteobacteria1" = "#FF7F00",
  "Firmicutes1" = "#f59b7b",
  "Firmicutes2" = "#f59b7b",
  "Bacteroidetes1" = "#fcc41e",
  "Bacteroidetes2" = "#fcc41e",
  "Acidobacteria1" = "#8d7cba",
  "Acidobacteria2" = "#8d7cba",
  "Chloroflexi1"="#bebfb4",
  "Chloroflexi2"="#bebfb4",
  "Planctomycetes1"="#65A2D2",
  "Planctomycetes2"="#65A2D2",
  "Gemmatimonadetes1"="pink",
  "Gemmatimonadetes2"="pink",
  "Spirochaetes1"="#974c01",
  "Spirochaetes2"="#974c01",
  "Verrucomicrobia1"="#e377c2",
  "Verrucomicrobia2"="#e377c2",
  "Cyanobacteria1"="#9d1d1d",
  "Cyanobacteria2"="#9d1d1d",
  "Deinococcus_Thermus1"="#57a4fd",
  "Deinococcus_Thermus2"="#57a4fd")
#chordDiagram(undirected_summary_industry, grid.col = grid.col, column.col = grid.col)
#pdf("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/HGT_industry弦图1.pdf", width = 5, height = 4)  # 设置 PDF 文件名与大小
chordDiagram(undirected_summary_industry, 
             grid.col = grid.col, 
             column.col = grid.col,
             transparency = 0.8)
#dev.off()   # 关闭绘图设备并保存 PDF

all_results2_urban <- subset(all_results2, group == "Urban")

sum(all_results2_urban$from_phylum == "Proteobacteria" &
      all_results2_urban$to_phylum   == "Proteobacteria")
sum(all_results2_urban$from_phylum == "Actinobacteria" &
      all_results2_urban$to_phylum   == "Actinobacteria")

library(dplyr)

undirected_summary_urban <- all_results2_urban %>%
  mutate(
    # 生成无向 pair（例如 A-B / B-A 都变成 A|B）
    pair = ifelse(from_phylum < to_phylum,
                  paste(from_phylum, to_phylum, sep = "|"),
                  paste(to_phylum, from_phylum, sep = "|"))
  ) %>%
  group_by(pair) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  tidyr::separate(pair, into = c("phylum1", "phylum2"), sep = "\\|")

undirected_summary_urban$phylum1 <- paste0(undirected_summary_urban$phylum1,"1")
undirected_summary_urban$phylum2 <- paste0( undirected_summary_urban$phylum2,"2")

colnames(undirected_summary_urban)[1]<-"from"
colnames(undirected_summary_urban)[2]<-"to"
colnames(undirected_summary_urban)[3]<-"value"


grid.col <- c(
  "Actinobacteria1" = "#81b21f",
  "Actinobacteria2" = "#81b21f",
  "Proteobacteria2" = "#FF7F00",
  "Proteobacteria1" = "#FF7F00",
  "Firmicutes1" = "#f59b7b",
  "Firmicutes2" = "#f59b7b",
  "Bacteroidetes1" = "#fcc41e",
  "Bacteroidetes2" = "#fcc41e",
  "Acidobacteria1" = "#8d7cba",
  "Acidobacteria2" = "#8d7cba",
  "Chloroflexi1"="#bebfb4",
  "Chloroflexi2"="#bebfb4",
  "Planctomycetes1"="#65A2D2",
  "Planctomycetes2"="#65A2D2",
  "Gemmatimonadetes1"="pink",
  "Gemmatimonadetes2"="pink",
  "Spirochaetes1"="#974c01",
  "Spirochaetes2"="#974c01",
  "Verrucomicrobia1"="#e377c2",
  "Verrucomicrobia2"="#e377c2",
  "Cyanobacteria1"="#9d1d1d",
  "Cyanobacteria2"="#9d1d1d",
  "Deinococcus_Thermus1"="#57a4fd",
  "Deinococcus_Thermus2"="#57a4fd")
#chordDiagram(undirected_summary_urban, grid.col = grid.col, column.col = grid.col)
#pdf("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/HGT_urban弦图1.pdf", width = 5, height = 4)  # 设置 PDF 文件名与大小
chordDiagram(undirected_summary_urban, 
             grid.col = grid.col, 
             column.col = grid.col,
             transparency = 0.8)
#dev.off()   # 关闭绘图设备并保存 PDF

grid.col <- c(
  "Actinobacteria1" = "#81b21f",
  "Actinobacteria2" = "#81b21f",
  "Proteobacteria2" = "#FF7F00",
  "Proteobacteria1" = "#FF7F00",
  "Firmicutes1" = "#f59b7b",
  "Firmicutes2" = "#f59b7b",
  "Bacteroidetes1" = "#fcc41e",
  "Bacteroidetes2" = "#fcc41e",
  "Acidobacteria1" = "#8d7cba",
  "Acidobacteria2" = "#8d7cba",
  "Chloroflexi1"="#bebfb4",
  "Chloroflexi2"="#bebfb4",
  "Planctomycetes1"="#65A2D2",
  "Planctomycetes2"="#65A2D2",
  "Gemmatimonadetes1"="pink",
  "Gemmatimonadetes2"="pink",
  "Spirochaetes1"="#974c01",
  "Spirochaetes2"="#974c01",
  "Verrucomicrobia1"="#e377c2",
  "Verrucomicrobia2"="#e377c2",
  "Cyanobacteria1"="#9d1d1d",
  "Cyanobacteria2"="#9d1d1d",
  "Deinococcus_Thermus1"="#57a4fd",
  "Deinococcus_Thermus2"="#57a4fd"
)
phylum_names <- unique(gsub("[12]$", "", names(grid.col)))
legend_colors <- sapply(phylum_names, function(x) grid.col[paste0(x, "1")])
plot.new()

# 画图例
legend(
  "topleft",
  legend = phylum_names,
  fill = legend_colors,
  border = NA,
  bty = "n",
  cex = 1,
  title = "Phylum"
)

dev.off()


#g-j
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(dplyr)
files <- list.files(
  "/Users/yangy/Documents/E/city/metagenome/HGT/diamond",
  pattern = "f6$",
  full.names = TRUE
)

all_counts <- lapply(files, function(file) {
  
  fname <- basename(file)
  sample_name <- str_extract(fname, "^[^_]+")   # 提取 _ 前的样品名
  
  # 如果文件为空，返回 NULL
  if (file.size(file) == 0) {
    message("⚠️ 空文件跳过: ", fname)
    return(NULL)
  }
  
  df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # 有些时候文件可能有行但列数不够（损坏），也跳过
  if (ncol(df) < 2) {
    message("⚠️ 文件格式异常跳过: ", fname)
    return(NULL)
  }
  
  colnames(df)[2] <- "gene"
  
  df %>%
    count(gene, name = sample_name)
})

# 去掉 NULL 文件
all_counts <- all_counts[!sapply(all_counts, is.null)]

# 合并所有样本计数
count_matrix <- reduce(all_counts, full_join, by = "gene") %>%
  replace(is.na(.), 0)


##Forest气泡图
hgt_arg<-count_matrix
row.names(hgt_arg)<-hgt_arg[,1]
hgt_arg<-hgt_arg[,-1]
gene<-hgt_arg
gene<-as.data.frame(t(gene))
##Forest
group <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
gene$group <- group$group[ match(rownames(gene), group$sample.id) ]
gene_forest <- gene[ gene$group == "Forest", ]

gene_forest<-as.data.frame(t(gene_forest))
gene_forest<-gene_forest[-235,]
gene_forest1 <- as.data.frame(lapply(gene_forest, as.numeric))
row.names(gene_forest1)<-row.names(gene_forest)
gene_forest<-gene_forest1
gene_forest$sum <- rowSums(gene_forest)

gene_forest$SARG.Seq.ID<-row.names(gene_forest)
gene_forest<-gene_forest[,c(23,24)]


type<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/type.csv",sep=",",header=T)




merged_gene <- merge(gene_forest, type[, c("SARG.Seq.ID", "Type", "HMM.category")], 
                     by = "SARG.Seq.ID", all.x = TRUE)

library(dplyr)
merged_gene_no0 <- merged_gene %>% filter(sum != 0)
merged_gene<-merged_gene_no0
merged_gene <- merged_gene[complete.cases(merged_gene), ]




summary_by_category <- merged_gene %>%
  group_by(HMM.category) %>%
  summarise(
    Type = first(Type),  # 保留该类对应的一个 Type 值
    across(where(is.numeric), sum, na.rm = TRUE),
    count = n()
  )
summary_by_category <- summary_by_category %>%
  filter(HMM.category != "")
summary_by_category$Network<-"Network"
colnames(summary_by_category)[3]<-"Value"
library(dplyr)

summary_by_category <- summary_by_category %>%
  rename(Level1 = Network,     # 重命名 Network 为 Level1
         Level2 = Type) %>%    # 重命名 Type 为 Level2
  relocate(Level1, .before = 1) %>%     # 把 Level1 移到第一列
  relocate(Level2, .after = Level1)     # 把 Level2 移到第二列
colnames(summary_by_category)[3]<-"Level3"
colnames(summary_by_category)[4]<-"Value"
summary_by_category<-summary_by_category[,1:4]

df<-summary_by_category
library(ggplot2)
library(dplyr)#AGrammarofDataManipulation
library(tidyr)#TidyMessyData
library(ggraph)#AnImplementationofGrammarofGraphicsforGraphsandNetworks#AnImplementationofGrammarofGraphicsforGraphsandNetworks
library(tidygraph)#ATidyAPIforGraphManipulation#ATidyAPIforGraphManipulation
library(reshape2)
df<-df[,-1]
df_node<-data.frame()
for (i in 1:(length(df) - 1)) {
  group <- colnames(df)[i]  # 获取当前列的名称
  data <- df[c(i, 3)] %>% 
    group_by(Level = !!sym(group)) %>% 
    summarise(size = sum(Value),
              level = i)  # 使用列名称作为 level
  df_node <- bind_rows(df_node, data)
}
df2 <- df[1:2]
df2$group <- df2$Level2
df3 <-unique(melt(df2, id.vars = c("group"))[c(1,3)])
colnames(df3)[2] <- "Level"
df_node <- merge(df_node, unique(df3), by = "Level")

##构建边文件
df_edge <- data.frame()
for (i in 1:(length(df) - 2)) {
  from_col <- colnames(df)[i]
  to_col <- colnames(df)[i + 1]
  data <- df %>%
    distinct(from=!!sym(from_col), to=!!sym(to_col))
  df_edge <- bind_rows(df_edge, data)
}
df_graph <- tbl_graph(df_node, df_edge)

unique(df_node$group)
type_colors <- c(
  "aminoglycoside"                      = "#ffca18",
  "fosfomycin"                          = "blue",
  "multidrug"                           = "#3be8b0",
  "beta_lactam"                         = "#fc636b",
  "antibacterial_fatty_acid"            = "#ffd984",
  "polymyxin"                           = "#92d0b0",
  "bacitracin"                          = "#f7a4a9",
  "macrolide-lincosamide-streptogramin" = "#a5a7d3",
  "novobiocin"                          = "#a4c3c5",
  "other_peptide_antibiotics"           = "#e0c17c",
  "rifamycin"                           = "#e0a679",
  "tetracenomycin_C"                    = "green",
  "tetracycline"                        = "orange",
  "vancomycin"                          = "pink"
)


level_to_no_fill <- c(
  "aminoglycoside", "fosfomycin", "multidrug" , "beta_lactam", 
  "antibacterial_fatty_acid", "polymyxin", "bacitracin", "macrolide-lincosamide-streptogramin", "novobiocin", 
  "other_peptide_antibiotics","rifamycin","tetracenomycin_C","tetracycline" ,"vancomycin")
p <- ggraph(df_graph, layout = 'circlepack', weight = size) +
  geom_node_circle(
    aes(
      fill = ifelse(Level %in% level_to_no_fill, "transparent", group),
      alpha = ifelse(Level %in% level_to_no_fill, 0, 1)
    ),
    show.legend = FALSE,
    color = NA,
    linewidth = 0
  ) +
  scale_fill_manual(values = type_colors) +
  geom_node_text(
    aes(label = ifelse(size > 0.04, Level, NA), size = size),
    color = "black",
    fontface = "bold",
    family = "serif",
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")
p

#ggsave("/Users/yangy/Documents/D/city/arg/figure/beta_diversity/hgt_arg_compposition_forest.pdf",p,height=6,width=8)

#cropland
gene_cropland <- gene[ gene$group == "Cropland", ]

gene_cropland<-as.data.frame(t(gene_cropland))
gene_cropland<-gene_cropland[-235,]
gene_cropland1 <- as.data.frame(lapply(gene_cropland, as.numeric))
row.names(gene_cropland1)<-row.names(gene_cropland)
gene_cropland<-gene_cropland1
gene_cropland$sum <- rowSums(gene_cropland)

gene_cropland$SARG.Seq.ID<-row.names(gene_cropland)
gene_cropland<-gene_cropland[,c(26,27)]


type<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/type.csv",sep=",",header=T)




merged_gene <- merge(gene_cropland, type[, c("SARG.Seq.ID", "Type", "HMM.category")], 
                     by = "SARG.Seq.ID", all.x = TRUE)

library(dplyr)
merged_gene_no0 <- merged_gene %>% filter(sum != 0)
merged_gene<-merged_gene_no0
merged_gene <- merged_gene[complete.cases(merged_gene), ]




summary_by_category <- merged_gene %>%
  group_by(HMM.category) %>%
  summarise(
    Type = first(Type),  # 保留该类对应的一个 Type 值
    across(where(is.numeric), sum, na.rm = TRUE),
    count = n()
  )
summary_by_category <- summary_by_category %>%
  filter(HMM.category != "")
summary_by_category$Network<-"Network"
colnames(summary_by_category)[3]<-"Value"
library(dplyr)

summary_by_category <- summary_by_category %>%
  rename(Level1 = Network,     # 重命名 Network 为 Level1
         Level2 = Type) %>%    # 重命名 Type 为 Level2
  relocate(Level1, .before = 1) %>%     # 把 Level1 移到第一列
  relocate(Level2, .after = Level1)     # 把 Level2 移到第二列
colnames(summary_by_category)[3]<-"Level3"
colnames(summary_by_category)[4]<-"Value"
summary_by_category<-summary_by_category[,1:4]

df<-summary_by_category
library(ggplot2)
library(dplyr)#AGrammarofDataManipulation
library(tidyr)#TidyMessyData
library(ggraph)#AnImplementationofGrammarofGraphicsforGraphsandNetworks#AnImplementationofGrammarofGraphicsforGraphsandNetworks
library(tidygraph)#ATidyAPIforGraphManipulation#ATidyAPIforGraphManipulation
library(reshape2)
df<-df[,-1]
df_node<-data.frame()
for (i in 1:(length(df) - 1)) {
  group <- colnames(df)[i]  # 获取当前列的名称
  data <- df[c(i, 3)] %>% 
    group_by(Level = !!sym(group)) %>% 
    summarise(size = sum(Value),
              level = i)  # 使用列名称作为 level
  df_node <- bind_rows(df_node, data)
}
df2 <- df[1:2]
df2$group <- df2$Level2
df3 <-unique(melt(df2, id.vars = c("group"))[c(1,3)])
colnames(df3)[2] <- "Level"
df_node <- merge(df_node, unique(df3), by = "Level")

##构建边文件
df_edge <- data.frame()
for (i in 1:(length(df) - 2)) {
  from_col <- colnames(df)[i]
  to_col <- colnames(df)[i + 1]
  data <- df %>%
    distinct(from=!!sym(from_col), to=!!sym(to_col))
  df_edge <- bind_rows(df_edge, data)
}
df_graph <- tbl_graph(df_node, df_edge)
unique(df_node$group)



type_colors <- c(
  "fosfomycin"                          = "blue",
  "multidrug"                           = "#3be8b0",
  "polymyxin"                           = "#92d0b0",
  "rifamycin"= "#e0a679",
  "bacitracin"                          = "#f7a4a9",
  "beta_lactam"="#fc636b",
  "macrolide-lincosamide-streptogramin" = "#a5a7d3",
  "novobiocin"                          = "#a4c3c5",
  "tetracycline"                        = "orange",
  "pleuromutilin_tiamulin" ="#c4ff0e",
  "vancomycin"                          = "pink")
level_to_no_fill <- c(
  "fosfomycin", "multidrug", "polymyxin" , "rifamycin", 
  "bacitracin", "beta_lactam", "macrolide-lincosamide-streptogramin", "novobiocin", "tetracycline" , 
  "pleuromutilin_tiamulin","vancomycin")

p <- ggraph(df_graph, layout = 'circlepack', weight = size) +
  geom_node_circle(
    aes(
      fill = ifelse(Level %in% level_to_no_fill, "transparent", group),
      alpha = ifelse(Level %in% level_to_no_fill, 0, 1)
    ),
    show.legend = FALSE,
    color = NA,
    linewidth = 0
  ) +
  scale_fill_manual(values = type_colors) +
  geom_node_text(
    aes(label = ifelse(size > 0.04, Level, NA), size = size),
    color = "black",
    fontface = "bold",
    family = "serif",
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")
p
#ggsave("/Users/yangy/Documents/D/city/arg/figure/beta_diversity/hgt_arg_compposition_cropland.pdf",p,height=6,width=8)

gene_industry <- gene[ gene$group == "Industry", ]

gene_industry<-as.data.frame(t(gene_industry))
gene_industry<-gene_industry[-235,]
gene_industry1 <- as.data.frame(lapply(gene_industry, as.numeric))
row.names(gene_industry1)<-row.names(gene_industry)
gene_industry<-gene_industry1
gene_industry$sum <- rowSums(gene_industry)

gene_industry$SARG.Seq.ID<-row.names(gene_industry)
gene_industry<-gene_industry[,c(27,28)]


type<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/type.csv",sep=",",header=T)




merged_gene <- merge(gene_industry, type[, c("SARG.Seq.ID", "Type", "HMM.category")], 
                     by = "SARG.Seq.ID", all.x = TRUE)

library(dplyr)
merged_gene_no0 <- merged_gene %>% filter(sum != 0)
merged_gene<-merged_gene_no0
merged_gene <- merged_gene[complete.cases(merged_gene), ]
summary_by_category <- merged_gene %>%
  group_by(HMM.category) %>%
  summarise(
    Type = first(Type),  # 保留该类对应的一个 Type 值
    across(where(is.numeric), sum, na.rm = TRUE),
    count = n()
  )
summary_by_category <- summary_by_category %>%
  filter(HMM.category != "")
summary_by_category$Network<-"Network"
colnames(summary_by_category)[3]<-"Value"
library(dplyr)

summary_by_category <- summary_by_category %>%
  rename(Level1 = Network,     # 重命名 Network 为 Level1
         Level2 = Type) %>%    # 重命名 Type 为 Level2
  relocate(Level1, .before = 1) %>%     # 把 Level1 移到第一列
  relocate(Level2, .after = Level1)     # 把 Level2 移到第二列
colnames(summary_by_category)[3]<-"Level3"
colnames(summary_by_category)[4]<-"Value"
summary_by_category<-summary_by_category[,1:4]

df<-summary_by_category
library(ggplot2)
library(dplyr)#AGrammarofDataManipulation
library(tidyr)#TidyMessyData
library(ggraph)#AnImplementationofGrammarofGraphicsforGraphsandNetworks#AnImplementationofGrammarofGraphicsforGraphsandNetworks
library(tidygraph)#ATidyAPIforGraphManipulation#ATidyAPIforGraphManipulation
library(reshape2)
df<-df[,-1]
df_node<-data.frame()
for (i in 1:(length(df) - 1)) {
  group <- colnames(df)[i]  # 获取当前列的名称
  data <- df[c(i, 3)] %>% 
    group_by(Level = !!sym(group)) %>% 
    summarise(size = sum(Value),
              level = i)  # 使用列名称作为 level
  df_node <- bind_rows(df_node, data)
}
df2 <- df[1:2]
df2$group <- df2$Level2
df3 <-unique(melt(df2, id.vars = c("group"))[c(1,3)])
colnames(df3)[2] <- "Level"
df_node <- merge(df_node, unique(df3), by = "Level")

##构建边文件
df_edge <- data.frame()
for (i in 1:(length(df) - 2)) {
  from_col <- colnames(df)[i]
  to_col <- colnames(df)[i + 1]
  data <- df %>%
    distinct(from=!!sym(from_col), to=!!sym(to_col))
  df_edge <- bind_rows(df_edge, data)
}
df_graph <- tbl_graph(df_node, df_edge)
unique(df_node$group)

type_colors <- c(
  "fosfomycin"                          = "blue",
  "multidrug"                           = "#3be8b0",
  "polymyxin"                           = "#92d0b0",
  "rifamycin"= "#e0a679",
  "bacitracin"                          = "#f7a4a9",
  "bicyclomycin"="#b97a56",
  "beta_lactam"="#fc636b",
  "macrolide-lincosamide-streptogramin" = "#a5a7d3",
  "chloramphenicol"="#5e9ceb",
    "trimethoprim"="#C49C94",
    "edeine" ="#00A087",
  "other_peptide_antibiotics"           = "#e0c17c",
  "mupirocin"="#3C5488",
    "novobiocin" ="#a4c3c5",
    "pleuromutilin_tiamulin"="#c4ff0e",
    "puromycin"="#DC0000",
    "streptothricin" ="#D6EAF8",
    "sulfonamide"="#F0E442",
    "tetracycline" = "orange",
    "vancomycin" ="pink")
level_to_no_fill <- c(
  "fosfomycin", "multidrug", "polymyxin" , "rifamycin", 
  "bacitracin", "bicyclomycin","beta_lactam", "macrolide-lincosamide-streptogramin","chloramphenicol",
  "trimethoprim","edeine","other_peptide_antibiotics","mupirocin","novobiocin","pleuromutilin_tiamulin",
  "puromycin","streptothricin","sulfonamide", "tetracycline" , "vancomycin")

p <- ggraph(df_graph, layout = 'circlepack', weight = size) +
  geom_node_circle(
    aes(
      fill = ifelse(Level %in% level_to_no_fill, "transparent", group),
      alpha = ifelse(Level %in% level_to_no_fill, 0, 1)
    ),
    show.legend = FALSE,
    color = NA,
    linewidth = 0
  ) +
  scale_fill_manual(values = type_colors) +
  geom_node_text(
    aes(label = ifelse(size > 0.04, Level, NA), size = size),
    color = "black",
    fontface = "bold",
    family = "serif",
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")
p
#ggsave("/Users/yangy/Documents/D/city/arg/figure/beta_diversity/hgt_arg_compposition_industry.pdf",p,height=6,width=8)

gene_urban <- gene[ gene$group == "Urban", ]

gene_urban<-as.data.frame(t(gene_urban))
gene_urban<-gene_urban[-235,]
gene_urban1 <- as.data.frame(lapply(gene_urban, as.numeric))
row.names(gene_urban1)<-row.names(gene_urban)
gene_urban<-gene_urban1
gene_urban$sum <- rowSums(gene_urban)

gene_urban$SARG.Seq.ID<-row.names(gene_urban)
gene_urban<-gene_urban[,c(27,28)]


type<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/type.csv",sep=",",header=T)




merged_gene <- merge(gene_urban, type[, c("SARG.Seq.ID", "Type", "HMM.category")], 
                     by = "SARG.Seq.ID", all.x = TRUE)

library(dplyr)
merged_gene_no0 <- merged_gene %>% filter(sum != 0)
merged_gene<-merged_gene_no0
merged_gene <- merged_gene[complete.cases(merged_gene), ]




summary_by_category <- merged_gene %>%
  group_by(HMM.category) %>%
  summarise(
    Type = first(Type),  # 保留该类对应的一个 Type 值
    across(where(is.numeric), sum, na.rm = TRUE),
    count = n()
  )
summary_by_category <- summary_by_category %>%
  filter(HMM.category != "")
summary_by_category$Network<-"Network"
colnames(summary_by_category)[3]<-"Value"
library(dplyr)

summary_by_category <- summary_by_category %>%
  rename(Level1 = Network,     # 重命名 Network 为 Level1
         Level2 = Type) %>%    # 重命名 Type 为 Level2
  relocate(Level1, .before = 1) %>%     # 把 Level1 移到第一列
  relocate(Level2, .after = Level1)     # 把 Level2 移到第二列
colnames(summary_by_category)[3]<-"Level3"
colnames(summary_by_category)[4]<-"Value"
summary_by_category<-summary_by_category[,1:4]

df<-summary_by_category
library(ggplot2)
library(dplyr)#AGrammarofDataManipulation
library(tidyr)#TidyMessyData
library(ggraph)#AnImplementationofGrammarofGraphicsforGraphsandNetworks#AnImplementationofGrammarofGraphicsforGraphsandNetworks
library(tidygraph)#ATidyAPIforGraphManipulation#ATidyAPIforGraphManipulation
library(reshape2)
df<-df[,-1]
df_node<-data.frame()
for (i in 1:(length(df) - 1)) {
  group <- colnames(df)[i]  # 获取当前列的名称
  data <- df[c(i, 3)] %>% 
    group_by(Level = !!sym(group)) %>% 
    summarise(size = sum(Value),
              level = i)  # 使用列名称作为 level
  df_node <- bind_rows(df_node, data)
}
df2 <- df[1:2]
df2$group <- df2$Level2
df3 <-unique(melt(df2, id.vars = c("group"))[c(1,3)])
colnames(df3)[2] <- "Level"
df_node <- merge(df_node, unique(df3), by = "Level")

##构建边文件
df_edge <- data.frame()
for (i in 1:(length(df) - 2)) {
  from_col <- colnames(df)[i]
  to_col <- colnames(df)[i + 1]
  data <- df %>%
    distinct(from=!!sym(from_col), to=!!sym(to_col))
  df_edge <- bind_rows(df_edge, data)
}
df_graph <- tbl_graph(df_node, df_edge)
unique(df_node$group)

type_colors <- c(
  "aminoglycoside"                      = "#ffca18",
  "beta_lactam"="#fc636b",
  "fosfomycin"                          = "blue",
  "quinolone"="#E377C2",
  "polymyxin"                           = "#92d0b0",
  "rifamycin"= "#e0a679",
  "bacitracin"                          = "#f7a4a9",
  "chloramphenicol"="#5e9ceb",
  "macrolide-lincosamide-streptogramin" = "#a5a7d3",
  "multidrug"                           = "#3be8b0",
  "other_peptide_antibiotics"           = "#e0c17c",
  "novobiocin"="#a4c3c5",
  "tetracycline"                        = "orange",
  "pleuromutilin_tiamulin"="#c4ff0e",
  "puromycin"="#DC0000",
  "sulfonamide"="#F0E442",
  "vancomycin"                          = "pink"
)



p <- ggraph(df_graph, layout = 'circlepack', weight = size) +
  geom_node_circle(
    aes(
      fill = ifelse(Level %in% level_to_no_fill, "transparent", group),
      alpha = ifelse(Level %in% level_to_no_fill, 0, 1)
    ),
    show.legend = FALSE,
    color = NA,
    linewidth = 0
  ) +
  scale_fill_manual(values = type_colors) +
  geom_node_text(
    aes(label = ifelse(size > 0.04, Level, NA), size = size),
    color = "black",
    fontface = "bold",
    family = "serif",
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")
p

#ggsave("/Users/yangy/Documents/D/city/arg/figure/beta_diversity/hgt_arg_compposition_urban.pdf",p,height=6,width=8)

#k
#trasfered ARGs 差异气泡图##shu liang
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(dplyr)
files <- list.files(
  "/Users/yangy/Documents/E/city/metagenome/HGT/diamond",
  pattern = "f6$",
  full.names = TRUE
)

all_counts <- lapply(files, function(file) {
  
  fname <- basename(file)
  sample_name <- str_extract(fname, "^[^_]+")   # 提取 _ 前的样品名
  
  # 如果文件为空，返回 NULL
  if (file.size(file) == 0) {
    message("⚠️ 空文件跳过: ", fname)
    return(NULL)
  }
  
  df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # 有些时候文件可能有行但列数不够（损坏），也跳过
  if (ncol(df) < 2) {
    message("⚠️ 文件格式异常跳过: ", fname)
    return(NULL)
  }
  
  colnames(df)[2] <- "gene"
  
  df %>%
    count(gene, name = sample_name)
})

# 去掉 NULL 文件
all_counts <- all_counts[!sapply(all_counts, is.null)]

# 合并所有样本计数
count_matrix <- reduce(all_counts, full_join, by = "gene") %>%
  replace(is.na(.), 0)


hgt_arg<-count_matrix
row.names(hgt_arg)<-hgt_arg[,1]
hgt_arg<-hgt_arg[,-1]
gene<-hgt_arg




type<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/type.csv",sep=",",header=T)


gene_df <- gene %>%
  tibble::rownames_to_column(var = "gene")
type_df <- type %>%
  rename(gene = 2, Type = 3) %>%
  select(gene, Type)
gene_with_type <- gene_df %>%
  left_join(type_df, by = "gene")


gene_with_type <- gene_with_type %>%
  tibble::column_to_rownames(var = "gene")
gene_with_type <- gene_with_type %>%
  filter(!is.na(Type))



gene<-gene_with_type 
library(dplyr)

type_sum <- gene %>%
  group_by(Type) %>%
  summarise(across(where(is.numeric), sum))
hgt<-as.data.frame(type_sum)

row.names(hgt)<-hgt[,1]
hgt<-hgt[,-1]

arg<-t(hgt)
arg <- arg[, colSums(arg, na.rm = TRUE) != 0]

arg <- arg[order(rownames(arg)), ]

group<-read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T,row.names=1)


sample_order <- rownames(arg)
group_filtered <- group[rownames(group) %in% sample_order, , drop = FALSE]
group_ordered <- group_filtered[sample_order, , drop = FALSE]
group<-group_ordered

identical(rownames(group),rownames(arg))

#添加分组信息
arg = cbind(arg,group$group) %>% as.data.frame()
colnames(arg)[25] = "group"
arg$group<-as.factor(arg$group)
# 把 character 列全部转为 numeric
arg[] <- lapply(arg, function(x) {
  if (is.character(x)) as.numeric(x) else x
})

# 按照分组求平均值
data_sum <- aggregate(. ~ group, data = arg, FUN = sum)
data_plot <- gather(data_sum,key = "env",value = "value",-"group")
data_plot$group = factor(data_plot$group,levels = c("Forest","Cropland","Industry","Urban"))




arg<-as.data.frame(arg)
arg$group<-group$group
data_analy = gather(arg,key = "env",value = "value",-"group")#利用for循环进行批量单因素方差分析
library(agricolae)
slist = unique(data_analy$env) %>% as.vector()
slist
result = data.frame()

for(i in 1:length(slist)){  
  data1 = filter(data_analy,env == slist[i])  # 执行单因素方差分析  
  dt <- data1 %>% group_by(group) %>%   
    summarise(mean=mean(value), sd = sd(value))  
  anova_result <- aov(value ~ group, data = data1)  # 查看方差分析结果 
  summary(anova_result)  
  anova_data <- summary(anova_result)[[1]]#获取F、df（组内与组间）和P  
  anova_data = anova_data[-2,]  
  anova_data$env = data1$env[1]  #进行多重比较  
  duncan = duncan.test(anova_result,"group")#Duncan  
  le = duncan$groups  
  le$Groups = rownames(le)  
  le = arrange(le,Groups)  
  re = cbind(dt,le)  
  re = re[,-(6)]   
  re = re[,-(4)]  
  re$env = data1$env[1]  
  re$result = paste0(round(re$mean,2),"±",round(re$sd,2),re$groups) 
  result = rbind(result,re)}
result
data_plot$groups = result$groups
library(dplyr)

data_plot <- data_plot %>% filter(env != "Others")
library(dplyr)
library(tidyr)

result = data.frame()

for (i in 1:length(slist)) {
  data1 <- filter(data_analy, env == slist[i])
  
  # 计算每组的均值和标准差
  dt <- data1 %>%
    group_by(group) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  
  # Forest 组数据
  forest_vals <- data1 %>% filter(group == "Forest") %>% pull(value)
  
  # 初始化 p 值列
  dt$p_value <- NA
  dt$sig <- ""
  
  # 只比较非 Forest 的组 vs Forest
  for (j in 1:nrow(dt)) {
    current_group <- dt$group[j]
    if (current_group != "Forest") {
      test_vals <- data1 %>% filter(group == current_group) %>% pull(value)
      
      # 用 tryCatch 来避免警告和错误
      test_result <- tryCatch({
        wilcox.test(test_vals, forest_vals)
      }, error = function(e) {
        # 如果出错，返回一个 NA 的 p 值
        list(p.value = NA)
      })
      
      # 如果 p.value 不是 NA，则继续处理
      if (!is.na(test_result$p.value)) {
        dt$p_value[j] <- test_result$p.value
        
        # 添加显著性标记
        if (test_result$p.value < 0.001) {
          dt$sig[j] <- "***"
        } else if (test_result$p.value < 0.01) {
          dt$sig[j] <- "**"
        } else if (test_result$p.value < 0.05) {
          dt$sig[j] <- "*"
        } else {
          dt$sig[j] <- ""
        }
      } else {
        # 如果 p 值为 NA，显著性标记为空
        dt$p_value[j] <- NA
        dt$sig[j] <- ""
      }
    } else {
      dt$p_value[j] <- NA
      dt$sig[j] <- ""
    }
  }
  
  # 添加 env 信息
  dt$env <- slist[i]
  dt$result <- paste0(round(dt$mean, 2), "±", round(dt$sd, 2), dt$sig)
  
  result <- rbind(result, dt)
}
# 查看结果
result



library(dplyr)
# 1. 提取 Urban 和 Forest 组，并计算差值
diff_order <- result %>%
  filter(group %in% c("Urban", "Forest")) %>%
  select(group, env, mean) %>%
  pivot_wider(names_from = group, values_from = mean) %>%
  mutate(diff = Urban - Forest)

# 2. 根据 Urban - Forest 的符号拆分排序
pos_env <- diff_order %>% filter(diff > 0) %>% arrange(desc(diff)) %>% pull(env)
neg_env <- diff_order %>% filter(diff < 0) %>% arrange(desc(diff)) %>% pull(env)

# 3. 合并顺序，先 Urban > Forest，后 Urban < Forest
env_order <- c(pos_env, neg_env)

# 4. 设置 factor 顺序
data_plot$env <- factor(data_plot$env, levels = env_order)

# 5. 合并显著性标记
data_plot <- left_join(data_plot, result %>% select(group, env, sig), by = c("group", "env"))
data_plot$label <- data_plot$sig

# 6. 去掉 Others，并替换长标签
data_plot <- data_plot %>% filter(env != "Others") %>%
  mutate(env = str_replace_all(env, c(
    "Antimicrobial activity/Competitive advantage" = "Antimicrobial activity/\nCompetitive advantage",
    "Nutritional/Metabolic factor" = "Nutritional/\nMetabolic factor",
    "Post-translational modification" = "Post-translational \nmodification"
  )))

# 7. 设置 group 顺序
data_plot$group <- factor(data_plot$group, levels = c("Forest", "Cropland", "Industry", "Urban"))

# 5. 合并显著性标记
data_plot <- left_join(data_plot, result %>% select(group, env, sig), by = c("group", "env"))
data_plot$label <- data_plot$sig


data_plot <- data_plot %>%
  filter(env != "Others") %>%
  mutate(env = str_replace_all(env, c(
    "Antimicrobial activity/Competitive advantage" = "Antimicrobial activity/\nCompetitive advantage",
    "Nutritional/Metabolic factor" = "Nutritional/\nMetabolic factor",
    "Post-translational modification" = "Post-translational \nmodification"
  )))


# 重新设定 factor 顺序
data_plot$env <- factor(data_plot$env, levels = env_order)


# 7. 设置 group 顺序
data_plot$group <- factor(data_plot$group, levels = c("Forest", "Cropland", "Industry", "Urban"))
colnames(data_plot)[5]<-"label"


levels(data_plot$env) <- gsub(
  "macrolide-lincosamide-streptogramin",
  "MLS",
  levels(data_plot$env)
)
levels(data_plot$env) <- gsub(
  "other_peptide_antibiotics",
  "other_peptide\n_antibiotics",
  levels(data_plot$env)
)

levels(data_plot$env) <- gsub(
  "antibacterial_fatty_acid",
  "antibacterial_\nfatty_acid",
  levels(data_plot$env)
)
levels(data_plot$env) <- gsub(
  "pleuromutilin_tiamulin",
  "pleuromutilin\ntiamulin",
  levels(data_plot$env)
)
p <- ggplot(data_plot, aes(env, group)) +
  geom_point(aes(fill = value, size = value), color = "black", shape = 21) +
  scale_fill_gradient(
    low = "#c0dedf", high = "#11c1b1",
    name = "HGT count"     # ← 颜色图例标题
  ) +
  scale_size(
    range = c(1, 20),
    name = "HGT count"         # ← 点大小图例标题
  ) +
  geom_text(aes(label = label), size = 6, family = "sans",
            position = position_dodge(0), vjust = 0.5) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.ticks.length = unit(3, "mm"),
    legend.margin = margin(20, unit = 'pt')
  )
p

#ggsave("/Users/yangy/Documents/D/city/arg/figure/beta_diversity/hgt_arg_type_count_差异气泡图型号.pdf",p,height=5,width=20)






#trasfered ARGs 差异气泡图
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(dplyr)
files <- list.files(
  "/Users/yangy/Documents/E/city/metagenome/HGT/diamond",
  pattern = "f6$",
  full.names = TRUE
)

all_counts <- lapply(files, function(file) {
  
  fname <- basename(file)
  sample_name <- str_extract(fname, "^[^_]+")   # 提取 _ 前的样品名
  
  # 如果文件为空，返回 NULL
  if (file.size(file) == 0) {
    message("⚠️ 空文件跳过: ", fname)
    return(NULL)
  }
  
  df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # 有些时候文件可能有行但列数不够（损坏），也跳过
  if (ncol(df) < 2) {
    message("⚠️ 文件格式异常跳过: ", fname)
    return(NULL)
  }
  
  colnames(df)[2] <- "gene"
  
  df %>%
    count(gene, name = sample_name)
})

# 去掉 NULL 文件
all_counts <- all_counts[!sapply(all_counts, is.null)]

# 合并所有样本计数
count_matrix <- reduce(all_counts, full_join, by = "gene") %>%
  replace(is.na(.), 0)


hgt_arg<-count_matrix
row.names(hgt_arg)<-hgt_arg[,1]
hgt_arg<-hgt_arg[,-1]
gene<-hgt_arg

arg<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 

gene_mat <- as.matrix(gene)
arg_mat  <- as.matrix(arg)
common_genes  <- intersect(rownames(gene_mat), rownames(arg_mat))
common_samples <- intersect(colnames(gene_mat), colnames(arg_mat))
arg_sub <- arg_mat[common_genes, common_samples, drop = FALSE]
filled_gene <- matrix(
  0,
  nrow = nrow(gene_mat),
  ncol = ncol(gene_mat),
  dimnames = list(rownames(gene_mat), colnames(gene_mat))
)
filled_gene[common_genes, common_samples] <- arg_sub
filled_gene <- as.data.frame(filled_gene)
gene<-filled_gene
type<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/type.csv",sep=",",header=T)


gene_df <- gene %>%
  tibble::rownames_to_column(var = "gene")
type_df <- type %>%
  rename(gene = 2, Type = 3) %>%
  select(gene, Type)
gene_with_type <- gene_df %>%
  left_join(type_df, by = "gene")


gene_with_type <- gene_with_type %>%
  tibble::column_to_rownames(var = "gene")
gene_with_type <- gene_with_type %>%
  filter(!is.na(Type))



gene<-gene_with_type 
library(dplyr)

type_sum <- gene %>%
  group_by(Type) %>%
  summarise(across(where(is.numeric), sum))
hgt<-as.data.frame(type_sum)

row.names(hgt)<-hgt[,1]
hgt<-hgt[,-1]

arg<-t(hgt)
arg <- arg[, colSums(arg, na.rm = TRUE) != 0]

arg <- arg[order(rownames(arg)), ]
data_scale <- scale(arg) %>% as.data.frame()
group<-read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T,row.names=1)

sample_order <- rownames(data_scale)
group_filtered <- group[rownames(group) %in% sample_order, , drop = FALSE]
group_ordered <- group_filtered[sample_order, , drop = FALSE]
group<-group_ordered

identical(rownames(group),rownames(data_scale))

#添加分组信息
data_scale = cbind(data_scale,group$group) %>% as.data.frame()
colnames(data_scale)[24] = "group"
# 按照分组求平均值
data_mean <- aggregate(. ~ group, data = data_scale, FUN = mean)
data_plot <- gather(data_mean,key = "env",value = "value",-"group")
data_plot$group = factor(data_plot$group,levels = c("Forest","Cropland","Industry","Urban"))




arg<-as.data.frame(arg)
arg$group<-group$group
data_analy = gather(arg,key = "env",value = "value",-"group")#利用for循环进行批量单因素方差分析
library(agricolae)
slist = unique(data_analy$env) %>% as.vector()
slist
result = data.frame()

for(i in 1:length(slist)){  
  data1 = filter(data_analy,env == slist[i])  # 执行单因素方差分析  
  dt <- data1 %>% group_by(group) %>%   
    summarise(mean=mean(value), sd = sd(value))  
  anova_result <- aov(value ~ group, data = data1)  # 查看方差分析结果 
  summary(anova_result)  
  anova_data <- summary(anova_result)[[1]]#获取F、df（组内与组间）和P  
  anova_data = anova_data[-2,]  
  anova_data$env = data1$env[1]  #进行多重比较  
  duncan = duncan.test(anova_result,"group")#Duncan  
  le = duncan$groups  
  le$Groups = rownames(le)  
  le = arrange(le,Groups)  
  re = cbind(dt,le)  
  re = re[,-(6)]   
  re = re[,-(4)]  
  re$env = data1$env[1]  
  re$result = paste0(round(re$mean,2),"±",round(re$sd,2),re$groups) 
  result = rbind(result,re)}
result
data_plot$groups = result$groups
library(dplyr)

data_plot <- data_plot %>% filter(env != "Others")




library(dplyr)
library(tidyr)

result = data.frame()

for (i in 1:length(slist)) {
  data1 <- filter(data_analy, env == slist[i])
  
  # 计算每组的均值和标准差
  dt <- data1 %>%
    group_by(group) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  
  # Forest 组数据
  forest_vals <- data1 %>% filter(group == "Forest") %>% pull(value)
  
  # 初始化 p 值列
  dt$p_value <- NA
  dt$sig <- ""
  
  # 只比较非 Forest 的组 vs Forest
  for (j in 1:nrow(dt)) {
    current_group <- dt$group[j]
    if (current_group != "Forest") {
      test_vals <- data1 %>% filter(group == current_group) %>% pull(value)
      
      # 用 tryCatch 来避免警告和错误
      test_result <- tryCatch({
        wilcox.test(test_vals, forest_vals)
      }, error = function(e) {
        # 如果出错，返回一个 NA 的 p 值
        list(p.value = NA)
      })
      
      # 如果 p.value 不是 NA，则继续处理
      if (!is.na(test_result$p.value)) {
        dt$p_value[j] <- test_result$p.value
        
        # 添加显著性标记
        if (test_result$p.value < 0.001) {
          dt$sig[j] <- "***"
        } else if (test_result$p.value < 0.01) {
          dt$sig[j] <- "**"
        } else if (test_result$p.value < 0.05) {
          dt$sig[j] <- "*"
        } else {
          dt$sig[j] <- ""
        }
      } else {
        # 如果 p 值为 NA，显著性标记为空
        dt$p_value[j] <- NA
        dt$sig[j] <- ""
      }
    } else {
      dt$p_value[j] <- NA
      dt$sig[j] <- ""
    }
  }
  
  # 添加 env 信息
  dt$env <- slist[i]
  dt$result <- paste0(round(dt$mean, 2), "±", round(dt$sd, 2), dt$sig)
  
  result <- rbind(result, dt)
}

# 查看结果
result
library(dplyr)
# 1. 提取 Urban 和 Forest 组，并计算差值
diff_order <- result %>%
  filter(group %in% c("Urban", "Forest")) %>%
  select(group, env, mean) %>%
  pivot_wider(names_from = group, values_from = mean) %>%
  mutate(diff = Urban - Forest)

# 2. 根据 Urban - Forest 的符号拆分排序
pos_env <- diff_order %>% filter(diff > 0) %>% arrange(desc(diff)) %>% pull(env)
neg_env <- diff_order %>% filter(diff < 0) %>% arrange(desc(diff)) %>% pull(env)

# 3. 合并顺序，先 Urban > Forest，后 Urban < Forest
env_order <- c(pos_env, neg_env)
env_order<-c("polymyxin","multidrug","bacitracin","quinolone","fosfomycin","other_peptide_antibiotics","tetracenomycin_C","chloramphenicol","sulfonamide","bicyclomycin",
             "streptothricin","antibacterial_fatty_acid","pleuromutilin_tiamulin","mupirocin","tetracycline","puromycin" ,"aminoglycoside",
             "beta_lactam","novobiocin","vancomycin","macrolide-lincosamide-streptogramin","rifamycin")

# 4. 设置 factor 顺序
data_plot$env <- factor(data_plot$env, levels = env_order)

# 5. 合并显著性标记
data_plot <- left_join(data_plot, result %>% select(group, env, sig), by = c("group", "env"))
data_plot$label <- data_plot$sig

# 6. 去掉 Others，并替换长标签
data_plot <- data_plot %>% filter(env != "Others") %>%
  mutate(env = str_replace_all(env, c(
    "Antimicrobial activity/Competitive advantage" = "Antimicrobial activity/\nCompetitive advantage",
    "Nutritional/Metabolic factor" = "Nutritional/\nMetabolic factor",
    "Post-translational modification" = "Post-translational \nmodification"
  )))

# 7. 设置 group 顺序
data_plot$group <- factor(data_plot$group, levels = c("Forest", "Cropland", "Industry", "Urban"))

# 5. 合并显著性标记
data_plot <- left_join(data_plot, result %>% select(group, env, sig), by = c("group", "env"))
data_plot$label <- data_plot$sig


data_plot <- data_plot %>%
  filter(env != "Others") %>%
  mutate(env = str_replace_all(env, c(
    "Antimicrobial activity/Competitive advantage" = "Antimicrobial activity/\nCompetitive advantage",
    "Nutritional/Metabolic factor" = "Nutritional/\nMetabolic factor",
    "Post-translational modification" = "Post-translational \nmodification"
  )))

# 同样也要替换 env_order 中的标签（保持一致）
env_order <- str_replace_all(env_order, c(
  "Antimicrobial activity/Competitive advantage" = "Antimicrobial activity/\nCompetitive advantage",
  "Nutritional/Metabolic factor" = "Nutritional/\nMetabolic factor",
  "Post-translational modification" = "Post-translational \nmodification"
))

# 重新设定 factor 顺序
data_plot$env <- factor(data_plot$env, levels = env_order)


# 7. 设置 group 顺序
data_plot$group <- factor(data_plot$group, levels = c("Forest", "Cropland", "Industry", "Urban"))
colnames(data_plot)[5]<-"label"


levels(data_plot$env) <- gsub(
  "macrolide-lincosamide-streptogramin",
  "macrolide-lincosamide\n-streptogramin",
  levels(data_plot$env)
)
levels(data_plot$env) <- gsub(
  "other_peptide_antibiotics",
  "other_peptide\n_antibiotics",
  levels(data_plot$env)
)

levels(data_plot$env) <- gsub(
  "antibacterial_fatty_acid",
  "antibacterial_\nfatty_acid",
  levels(data_plot$env)
)
levels(data_plot$env) <- gsub(
  "pleuromutilin_tiamulin",
  "pleuromutilin\ntiamulin",
  levels(data_plot$env)
)
p <- ggplot(data_plot, aes(env, group)) +
  geom_point(aes(fill = value, size = value), color = "black", shape = 21) +
  scale_fill_gradient(low = "#f0fbfa", high = "#11c1b1") +
  scale_size(range = c(1, 16)) +
  geom_text(aes(label = label, y = group), size = 6, family = "sans",
            position = position_dodge(0), vjust = 0.5) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),   # 图面边框
    plot.border = element_rect(colour = "black", fill = NA, size = 1),    # 整张图外框
    axis.title.x = element_text(size = 16, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 14, colour = "black"),
    axis.text.y = element_text(face = "bold", size = 14, colour = "black"),
    axis.ticks.length = unit(3, "mm"),
    legend.margin = margin(20, unit = 'pt')
  ) +
  guides(
    size = guide_legend(title = "standardized value"),
    fill = guide_legend(title = expression("standardized value"))
  )

p

#ggsave("/Users/yangy/Documents/D/city/arg/figure/beta_diversity/hgt_arg_type_差异气泡图型号.pdf",p,height=5,width=18)









