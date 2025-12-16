#Fig3
#effect size on alpha diversity
library(metafor)
library(ggplot2)
library(glmulti)
library(betapart)
library(agricolae)
library(ggplot2)
library(patchwork)
library(vegan)
library(reshape2)
rm(list=ls())

otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E//city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]
identical(rownames(otu), group$sample.id)
group$bac_shannon <- vegan::diversity(otu, index = "shannon")


high <-read.csv("/Users/yangy/Documents/E//city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
high <-high [,-1]
high<-as.data.frame(t(high))
high$sample.id <- rownames(high)
high <- high[order(high$sample.id),]
identical(rownames(high), group$sample.id)
high<-high[,-74]
group$high_shannon <- vegan::diversity(high, index = "shannon")




arg<-read.csv("/Users/yangy/Documents/E//city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 
arg<-as.data.frame(t(arg))
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-4544]
group$arg_shannon <- vegan::diversity(arg, index = "shannon")

vfg<-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
vfg<-as.data.frame(t(vfg))
vfg$sample.id <- rownames(vfg)
vfg <- vfg[order(vfg$sample.id),]
identical(rownames(vfg), group$sample.id)
vfg<-vfg[,-5410]
group$vfg_shannon <- vegan::diversity(vfg, index = "shannon")

mge<-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
mge<-as.data.frame(t(mge))
mge$sample.id <- rownames(mge)
mge <- mge[order(mge$sample.id),]
identical(rownames(mge), group$sample.id)
mge<-mge[,-1432]
group$mge_shannon <- vegan::diversity(mge, index = "shannon")


data<-group
colnames(data)[3]<-"Group"
data$Group[data$Group == "Forest"] <- "CK"
colnames(data)[3]<-"TR"
library(tidyverse) 
library(metafor) 

# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(bac_shannon), # 计算生物量均值
            SD = sd(bac_shannon), 
            # 计算生物量标准差
            N = n()) # 计算样本量
treatment_data <- summary_data %>% filter(TR !="CK") # 处理组数据            
control_data <- summary_data %>% filter(TR =="CK") # 对照组数据


# 获取对照组的统计量
C_mean <- control_data$Mean# 对照组均值
C_sd <- control_data$SD# 对照组标准差
C_n <- control_data$N# 对照组样本量

# 为每个处理组添加对应的对照组统计量
d1 <- treatment_data %>%
  mutate(
    CMean = C_mean, # 添加对照组均值
    CSD = C_sd, # 添加对照组标准差
    CN = C_n # 添加对照组样本量 
  )
print(d1)
# 计算效应值 (响应比 Response Ratio)# ROM = Log Response Ratio，对数响应比
d2 <- escalc(measure ="ROM",
             data = d1,
             m1i = Mean, sd1i = SD, n1i = N, # 处理组参数
             m2i = CMean, sd2i = CSD, n2i = CN)# 对照组参数
print(d2)
# 1. 计算整体效应值（随机效应模型）
overall_effect <- rma(yi, vi, data = d2, method = "REML")
summary(overall_effect)
# 2. 计算各处理组的效应值
subgroup_effects <- list()
for(tr in unique(d2$TR)) {
  subgroup_data <- d2 %>% filter(TR == tr)
  subgroup_effects[[tr]] <- rma(yi, vi, data = subgroup_data, method = "REML")
}
# 3. 整理结果用于绘图
results_df <- data.frame(
  Subgroup = character(),
  estimate = numeric(),
  se = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)
# 添加各处理组结果
for(tr in names(subgroup_effects)) {
  result <- subgroup_effects[[tr]]
  results_df <- rbind(results_df, data.frame(
    Subgroup = tr,
    estimate = result$b[1],   # 效应值估计
    se = result$se[1],        # 标准误
    pval = result$pval[1]     # p值
  ))
}
# 添加整体效应结果
results_df <- rbind(results_df, data.frame(
  Subgroup = "Overall",
  estimate = overall_effect$b[1],
  se = overall_effect$se[1],
  pval = overall_effect$pval[1]
))
# 添加显著性标识
results_df$sig <- ifelse(results_df$pval < 0.001, "***",
                         ifelse(results_df$pval < 0.01, "**",
                                ifelse(results_df$pval < 0.05, "*", "")))
print(results_df)

library(ggplot2)
library(ggthemes) 

plot_data <- results_df
plot_data <- plot_data %>%
  filter(Subgroup != "Overall")

plot_data$Subgroup <- ifelse(plot_data$Subgroup == "Cropland", "Cropland vs Forest",
                             ifelse(plot_data$Subgroup == "Industry", "Industry vs Forest",
                                    ifelse(plot_data$Subgroup == "Urban", "Urban vs Forest", plot_data$Subgroup)))



plot_data$index <- factor(plot_data$Subgroup, levels = rev(c("Cropland vs Forest", "Industry vs Forest", "Urban vs Forest")))
# 判断显著性
plot_data$significant <- ifelse(plot_data$pval < 0.05, "Significant", "Non-significant")
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)
# 为显著与非显著设置填充色
plot_data$fill_color <- ifelse(
  plot_data$significant == "Significant",
  custom_colors[plot_data$Subgroup],
  "white"
)
plot_data$point_shape <- ifelse(plot_data$significant == "Significant", 16, 21)
plot_data
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)

plot_data_bacteria<-plot_data
plot_data_bacteria$group<-"bacteria"




# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(arg_shannon), # 计算生物量均值
            SD = sd(arg_shannon), 
            # 计算生物量标准差
            N = n()) # 计算样本量
treatment_data <- summary_data %>% filter(TR !="CK") # 处理组数据            
control_data <- summary_data %>% filter(TR =="CK") # 对照组数据


# 获取对照组的统计量
C_mean <- control_data$Mean# 对照组均值
C_sd <- control_data$SD# 对照组标准差
C_n <- control_data$N# 对照组样本量

# 为每个处理组添加对应的对照组统计量
d1 <- treatment_data %>%
  mutate(
    CMean = C_mean, # 添加对照组均值
    CSD = C_sd, # 添加对照组标准差
    CN = C_n # 添加对照组样本量 
  )
print(d1)
# 计算效应值 (响应比 Response Ratio)# ROM = Log Response Ratio，对数响应比
d2 <- escalc(measure ="ROM",
             data = d1,
             m1i = Mean, sd1i = SD, n1i = N, # 处理组参数
             m2i = CMean, sd2i = CSD, n2i = CN)# 对照组参数
print(d2)
# 1. 计算整体效应值（随机效应模型）
overall_effect <- rma(yi, vi, data = d2, method = "REML")
summary(overall_effect)
# 2. 计算各处理组的效应值
subgroup_effects <- list()
for(tr in unique(d2$TR)) {
  subgroup_data <- d2 %>% filter(TR == tr)
  subgroup_effects[[tr]] <- rma(yi, vi, data = subgroup_data, method = "REML")
}
# 3. 整理结果用于绘图
results_df <- data.frame(
  Subgroup = character(),
  estimate = numeric(),
  se = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)
# 添加各处理组结果
for(tr in names(subgroup_effects)) {
  result <- subgroup_effects[[tr]]
  results_df <- rbind(results_df, data.frame(
    Subgroup = tr,
    estimate = result$b[1],   # 效应值估计
    se = result$se[1],        # 标准误
    pval = result$pval[1]     # p值
  ))
}
# 添加整体效应结果
results_df <- rbind(results_df, data.frame(
  Subgroup = "Overall",
  estimate = overall_effect$b[1],
  se = overall_effect$se[1],
  pval = overall_effect$pval[1]
))
# 添加显著性标识
results_df$sig <- ifelse(results_df$pval < 0.001, "***",
                         ifelse(results_df$pval < 0.01, "**",
                                ifelse(results_df$pval < 0.05, "*", "")))
print(results_df)

library(ggplot2)
library(ggthemes) 

plot_data <- results_df
plot_data <- plot_data %>%
  filter(Subgroup != "Overall")

plot_data$Subgroup <- ifelse(plot_data$Subgroup == "Cropland", "Cropland vs Forest",
                             ifelse(plot_data$Subgroup == "Industry", "Industry vs Forest",
                                    ifelse(plot_data$Subgroup == "Urban", "Urban vs Forest", plot_data$Subgroup)))



plot_data$index <- factor(plot_data$Subgroup, levels = rev(c("Cropland vs Forest", "Industry vs Forest", "Urban vs Forest")))
# 判断显著性
plot_data$significant <- ifelse(plot_data$pval < 0.05, "Significant", "Non-significant")
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)
# 为显著与非显著设置填充色
plot_data$fill_color <- ifelse(
  plot_data$significant == "Significant",
  custom_colors[plot_data$Subgroup],
  "white"
)
plot_data$point_shape <- ifelse(plot_data$significant == "Significant", 16, 21)
plot_data
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)

plot_data_arg<-plot_data
plot_data_arg$group<-"arg"





# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(high_shannon), # 计算生物量均值
            SD = sd(high_shannon), 
            # 计算生物量标准差
            N = n()) # 计算样本量
treatment_data <- summary_data %>% filter(TR !="CK") # 处理组数据            
control_data <- summary_data %>% filter(TR =="CK") # 对照组数据


# 获取对照组的统计量
C_mean <- control_data$Mean# 对照组均值
C_sd <- control_data$SD# 对照组标准差
C_n <- control_data$N# 对照组样本量

# 为每个处理组添加对应的对照组统计量
d1 <- treatment_data %>%
  mutate(
    CMean = C_mean, # 添加对照组均值
    CSD = C_sd, # 添加对照组标准差
    CN = C_n # 添加对照组样本量 
  )
print(d1)
# 计算效应值 (响应比 Response Ratio)# ROM = Log Response Ratio，对数响应比
d2 <- escalc(measure ="ROM",
             data = d1,
             m1i = Mean, sd1i = SD, n1i = N, # 处理组参数
             m2i = CMean, sd2i = CSD, n2i = CN)# 对照组参数
print(d2)
# 1. 计算整体效应值（随机效应模型）
overall_effect <- rma(yi, vi, data = d2, method = "REML")
summary(overall_effect)
# 2. 计算各处理组的效应值
subgroup_effects <- list()
for(tr in unique(d2$TR)) {
  subgroup_data <- d2 %>% filter(TR == tr)
  subgroup_effects[[tr]] <- rma(yi, vi, data = subgroup_data, method = "REML")
}
# 3. 整理结果用于绘图
results_df <- data.frame(
  Subgroup = character(),
  estimate = numeric(),
  se = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)
# 添加各处理组结果
for(tr in names(subgroup_effects)) {
  result <- subgroup_effects[[tr]]
  results_df <- rbind(results_df, data.frame(
    Subgroup = tr,
    estimate = result$b[1],   # 效应值估计
    se = result$se[1],        # 标准误
    pval = result$pval[1]     # p值
  ))
}
# 添加整体效应结果
results_df <- rbind(results_df, data.frame(
  Subgroup = "Overall",
  estimate = overall_effect$b[1],
  se = overall_effect$se[1],
  pval = overall_effect$pval[1]
))
# 添加显著性标识
results_df$sig <- ifelse(results_df$pval < 0.001, "***",
                         ifelse(results_df$pval < 0.01, "**",
                                ifelse(results_df$pval < 0.05, "*", "")))
print(results_df)

library(ggplot2)
library(ggthemes) 

plot_data <- results_df
plot_data <- plot_data %>%
  filter(Subgroup != "Overall")

plot_data$Subgroup <- ifelse(plot_data$Subgroup == "Cropland", "Cropland vs Forest",
                             ifelse(plot_data$Subgroup == "Industry", "Industry vs Forest",
                                    ifelse(plot_data$Subgroup == "Urban", "Urban vs Forest", plot_data$Subgroup)))



plot_data$index <- factor(plot_data$Subgroup, levels = rev(c("Cropland vs Forest", "Industry vs Forest", "Urban vs Forest")))
# 判断显著性
plot_data$significant <- ifelse(plot_data$pval < 0.05, "Significant", "Non-significant")
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)
# 为显著与非显著设置填充色
plot_data$fill_color <- ifelse(
  plot_data$significant == "Significant",
  custom_colors[plot_data$Subgroup],
  "white"
)
plot_data$point_shape <- ifelse(plot_data$significant == "Significant", 16, 21)
plot_data
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)

plot_data_high<-plot_data
plot_data_high$group<-"high"




##vfg
# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(vfg_shannon), # 计算生物量均值
            SD = sd(vfg_shannon), 
            # 计算生物量标准差
            N = n()) # 计算样本量
treatment_data <- summary_data %>% filter(TR !="CK") # 处理组数据            
control_data <- summary_data %>% filter(TR =="CK") # 对照组数据


# 获取对照组的统计量
C_mean <- control_data$Mean# 对照组均值
C_sd <- control_data$SD# 对照组标准差
C_n <- control_data$N# 对照组样本量

# 为每个处理组添加对应的对照组统计量
d1 <- treatment_data %>%
  mutate(
    CMean = C_mean, # 添加对照组均值
    CSD = C_sd, # 添加对照组标准差
    CN = C_n # 添加对照组样本量 
  )
print(d1)
# 计算效应值 (响应比 Response Ratio)# ROM = Log Response Ratio，对数响应比
d2 <- escalc(measure ="ROM",
             data = d1,
             m1i = Mean, sd1i = SD, n1i = N, # 处理组参数
             m2i = CMean, sd2i = CSD, n2i = CN)# 对照组参数
print(d2)
# 1. 计算整体效应值（随机效应模型）
overall_effect <- rma(yi, vi, data = d2, method = "REML")
summary(overall_effect)
# 2. 计算各处理组的效应值
subgroup_effects <- list()
for(tr in unique(d2$TR)) {
  subgroup_data <- d2 %>% filter(TR == tr)
  subgroup_effects[[tr]] <- rma(yi, vi, data = subgroup_data, method = "REML")
}
# 3. 整理结果用于绘图
results_df <- data.frame(
  Subgroup = character(),
  estimate = numeric(),
  se = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)
# 添加各处理组结果
for(tr in names(subgroup_effects)) {
  result <- subgroup_effects[[tr]]
  results_df <- rbind(results_df, data.frame(
    Subgroup = tr,
    estimate = result$b[1],   # 效应值估计
    se = result$se[1],        # 标准误
    pval = result$pval[1]     # p值
  ))
}
# 添加整体效应结果
results_df <- rbind(results_df, data.frame(
  Subgroup = "Overall",
  estimate = overall_effect$b[1],
  se = overall_effect$se[1],
  pval = overall_effect$pval[1]
))
# 添加显著性标识
results_df$sig <- ifelse(results_df$pval < 0.001, "***",
                         ifelse(results_df$pval < 0.01, "**",
                                ifelse(results_df$pval < 0.05, "*", "")))
print(results_df)

library(ggplot2)
library(ggthemes) 

plot_data <- results_df
plot_data <- plot_data %>%
  filter(Subgroup != "Overall")

plot_data$Subgroup <- ifelse(plot_data$Subgroup == "Cropland", "Cropland vs Forest",
                             ifelse(plot_data$Subgroup == "Industry", "Industry vs Forest",
                                    ifelse(plot_data$Subgroup == "Urban", "Urban vs Forest", plot_data$Subgroup)))



plot_data$index <- factor(plot_data$Subgroup, levels = rev(c("Cropland vs Forest", "Industry vs Forest", "Urban vs Forest")))
# 判断显著性
plot_data$significant <- ifelse(plot_data$pval < 0.05, "Significant", "Non-significant")
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)
# 为显著与非显著设置填充色
plot_data$fill_color <- ifelse(
  plot_data$significant == "Significant",
  custom_colors[plot_data$Subgroup],
  "white"
)
plot_data$point_shape <- ifelse(plot_data$significant == "Significant", 16, 21)
plot_data
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)

plot_data_vfg<-plot_data
plot_data_vfg$group<-"vfg"

#mge
# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(mge_shannon), # 计算生物量均值
            SD = sd(mge_shannon), 
            # 计算生物量标准差
            N = n()) # 计算样本量
treatment_data <- summary_data %>% filter(TR !="CK") # 处理组数据            
control_data <- summary_data %>% filter(TR =="CK") # 对照组数据


# 获取对照组的统计量
C_mean <- control_data$Mean# 对照组均值
C_sd <- control_data$SD# 对照组标准差
C_n <- control_data$N# 对照组样本量

# 为每个处理组添加对应的对照组统计量
d1 <- treatment_data %>%
  mutate(
    CMean = C_mean, # 添加对照组均值
    CSD = C_sd, # 添加对照组标准差
    CN = C_n # 添加对照组样本量 
  )
print(d1)
# 计算效应值 (响应比 Response Ratio)# ROM = Log Response Ratio，对数响应比
d2 <- escalc(measure ="ROM",
             data = d1,
             m1i = Mean, sd1i = SD, n1i = N, # 处理组参数
             m2i = CMean, sd2i = CSD, n2i = CN)# 对照组参数
print(d2)
# 1. 计算整体效应值（随机效应模型）
overall_effect <- rma(yi, vi, data = d2, method = "REML")
summary(overall_effect)
# 2. 计算各处理组的效应值
subgroup_effects <- list()
for(tr in unique(d2$TR)) {
  subgroup_data <- d2 %>% filter(TR == tr)
  subgroup_effects[[tr]] <- rma(yi, vi, data = subgroup_data, method = "REML")
}
# 3. 整理结果用于绘图
results_df <- data.frame(
  Subgroup = character(),
  estimate = numeric(),
  se = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)
# 添加各处理组结果
for(tr in names(subgroup_effects)) {
  result <- subgroup_effects[[tr]]
  results_df <- rbind(results_df, data.frame(
    Subgroup = tr,
    estimate = result$b[1],   # 效应值估计
    se = result$se[1],        # 标准误
    pval = result$pval[1]     # p值
  ))
}
# 添加整体效应结果
results_df <- rbind(results_df, data.frame(
  Subgroup = "Overall",
  estimate = overall_effect$b[1],
  se = overall_effect$se[1],
  pval = overall_effect$pval[1]
))
# 添加显著性标识
results_df$sig <- ifelse(results_df$pval < 0.001, "***",
                         ifelse(results_df$pval < 0.01, "**",
                                ifelse(results_df$pval < 0.05, "*", "")))
print(results_df)

library(ggplot2)
library(ggthemes) 

plot_data <- results_df
plot_data <- plot_data %>%
  filter(Subgroup != "Overall")

plot_data$Subgroup <- ifelse(plot_data$Subgroup == "Cropland", "Cropland vs Forest",
                             ifelse(plot_data$Subgroup == "Industry", "Industry vs Forest",
                                    ifelse(plot_data$Subgroup == "Urban", "Urban vs Forest", plot_data$Subgroup)))



plot_data$index <- factor(plot_data$Subgroup, levels = rev(c("Cropland vs Forest", "Industry vs Forest", "Urban vs Forest")))
# 判断显著性
plot_data$significant <- ifelse(plot_data$pval < 0.05, "Significant", "Non-significant")
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)
# 为显著与非显著设置填充色
plot_data$fill_color <- ifelse(
  plot_data$significant == "Significant",
  custom_colors[plot_data$Subgroup],
  "white"
)
plot_data$point_shape <- ifelse(plot_data$significant == "Significant", 16, 21)
plot_data
custom_colors <- c(
  "Cropland vs Forest" = "#fdb218",       
  "Industry vs Forest" = "#5e9ceb", 
  "Urban vs Forest" = "#ff7f27"
)

plot_data_mge<-plot_data
plot_data_mge$group<-"mge"






plot<-rbind(plot_data_bacteria,plot_data_arg,plot_data_high,plot_data_vfg,plot_data_mge)


##urban effect
plot_urban <- plot %>%
  dplyr::filter(Subgroup == "Urban vs Forest")
unique(plot_urban$group)

custom_colors <- c(
  "bacteria" = "#FA9A38",       
  "arg" = "#6B30C5", 
  "high" = "#097369",
  "vfg"="#A13B46",
  "mge"="#777777"
)

plot_urban$group <- factor(plot_urban$group,
                           levels = rev(c("bacteria", "arg", "high","vfg","mge")))


format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 科学计数法
  } else {
    sprintf("%.3f", round(p, 3))          # 三位小数
  }
}

plot_urban$sig1 <- sprintf(
  "estimate = %.3f, P = %s",
  round(plot_urban$estimate, 3),
  sapply(plot_urban$pval, format_p)
)



p <- ggplot(plot_urban, aes(x = estimate, y = group)) +
  
  # ✅ ARG行浅灰色背景
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1.5, ymax = 2.5,
           fill = "grey90") +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 3.5, ymax = 4.5,
           fill = "grey90") +
  # 95% CI
  geom_errorbarh(aes(
    xmax = estimate + 1.96 * se,
    xmin = estimate - 1.96 * se,
    color = group,
    linetype = "95% CI"
  ),
  size = 1.5, height = 0, alpha = 0.4) +
  
  # SE
  geom_errorbarh(aes(
    xmax = estimate + se,
    xmin = estimate - se,
    color = group,
    linetype = "±SE"
  ),
  size = 2.5, height = 0, alpha = 0.7) +
  
  geom_point(aes(
    fill = fill_color,
    color = group,
    shape = significant
  ),
  size = 6, stroke = 1) +
  
  geom_text(aes(label = sig1), x = -1, size = 6, hjust = 0,nudge_y = -0.15 ) +
  
  scale_color_manual(values = custom_colors, name = "") +
  scale_fill_identity() +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Non-significant" = 21)) +
  scale_linetype_manual(name = "Error bars",
                        values = c("95% CI" = "solid", "±SE" = "solid"),
                        guide = guide_legend(override.aes = list(size = c(1.5, 2.5),
                                                                 alpha = c(0.4,0.7),
                                                                 color = "black"))) +
  
  scale_y_discrete(labels = c("bacteria" = "Bacterial community",
                              "arg" = "Resistome",
                              "high" = "High risk resistome",
                              "vfg"="Virulome",
                              "mge"="Mobilome")) +
  
  scale_x_continuous(
    limits = c(-2, 2),
    labels = function(x) sprintf("%.1f", x)
  ) +
  
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 0.6) +
  
  labs(x = ' ', y = '', title = "Urban effect on H'") +
  theme_few() +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    axis.title.x = element_blank()
  )

p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/effect_size_alpha_diversity_urban_effect_vfg_mge版本.pdf",p,width=8,height=6)



###industry
plot_industry <- plot %>%
  dplyr::filter(Subgroup == "Industry vs Forest")
unique(plot_industry$group)

custom_colors <- c(
  "bacteria" = "#FA9A38",       
  "arg" = "#6B30C5", 
  "high" = "#097369",
  "vfg"="#A13B46",
  "mge"="#777777"
)
plot_industry$group <- factor(plot_industry$group,
                              levels = rev(c("bacteria", "arg", "high","vfg","mge")))

plot_industry$sig1 <- sprintf(
  "estimate = %.3f, P = %s",
  round(plot_industry$estimate, 3),
  sapply(plot_industry$pval, format_p)
)






p <- ggplot(plot_industry, aes(x = estimate, y = group)) +
  
  # ✅ ARG行浅灰色背景
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1.5, ymax = 2.5,
           fill = "grey90") +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 3.5, ymax = 4.5,
           fill = "grey90") +
  # 95% CI
  geom_errorbarh(aes(
    xmax = estimate + 1.96 * se,
    xmin = estimate - 1.96 * se,
    color = group,
    linetype = "95% CI"
  ),
  size = 1.5, height = 0, alpha = 0.4) +
  
  # SE
  geom_errorbarh(aes(
    xmax = estimate + se,
    xmin = estimate - se,
    color = group,
    linetype = "±SE"
  ),
  size = 2.5, height = 0, alpha = 0.7) +
  
  geom_point(aes(
    fill = fill_color,
    color = group,
    shape = significant
  ),
  size = 6, stroke = 1) +
  
  geom_text(aes(label = sig1), x = -1, size = 6, hjust = 0,nudge_y = -0.15 ) +
  
  scale_color_manual(values = custom_colors, name = "") +
  scale_fill_identity() +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Non-significant" = 21)) +
  scale_linetype_manual(name = "Error bars",
                        values = c("95% CI" = "solid", "±SE" = "solid"),
                        guide = guide_legend(override.aes = list(size = c(1.5, 2.5),
                                                                 alpha = c(0.4,0.7),
                                                                 color = "black"))) +
  
  scale_y_discrete(labels = c("bacteria" = "Bacterial community",
                              "arg" = "Resistome",
                              "high" = "High risk resistome",
                              "vfg"="Virulome",
                              "mge"="Mobilome")) +
  
  scale_x_continuous(
    limits = c(-2, 2),
    labels = function(x) sprintf("%.1f", x)
  ) +
  
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 0.6) +
  
  labs(x = ' ', y = '', title = "Industry effect on H'") +
  theme_few() +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    axis.title.x = element_blank()
  )

p





#ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/effect_size_alpha_diversity_industry_effect_vfg_mge版本.pdf",p,width=8,height=6)


##cropland
plot_cropland <- plot %>%
  dplyr::filter(Subgroup == "Cropland vs Forest")
unique(plot_cropland$group)

custom_colors <- c(
  "bacteria" = "#FA9A38",       
  "arg" = "#6B30C5", 
  "high" = "#097369",
  "vfg"="#A13B46",
  "mge"="#777777"
)
plot_cropland$group <- factor(plot_cropland$group,
                              levels = rev(c("bacteria", "arg", "high","vfg","mge")))

plot_cropland$sig1 <- sprintf(
  "estimate = %.3f, P = %s",
  round(plot_cropland$estimate, 3),
  sapply(plot_cropland$pval, format_p)
)


p <- ggplot(plot_cropland, aes(x = estimate, y = group)) +
  
  # ✅ ARG行浅灰色背景
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1.5, ymax = 2.5,
           fill = "grey90") +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 3.5, ymax = 4.5,
           fill = "grey90") +
  
  # 95% CI
  geom_errorbarh(aes(
    xmax = estimate + 1.96 * se,
    xmin = estimate - 1.96 * se,
    color = group,
    linetype = "95% CI"
  ),
  size = 1.5, height = 0, alpha = 0.4) +
  
  # SE
  geom_errorbarh(aes(
    xmax = estimate + se,
    xmin = estimate - se,
    color = group,
    linetype = "±SE"
  ),
  size = 2.5, height = 0, alpha = 0.7) +
  
  geom_point(aes(
    fill = fill_color,
    color = group,
    shape = significant
  ),
  size = 6, stroke = 1) +
  
  geom_text(aes(label = sig1), x = -1, size = 6, hjust = 0,nudge_y = -0.15 ) +
  
  scale_color_manual(values = custom_colors, name = "") +
  scale_fill_identity() +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Non-significant" = 21)) +
  scale_linetype_manual(name = "Error bars",
                        values = c("95% CI" = "solid", "±SE" = "solid"),
                        guide = guide_legend(override.aes = list(size = c(1.5, 2.5),
                                                                 alpha = c(0.4,0.7),
                                                                 color = "black"))) +
  
  scale_y_discrete(labels = c("bacteria" = "Bacterial community",
                              "arg" = "Resistome",
                              "high" = "High risk resistome",
                              "vfg"="Virulome",
                              "mge"="Mobilome")) +
  
  scale_x_continuous(
    limits = c(-2, 2),
    labels = function(x) sprintf("%.1f", x)
  ) +
  
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 0.6) +
  
  labs(x = ' ', y = '', title = "Cropland effect on H'") +
  theme_few() +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    axis.title.x = element_blank()
  )

p
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/effect_size_alpha_diversity_cropland_effect_vfg_版本.pdf",p,width=8,height=6)


#d-h heatmap 
library(ggplot2)
library(patchwork)
library(vegan)
library(ggplot2)
library(picante)
library(agricolae)
library(dplyr)
library(agricolae)
library(FSA)

# 清空环境
rm(list=ls())

##计算bac的Shannon
otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]
identical(rownames(otu), group$sample.id)
group$shannon_bac <- vegan::diversity(otu, index = "shannon")

##计算ARG的Shannon
arg<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 
arg<-as.data.frame(t(arg))


arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-4544]
group$shannon_arg <- vegan::diversity(arg, index = "shannon")

##计算high_risk_arg的Shannon
arg <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
arg<-arg[,-1]
arg<-as.data.frame(t(arg))

arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)

arg<-arg[,-74]
group$shannon_high <- vegan::diversity(arg, index = "shannon")

#计算vfg的Shannon
vfg<-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
vfg<-as.data.frame(t(vfg))
vfg$sample.id <- rownames(vfg)
vfg <- vfg[order(vfg$sample.id),]
identical(rownames(vfg), group$sample.id)
vfg<-vfg[,-5410]
group$shannon_vfg <- vegan::diversity(vfg, index = "shannon")

mge<-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
mge<-as.data.frame(t(mge))
mge$sample.id <- rownames(mge)
mge <- mge[order(mge$sample.id),]
identical(rownames(mge), group$sample.id)
mge<-mge[,-1432]
group$shannon_mge <- vegan::diversity(mge, index = "shannon")




##env与bac_Shannon,ARG_shannon和high_risk_arg_Shannon之间的相关性

environment_vars <- colnames(group)[5:19]  # 环境变量列从第4列开始
library(psych)

environment_vars <- c("pH", "Latitude", "Longitude", "MAT", "MAP", "Aridity.index", 
                      "Per.capita.GDP", "population.density", "PM2.5", "SO2", 
                      "PM10", "NO2", "CO", "O3","HDI")



microbiome_vars <- c("shannon_bac", "shannon_arg", "shannon_high","shannon_vfg","shannon_mge")
results <- data.frame()

for (group_name in unique(group$group)) {
  group_data <- subset(group, group == group_name)
  
  # 提取微生物和环境变量矩阵
  micro_mat <- group_data[, microbiome_vars]
  env_mat   <- group_data[, environment_vars]
  
  # 一次计算所有 Spearman 相关
  cor_res <- corr.test(micro_mat, env_mat, use="pairwise",method="spearman",adjust="fdr",alpha=.05,ci=FALSE)
  
  # 把结果矩阵“展开”成长格式（tidy）
  res_df <- reshape2::melt(cor_res$r, varnames = c("Microbial_Feature", "Env"), value.name = "Spearman")
  p_df   <- reshape2::melt(cor_res$p, varnames = c("Microbial_Feature", "Env"), value.name = "P_value")
  
  # 合并 r 和 p
  res_all <- merge(res_df, p_df, by = c("Microbial_Feature", "Env"))
  res_all$Group <- group_name
  
  results <- rbind(results, res_all)
}




results$sig <- ""
results$sig[results$P_value < 0.05 & results$P_value >= 0.01] <- "*"
results$sig[results$P_value < 0.01 & results$P_value >= 0.001] <- "**"
results$sig[results$P_value < 0.001] <- "***"
library(dplyr)
# 显著性列
results$significance <- ifelse(is.na(results$sig), "", results$sig)
results1<- results %>% filter(Group != "Total")

# 排序分组
results1$Group <- factor(results1$Group, levels = c("Forest", "Cropland", "Industry", "Urban"))

title_mapping <- list(
  "Latitude" = "Latitude",
  "Longitude" = "Longitude",
  "MAT" = "Mean annual\ntemperature",   # ← 改这里
  "MAP" = "Mean annual\nprecipitation",
  "Aridity.index"="Aridity index",
  "Per.capita.GDP" = "GDP per capita",
  "population.density" = "Population density",
  "PM2.5" = "Air PM2.5",
  "SO2" = "Air SO2",
  "PM10" = "Air PM10",
  "NO2" = "Air NO2",
  "CO"="Air CO",
  "O3"="Air O3",
  "HDI"="HDI",
  "pH"="Soil pH"
)



results1$Env_label <- unlist(title_mapping[as.character(results1$Env)])
# Apply recoding using case_when
results1 <- results1 %>%
  mutate(Microbial_Feature = case_when(
    Microbial_Feature == "shannon_bac" ~ "Bacterial diversity",
    Microbial_Feature == "shannon_arg" ~ "Resistome diversity",
    Microbial_Feature == "shannon_high" ~ "High risk resistome diversity",
    Microbial_Feature == "shannon_vfg" ~ "VFGs diversity",
    Microbial_Feature == "shannon_mge" ~ "MGEs diversity",
    TRUE ~ Microbial_Feature  # Keep other values unchanged
  ))
results1$Microbial_Feature <- factor(
  results1$Microbial_Feature,
  levels = c("Bacterial diversity","Resistome diversity","High risk resistome diversity","VFGs diversity","MGEs diversity")
)
#构建两行标签（仅对显著性点）
results1$label_text <- ifelse(
  results1$significance != "",
  paste0(round(results1$Spearman, 2), "\n", results1$significance),
  NA
)
facet_labels <- c(
  "Bacterial diversity" = "Bacterial H'",
  "Resistome diversity" = "Resistome H'",
  "High risk resistome diversity" = "High risk resistome H'",
  "VFGs diversity"="Virulome H'",
  "MGEs diversity"="Mobilome H'"# 这里换行
)

results1$Env_label <- factor(
  results1$Env_label,
  levels = c(
    "HDI",
    "Population density",
    "GDP per capita",
    "Air O3",
    "Air CO",
    "Air NO2",
    "Air PM10",
    "Air SO2",
    "Air PM2.5",
    "Soil pH",
    "Aridity index",
    "Mean annual\nprecipitation",
    "Mean annual\ntemperature",  # ← 这里也要换行
    "Longitude",
    "Latitude"
  )
)



p_all <- ggplot(results1, aes(x = Group, y = Env_label, fill = Spearman)) + 
  geom_tile() +
  geom_text(aes(label = label_text), size = 5, na.rm = TRUE, lineheight = 0.6) +
  scale_fill_gradientn(colors = c('#00B554', 'white', '#B586CA'), limits = c(-1, 1)) +
  facet_wrap(~Microbial_Feature, scales = "fixed", ncol = 6, strip.position = "top",
             labeller = labeller(Microbial_Feature = facet_labels)) +
  labs(fill = "Spearman ρ", title = " ", x = NULL, y = NULL) +
  theme(
    strip.text = element_text(angle = 0, size = 14, face = "bold"),
    strip.placement = "outside",
    axis.text.x = element_text(hjust = 1, color = "black", face = "bold", size = 12, angle =60),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.spacing = unit(1, "lines"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

p_all

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/shannon_diversity/shannon与env关系加上vfg和mge版本.pdf",p_all, width = 15, height = 7, device = cairo_pdf)

#ipH difference and diversity difference
##difference 
library(ggplot2)
library(patchwork)
library(vegan)
library(ggplot2)
library(picante)
library(agricolae)
library(dplyr)
library(agricolae)
library(FSA)

# 清空环境
rm(list=ls())
otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]
identical(rownames(otu), group$sample.id)
group$shannon <- vegan::diversity(otu, index = "shannon")
format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 使用科学计数法
  } else {
    round(p, 3)
  }
}


calculate_pH_change <- function(group, Location_label) {
  # Filter samples for the selected locations
  sam_alp_filtered <- filter(group, Location %in% Location_label)
  
  # Split data into groups based on the 'group' column
  alp_forest <- filter(sam_alp_filtered, group %in% "Forest")
  alp_cropland <- filter(sam_alp_filtered, group %in% "Cropland")
  alp_industry <- filter(sam_alp_filtered, group %in% "Industry")
  alp_urban <- filter(sam_alp_filtered, group %in% "Urban")
  
  # Create an empty data frame to store all comparison results
  all_comparisons <- data.frame()
  
  # Calculate the average genome size for the Forest group in each Location
  forest_mean_size <- alp_forest %>%
    group_by(Location) %>%
    summarise(mean_forest_size = mean(pH, na.rm = TRUE))
  
  # For each Location, compare the mean of Forest group to each other group
  for (location in unique(sam_alp_filtered$Location)) {
    forest_size <- forest_mean_size$mean_forest_size[forest_mean_size$Location == location]
    
    # Forest vs Cropland (use original values for Cropland)
    diff_cropland <- alp_cropland %>%
      filter(Location == location) %>%
      mutate(diff = pH - forest_size)
    all_comparisons <- rbind(all_comparisons, 
                             data.frame(Location = location,
                                        group_comparison = "forest_cropland",
                                        pH_diff = diff_cropland$diff,
                                        pH_mean_forest_size = forest_size))
    
    # Forest vs Industry (use original values for Industry)
    diff_industry <- alp_industry %>%
      filter(Location == location) %>%
      mutate(diff = pH - forest_size)
    all_comparisons <- rbind(all_comparisons, 
                             data.frame(Location = location,
                                        group_comparison = "forest_industry",
                                        pH_diff = diff_industry$diff,
                                        pH_mean_forest_size = forest_size))
    
    # Forest vs Urban (use original values for Urban)
    diff_urban <- alp_urban %>%
      filter(Location == location) %>%
      mutate(diff = pH - forest_size)
    all_comparisons <- rbind(all_comparisons, 
                             data.frame(Location = location,
                                        group_comparison = "forest_urban",
                                        pH_diff = diff_urban$diff,
                                        pH_mean_forest_size = forest_size))
  }
  
  # Sort by 'group_comparison' (if you want to sort by a different column, change the name of the column)
  all_comparisons <- all_comparisons %>%
    arrange(group_comparison)  # Sorting by group_comparison (you can change to another column if needed)
  
  # Return the final comparison results data frame without calculating change
  return(all_comparisons)
}


Beijing_pH_change <- calculate_pH_change(group, "Beijing")
Changchun_pH_change <- calculate_pH_change(group, "Changchun")
Changsha_pH_change <- calculate_pH_change(group, "Changsha")
Chengdu_pH_change <- calculate_pH_change(group, "Chengdu")
Fuzhou_pH_change <- calculate_pH_change(group, "Fuzhou")
Guangzhou_pH_change <- calculate_pH_change(group, "Guangzhou")
Guiyang_pH_change <- calculate_pH_change(group, "Guiyang")
Hangzhou_pH_change <- calculate_pH_change(group, "Hangzhou")
Harbin_pH_change <- calculate_pH_change(group, "Harbin")
Hefei_pH_change <- calculate_pH_change(group, "Hefei")
Hohhot_pH_change <- calculate_pH_change(group, "Hohhot")
Jinan_pH_change <- calculate_pH_change(group, "Jinan")
Kunming_pH_change <- calculate_pH_change(group, "Kunming")
Lanzhou_pH_change <- calculate_pH_change(group, "Lanzhou")
Nanjing_pH_change <- calculate_pH_change(group, "Nanjing")
Nanning_pH_change <- calculate_pH_change(group, "Nanning")
Shanghai_pH_change <- calculate_pH_change(group, "Shanghai")
Shenyang_pH_change <- calculate_pH_change(group, "Shenyang")
Shijiazhuang_pH_change <- calculate_pH_change(group, "Shijiazhuang")
Taiyuan_pH_change <- calculate_pH_change(group, "Taiyuan")
Tianjin_pH_change <- calculate_pH_change(group, "Tianjin")
Urumqi_pH_change <- calculate_pH_change(group, "Urumqi")
Wuhan_pH_change <- calculate_pH_change(group, "Wuhan")
Xian_pH_change <- calculate_pH_change(group, "Xian")
Xining_pH_change <- calculate_pH_change(group, "Xining")
yinchuan_pH_change <- calculate_pH_change(group, "Yinchuan")
Zhengzhou_pH_change <- calculate_pH_change(group, "Zhengzhou")

pH_change<-rbind(Beijing_pH_change,Changchun_pH_change,Changsha_pH_change,
                 Chengdu_pH_change, Fuzhou_pH_change,Guangzhou_pH_change,Guiyang_pH_change
                 ,Hangzhou_pH_change,Harbin_pH_change,Hefei_pH_change,Hohhot_pH_change,
                 Jinan_pH_change,Kunming_pH_change,Lanzhou_pH_change,Nanjing_pH_change,
                 Nanning_pH_change,Shanghai_pH_change,Shenyang_pH_change,Shijiazhuang_pH_change,
                 Taiyuan_pH_change,Tianjin_pH_change,Urumqi_pH_change,
                 Wuhan_pH_change,Xian_pH_change,Xining_pH_change,yinchuan_pH_change,Zhengzhou_pH_change)
colnames(pH_change)[2]<-"group"
pH_change$group <- gsub("forest_cropland", "Cropland vs Forest", pH_change$group)
pH_change$group <- gsub("forest_industry", "Industry vs Forest", pH_change$group)
pH_change$group <- gsub("forest_urban", "Urban vs Forest", pH_change$group)
colnames(pH_change)[3]<-"pH"



##计算shannon change
calculate_shannon_change <- function(df, Location_label) {
  # Filter samples for the selected locations
  sam_alp_filtered <- dplyr::filter(df, Location %in% Location_label)
  
  # Split data into groups based on the 'group' column
  alp_forest <- dplyr::filter(sam_alp_filtered, group %in% "Forest")
  alp_cropland <- dplyr::filter(sam_alp_filtered, group %in% "Cropland")
  alp_industry <- dplyr::filter(sam_alp_filtered, group %in% "Industry")
  alp_urban <- dplyr::filter(sam_alp_filtered, group %in% "Urban")
  
  # Create an empty data frame to store all comparison results
  all_comparisons <- data.frame()
  
  # Calculate the average genome size for the Forest group in each Location
  forest_mean_size <- alp_forest %>%
    dplyr::group_by(Location) %>%
    dplyr::summarise(mean_forest_size = mean(shannon, na.rm = TRUE))
  
  # For each Location, compare the mean of Forest group to each other group
  for (location in unique(sam_alp_filtered$Location)) {
    forest_size <- forest_mean_size$mean_forest_size[forest_mean_size$Location == location]
    
    # Forest vs Cropland
    diff_cropland <- alp_cropland %>%
      dplyr::filter(Location == location) %>%
      dplyr::mutate(diff = shannon - forest_size)
    all_comparisons <- rbind(all_comparisons,
                             data.frame(Location = location,
                                        group_comparison = "forest_cropland",
                                        diff = diff_cropland$diff,
                                        mean_forest_size = forest_size))
    
    # Forest vs Industry
    diff_industry <- alp_industry %>%
      dplyr::filter(Location == location) %>%
      dplyr::mutate(diff = shannon - forest_size)
    all_comparisons <- rbind(all_comparisons,
                             data.frame(Location = location,
                                        group_comparison = "forest_industry",
                                        diff = diff_industry$diff,
                                        mean_forest_size = forest_size))
    
    # Forest vs Urban
    diff_urban <- alp_urban %>%
      dplyr::filter(Location == location) %>%
      dplyr::mutate(diff = shannon - forest_size)
    all_comparisons <- rbind(all_comparisons,
                             data.frame(Location = location,
                                        group_comparison = "forest_urban",
                                        diff = diff_urban$diff,
                                        mean_forest_size = forest_size))
  }
  
  # Sort by group_comparison
  all_comparisons <- all_comparisons %>%
    dplyr::arrange(group_comparison)
  
  return(all_comparisons)
}

Beijing_shannon_change <- calculate_shannon_change(group, "Beijing")
Changchun_shannon_change <- calculate_shannon_change(group, "Changchun")
Changsha_shannon_change <- calculate_shannon_change(group, "Changsha")
Chengdu_shannon_change <- calculate_shannon_change(group, "Chengdu")
Fuzhou_shannon_change <- calculate_shannon_change(group, "Fuzhou")
Guangzhou_shannon_change <- calculate_shannon_change(group, "Guangzhou")
Guiyang_shannon_change <- calculate_shannon_change(group, "Guiyang")
Hangzhou_shannon_change <- calculate_shannon_change(group, "Hangzhou")
Harbin_shannon_change <- calculate_shannon_change(group, "Harbin")
Hefei_shannon_change <- calculate_shannon_change(group, "Hefei")
Hohhot_shannon_change <- calculate_shannon_change(group, "Hohhot")
Jinan_shannon_change <- calculate_shannon_change(group, "Jinan")
Kunming_shannon_change <- calculate_shannon_change(group, "Kunming")
Lanzhou_shannon_change <- calculate_shannon_change(group, "Lanzhou")
Nanjing_shannon_change <- calculate_shannon_change(group, "Nanjing")
Nanning_shannon_change <- calculate_shannon_change(group, "Nanning")
Shanghai_shannon_change <- calculate_shannon_change(group, "Shanghai")
Shenyang_shannon_change <- calculate_shannon_change(group, "Shenyang")
Shijiazhuang_shannon_change <- calculate_shannon_change(group, "Shijiazhuang")
Taiyuan_shannon_change <- calculate_shannon_change(group, "Taiyuan")
Tianjin_shannon_change <- calculate_shannon_change(group, "Tianjin")
Urumqi_shannon_change <- calculate_shannon_change(group, "Urumqi")
Wuhan_shannon_change <- calculate_shannon_change(group, "Wuhan")
Xian_shannon_change <- calculate_shannon_change(group, "Xian")
Xining_shannon_change <- calculate_shannon_change(group, "Xining")
yinchuan_shannon_change <- calculate_shannon_change(group, "Yinchuan")
Zhengzhou_shannon_change <- calculate_shannon_change(group, "Zhengzhou")

bac_shannon_change<-rbind(Beijing_shannon_change,Changchun_shannon_change,Changsha_shannon_change,
                          Chengdu_shannon_change, Fuzhou_shannon_change,Guangzhou_shannon_change,Guiyang_shannon_change
                          ,Hangzhou_shannon_change,Harbin_shannon_change,Hefei_shannon_change,Hohhot_shannon_change,
                          Jinan_shannon_change,Kunming_shannon_change,Lanzhou_shannon_change,Nanjing_shannon_change,
                          Nanning_shannon_change,Shanghai_shannon_change,Shenyang_shannon_change,Shijiazhuang_shannon_change,
                          Taiyuan_shannon_change,Tianjin_shannon_change,Urumqi_shannon_change,
                          Wuhan_shannon_change,Xian_shannon_change,Xining_shannon_change,yinchuan_shannon_change,Zhengzhou_shannon_change)

colnames(bac_shannon_change)[2]<-"group"
bac_shannon_change$group<-as.factor(bac_shannon_change$group)

bac_shannon_change$group <- gsub("forest_cropland", "Cropland vs Forest", bac_shannon_change$group)
bac_shannon_change$group <- gsub("forest_industry", "Industry vs Forest", bac_shannon_change$group)
bac_shannon_change$group <- gsub("forest_urban", "Urban vs Forest", bac_shannon_change$group)
colnames(bac_shannon_change)[3]<-"shannon"

bac_pH_shannon<-cbind(bac_shannon_change,pH_change)
colnames(bac_pH_shannon)[2]<-"group"
bac_pH_shannon<-bac_pH_shannon[,-c(5,6)]
bac_pH_shannon$group1<-"bacteria"

##计算arg 的shannon difference
arg<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-4544]
group$shannon <- vegan::diversity(arg, index = "shannon")


Beijing_shannon_change <- calculate_shannon_change(group, "Beijing")
Changchun_shannon_change <- calculate_shannon_change(group, "Changchun")
Changsha_shannon_change <- calculate_shannon_change(group, "Changsha")
Chengdu_shannon_change <- calculate_shannon_change(group, "Chengdu")
Fuzhou_shannon_change <- calculate_shannon_change(group, "Fuzhou")
Guangzhou_shannon_change <- calculate_shannon_change(group, "Guangzhou")
Guiyang_shannon_change <- calculate_shannon_change(group, "Guiyang")
Hangzhou_shannon_change <- calculate_shannon_change(group, "Hangzhou")
Harbin_shannon_change <- calculate_shannon_change(group, "Harbin")
Hefei_shannon_change <- calculate_shannon_change(group, "Hefei")
Hohhot_shannon_change <- calculate_shannon_change(group, "Hohhot")
Jinan_shannon_change <- calculate_shannon_change(group, "Jinan")
Kunming_shannon_change <- calculate_shannon_change(group, "Kunming")
Lanzhou_shannon_change <- calculate_shannon_change(group, "Lanzhou")
Nanjing_shannon_change <- calculate_shannon_change(group, "Nanjing")
Nanning_shannon_change <- calculate_shannon_change(group, "Nanning")
Shanghai_shannon_change <- calculate_shannon_change(group, "Shanghai")
Shenyang_shannon_change <- calculate_shannon_change(group, "Shenyang")
Shijiazhuang_shannon_change <- calculate_shannon_change(group, "Shijiazhuang")
Taiyuan_shannon_change <- calculate_shannon_change(group, "Taiyuan")
Tianjin_shannon_change <- calculate_shannon_change(group, "Tianjin")
Urumqi_shannon_change <- calculate_shannon_change(group, "Urumqi")
Wuhan_shannon_change <- calculate_shannon_change(group, "Wuhan")
Xian_shannon_change <- calculate_shannon_change(group, "Xian")
Xining_shannon_change <- calculate_shannon_change(group, "Xining")
yinchuan_shannon_change <- calculate_shannon_change(group, "Yinchuan")
Zhengzhou_shannon_change <- calculate_shannon_change(group, "Zhengzhou")

arg_shannon_change<-rbind(Beijing_shannon_change,Changchun_shannon_change,Changsha_shannon_change,
                          Chengdu_shannon_change, Fuzhou_shannon_change,Guangzhou_shannon_change,Guiyang_shannon_change
                          ,Hangzhou_shannon_change,Harbin_shannon_change,Hefei_shannon_change,Hohhot_shannon_change,
                          Jinan_shannon_change,Kunming_shannon_change,Lanzhou_shannon_change,Nanjing_shannon_change,
                          Nanning_shannon_change,Shanghai_shannon_change,Shenyang_shannon_change,Shijiazhuang_shannon_change,
                          Taiyuan_shannon_change,Tianjin_shannon_change,Urumqi_shannon_change,
                          Wuhan_shannon_change,Xian_shannon_change,Xining_shannon_change,yinchuan_shannon_change,Zhengzhou_shannon_change)


colnames(arg_shannon_change)[2]<-"group"
arg_shannon_change$group<-as.factor(arg_shannon_change$group)

arg_shannon_change$group <- gsub("forest_cropland", "Cropland vs Forest", arg_shannon_change$group)
arg_shannon_change$group <- gsub("forest_industry", "Industry vs Forest", arg_shannon_change$group)
arg_shannon_change$group <- gsub("forest_urban", "Urban vs Forest", arg_shannon_change$group)
colnames(arg_shannon_change)[3]<-"shannon"

arg_pH_shannon<-cbind(arg_shannon_change,pH_change)
colnames(arg_pH_shannon)[2]<-"group"
arg_pH_shannon<-arg_pH_shannon[,-c(5,6)]
arg_pH_shannon$group1<-"arg"

##计算high
arg <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
arg<-arg[,-1]
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)

arg<-arg[,-74]
group$shannon <- vegan::diversity(arg, index = "shannon")

Beijing_shannon_change <- calculate_shannon_change(group, "Beijing")
Changchun_shannon_change <- calculate_shannon_change(group, "Changchun")
Changsha_shannon_change <- calculate_shannon_change(group, "Changsha")
Chengdu_shannon_change <- calculate_shannon_change(group, "Chengdu")
Fuzhou_shannon_change <- calculate_shannon_change(group, "Fuzhou")
Guangzhou_shannon_change <- calculate_shannon_change(group, "Guangzhou")
Guiyang_shannon_change <- calculate_shannon_change(group, "Guiyang")
Hangzhou_shannon_change <- calculate_shannon_change(group, "Hangzhou")
Harbin_shannon_change <- calculate_shannon_change(group, "Harbin")
Hefei_shannon_change <- calculate_shannon_change(group, "Hefei")
Hohhot_shannon_change <- calculate_shannon_change(group, "Hohhot")
Jinan_shannon_change <- calculate_shannon_change(group, "Jinan")
Kunming_shannon_change <- calculate_shannon_change(group, "Kunming")
Lanzhou_shannon_change <- calculate_shannon_change(group, "Lanzhou")
Nanjing_shannon_change <- calculate_shannon_change(group, "Nanjing")
Nanning_shannon_change <- calculate_shannon_change(group, "Nanning")
Shanghai_shannon_change <- calculate_shannon_change(group, "Shanghai")
Shenyang_shannon_change <- calculate_shannon_change(group, "Shenyang")
Shijiazhuang_shannon_change <- calculate_shannon_change(group, "Shijiazhuang")
Taiyuan_shannon_change <- calculate_shannon_change(group, "Taiyuan")
Tianjin_shannon_change <- calculate_shannon_change(group, "Tianjin")
Urumqi_shannon_change <- calculate_shannon_change(group, "Urumqi")
Wuhan_shannon_change <- calculate_shannon_change(group, "Wuhan")
Xian_shannon_change <- calculate_shannon_change(group, "Xian")
Xining_shannon_change <- calculate_shannon_change(group, "Xining")
yinchuan_shannon_change <- calculate_shannon_change(group, "Yinchuan")
Zhengzhou_shannon_change <- calculate_shannon_change(group, "Zhengzhou")

high_shannon_change<-rbind(Beijing_shannon_change,Changchun_shannon_change,Changsha_shannon_change,
                           Chengdu_shannon_change, Fuzhou_shannon_change,Guangzhou_shannon_change,Guiyang_shannon_change
                           ,Hangzhou_shannon_change,Harbin_shannon_change,Hefei_shannon_change,Hohhot_shannon_change,
                           Jinan_shannon_change,Kunming_shannon_change,Lanzhou_shannon_change,Nanjing_shannon_change,
                           Nanning_shannon_change,Shanghai_shannon_change,Shenyang_shannon_change,Shijiazhuang_shannon_change,
                           Taiyuan_shannon_change,Tianjin_shannon_change,Urumqi_shannon_change,
                           Wuhan_shannon_change,Xian_shannon_change,Xining_shannon_change,yinchuan_shannon_change,Zhengzhou_shannon_change)


colnames(high_shannon_change)[2]<-"group"
high_shannon_change$group<-as.factor(high_shannon_change$group)

high_shannon_change$group <- gsub("forest_cropland", "Cropland vs Forest", high_shannon_change$group)
high_shannon_change$group <- gsub("forest_industry", "Industry vs Forest", high_shannon_change$group)
high_shannon_change$group <- gsub("forest_urban", "Urban vs Forest", high_shannon_change$group)
colnames(high_shannon_change)[3]<-"shannon"

high_pH_shannon<-cbind(high_shannon_change,pH_change)
colnames(high_pH_shannon)[2]<-"group"
high_pH_shannon<-high_pH_shannon[,-c(5,6)]
high_pH_shannon$group1<-"high"

##计算bac

group<-bac_pH_shannon


results <- data.frame(Group = character(), 
                      Spearman_Correlation = numeric(), 
                      P_value = numeric(), 
                      stringsAsFactors = FALSE)

# 按 group 列分组计算
for (group_name in unique(group$group)) {
  # 提取当前组的数据
  group_data <- subset(group, group == group_name)
  
  # 计算 Spearman 相关系数和 p 值，并抑制警告
  cor_test <- suppressWarnings(
    cor.test(group_data$pH, group_data$shannon, method = "spearman")
  )
  
  # 将结果添加到 results 数据框
  results <- rbind(results, data.frame(Group = group_name,
                                       Spearman_Correlation = cor_test$estimate,
                                       P_value = cor_test$p.value))
}
results

format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 使用科学计数法
  } else {
    round(p, 3)
  }
}

cropland_P_value <-format_p(results$P_value[results$Group == "Cropland vs Forest"])
cropland_r_value <-round((results$Spearman_Correlation[results$Group == "Cropland vs Forest"]),3)
label_Cropland<-paste0("Cropland vs Forest: ρ = ", cropland_r_value, ", p = ", cropland_P_value)


industry_P_value <-format_p(results$P_value[results$Group == "Industry vs Forest"])
industry_r_value <-round((results$Spearman_Correlation[results$Group == "Industry vs Forest"]),3)
label_industry<-paste0("Industry vs Forest: ρ = ", industry_r_value, ", p = ", industry_P_value)

urban_P_value <-format_p(results$P_value[results$Group == "Urban vs Forest"])
urban_r_value <-round((results$Spearman_Correlation[results$Group == "Urban vs Forest"]),3)
label_urban<-paste0("Urban vs Forest: ρ = ", urban_r_value, ", p = ", urban_P_value)
ecosystem_labels <- paste0(
  label_Cropland, "\n",
  label_industry, "\n",
  label_urban
)





ecosystem_order<-c("Cropland vs Forest","Industry vs Forest","Urban vs Forest")
group <- group %>%
  mutate(R = case_when(
    group == "Urban vs Forest" ~ urban_r_value,
    group == "Industry vs Forest" ~ industry_r_value,
    group == "Cropland vs Forest" ~ cropland_r_value
  ))



cropland_P_value <- as.numeric(cropland_P_value)
industry_P_value <- as.numeric(industry_P_value)
urban_P_value <- as.numeric(urban_P_value)

group <- group %>%
  mutate(P = case_when(
    group == "Urban vs Forest" ~ urban_P_value,
    group == "Industry vs Forest" ~ industry_P_value,
    group == "Cropland vs Forest" ~ cropland_P_value
  ))

group <- group %>%
  mutate(linetype_group = case_when(
    P < 0.05 ~ "solid",
    P >= 0.05 ~ "dashed"))

ecosystem_labels 
label_data <- data.frame(
  group = c("Cropland vs Forest", "Industry vs Forest", "Urban vs Forest"),
  label = c(
    "ρ = 0.485, p = 0.011",
    "ρ = 0.676, p = 1.60e-04",
    "ρ = 0.562, p = 2.70e-03"
  ),
  x = 2,                  # x 坐标
  y = c(2.2, 2.0,1.8),  # y 坐标手动错开
  color = c("#00FFFF","#5e9ceb","#ff7f27")  # 对应 group 颜色
)

group_bac<-group

#############计算arg

group<-arg_pH_shannon


results <- data.frame(Group = character(), 
                      Spearman_Correlation = numeric(), 
                      P_value = numeric(), 
                      stringsAsFactors = FALSE)

# 按 group 列分组计算
for (group_name in unique(group$group)) {
  # 提取当前组的数据
  group_data <- subset(group, group == group_name)
  
  # 计算 Spearman 相关系数和 p 值，并抑制警告
  cor_test <- suppressWarnings(
    cor.test(group_data$pH, group_data$shannon, method = "spearman")
  )
  
  # 将结果添加到 results 数据框
  results <- rbind(results, data.frame(Group = group_name,
                                       Spearman_Correlation = cor_test$estimate,
                                       P_value = cor_test$p.value))
}
results

format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 使用科学计数法
  } else {
    round(p, 3)
  }
}

cropland_P_value <-format_p(results$P_value[results$Group == "Cropland vs Forest"])
cropland_r_value <-round((results$Spearman_Correlation[results$Group == "Cropland vs Forest"]),3)
label_Cropland<-paste0("Cropland vs Forest: ρ = ", cropland_r_value, ", p = ", cropland_P_value)


industry_P_value <-format_p(results$P_value[results$Group == "Industry vs Forest"])
industry_r_value <-round((results$Spearman_Correlation[results$Group == "Industry vs Forest"]),3)
label_industry<-paste0("Industry vs Forest: ρ = ", industry_r_value, ", p = ", industry_P_value)

urban_P_value <-format_p(results$P_value[results$Group == "Urban vs Forest"])
urban_r_value <-round((results$Spearman_Correlation[results$Group == "Urban vs Forest"]),3)
label_urban<-paste0("Urban vs Forest: ρ = ", urban_r_value, ", p = ", urban_P_value)
ecosystem_labels <- paste0(
  label_Cropland, "\n",
  label_industry, "\n",
  label_urban
)





ecosystem_order<-c("Cropland vs Forest","Industry vs Forest","Urban vs Forest")
group <- group %>%
  mutate(R = case_when(
    group == "Urban vs Forest" ~ urban_r_value,
    group == "Industry vs Forest" ~ industry_r_value,
    group == "Cropland vs Forest" ~ cropland_r_value
  ))



cropland_P_value <- as.numeric(cropland_P_value)
industry_P_value <- as.numeric(industry_P_value)
urban_P_value <- as.numeric(urban_P_value)

group <- group %>%
  mutate(P = case_when(
    group == "Urban vs Forest" ~ urban_P_value,
    group == "Industry vs Forest" ~ industry_P_value,
    group == "Cropland vs Forest" ~ cropland_P_value
  ))

group <- group %>%
  mutate(linetype_group = case_when(
    P < 0.05 ~ "solid",
    P >= 0.05 ~ "dashed"))

ecosystem_labels 
label_data <- data.frame(
  group = c("Cropland vs Forest", "Industry vs Forest", "Urban vs Forest"),
  label = c(
    "ρ = 0.39, p = 0.045",
    "ρ = 0.579, p = 1.86e-03",
    "ρ = 0.392, p =  0.044"
  ),
  x = 2,                  # x 坐标
  y = c(2.2, 2.0,1.8),  # y 坐标手动错开
  color = c("#00FFFF","#5e9ceb","#ff7f27")  # 对应 group 颜色
)

group_arg<-group

##计算high risk
group<-high_pH_shannon


results <- data.frame(Group = character(), 
                      Spearman_Correlation = numeric(), 
                      P_value = numeric(), 
                      stringsAsFactors = FALSE)

# 按 group 列分组计算
for (group_name in unique(group$group)) {
  # 提取当前组的数据
  group_data <- subset(group, group == group_name)
  
  # 计算 Spearman 相关系数和 p 值，并抑制警告
  cor_test <- suppressWarnings(
    cor.test(group_data$pH, group_data$shannon, method = "spearman")
  )
  
  # 将结果添加到 results 数据框
  results <- rbind(results, data.frame(Group = group_name,
                                       Spearman_Correlation = cor_test$estimate,
                                       P_value = cor_test$p.value))
}
results

format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 使用科学计数法
  } else {
    round(p, 3)
  }
}

cropland_P_value <-format_p(results$P_value[results$Group == "Cropland vs Forest"])
cropland_r_value <-round((results$Spearman_Correlation[results$Group == "Cropland vs Forest"]),3)
label_Cropland<-paste0("Cropland vs Forest: ρ = ", cropland_r_value, ", p = ", cropland_P_value)


industry_P_value <-format_p(results$P_value[results$Group == "Industry vs Forest"])
industry_r_value <-round((results$Spearman_Correlation[results$Group == "Industry vs Forest"]),3)
label_industry<-paste0("Industry vs Forest: ρ = ", industry_r_value, ", p = ", industry_P_value)

urban_P_value <-format_p(results$P_value[results$Group == "Urban vs Forest"])
urban_r_value <-round((results$Spearman_Correlation[results$Group == "Urban vs Forest"]),3)
label_urban<-paste0("Urban vs Forest: ρ = ", urban_r_value, ", p = ", urban_P_value)
ecosystem_labels <- paste0(
  label_Cropland, "\n",
  label_industry, "\n",
  label_urban
)





ecosystem_order<-c("Cropland vs Forest","Industry vs Forest","Urban vs Forest")
group <- group %>%
  mutate(R = case_when(
    group == "Urban vs Forest" ~ urban_r_value,
    group == "Industry vs Forest" ~ industry_r_value,
    group == "Cropland vs Forest" ~ cropland_r_value
  ))



cropland_P_value <- as.numeric(cropland_P_value)
industry_P_value <- as.numeric(industry_P_value)
urban_P_value <- as.numeric(urban_P_value)

group <- group %>%
  mutate(P = case_when(
    group == "Urban vs Forest" ~ urban_P_value,
    group == "Industry vs Forest" ~ industry_P_value,
    group == "Cropland vs Forest" ~ cropland_P_value
  ))

group <- group %>%
  mutate(linetype_group = case_when(
    P < 0.05 ~ "solid",
    P >= 0.05 ~ "dashed"))

ecosystem_labels 


group_high<-group
group_high$group1<-"high"
group_arg$group1<-"arg"

group<-rbind(group_bac,group_arg,group_high)

library(dplyr)

group_cropland <- group %>% 
  filter(group == "Cropland vs Forest")

group_industry <- group %>% 
  filter(group == "Industry vs Forest")

group_urban <- group %>% 
  filter(group == "Urban vs Forest")

#####cropland difference between forest

group_order<-c("bacteria","arg","high")

label_data_cropland <- data.frame(
  group1 = c("bacteria", "arg", "high"),
  label = c(
    "Bacteria: ρ =  0.485, P = 0.011",
    "Resistome: ρ = 0.390, P = 0.045",
    "High risk resistome: ρ = -0.384, P =  0.049"
  ),
  x = -0.5,                  # x 坐标
  y = c(2.4, 2.2,2.0),  # y 坐标手动错开
  color = c("#00FFFF","#5e9ceb","#ff7f27")  # 对应 group 颜色
)


p <- ggplot(group_cropland, aes(x = pH, y = shannon, color = group1)) +
  geom_jitter(size = 3, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = " ",
    x = "Soil pH differnece\nbetween cropland and forest", y = "H' difference\nbetween cropland and forest"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # 主网格线
    panel.grid.minor = element_line(color = "grey90", size = 0.3)   # 次网格线
  ) +
  theme(legend.position = "none")+
  # 对每个 group 绘制趋势线
  geom_smooth(
    aes(x = pH, y = shannon, color = group1, linetype = linetype_group, size = abs(R)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3,
    show.legend = FALSE
  )+
  #scale_size_continuous(
  #range = c(1, 4),
  #guide = "none"
  #) +
  scale_color_manual(
    name = " ",
    limits = group_order,
    values = c("#FA9A38","#6B30C5","#097369")
  ) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(color = "black", size = 16, face = "bold"),
    axis.text.x = element_text(size = 18, color = "black", vjust = 0.5, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm")
  ) +
  geom_text(
    data = label_data_cropland,
    aes(x = x, y = y, label = label, color = group1),
    hjust = 0,
    size = 5,
    show.legend = FALSE
  ) +
  guides(color = guide_legend(
    override.aes = list(alpha = 1)  # 图例中的点/线不透明
  ))



p

ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/pH_difference_shannon_difference_cropland_combine.pdf",p,height=6,width=6, device = cairo_pdf
)


##industry

group_order<-c("bacteria","arg","high")

label_data_industry <- data.frame(
  group1 = c("bacteria", "arg", "high"),
  label = c(
    "Bacteria: ρ =  0.676, P = 0.0001",
    "Resistome: ρ = 0.579, P = 0.002",
    "High risk resistome: ρ = -0.382, P =  0.050"
  ),
  x = 0.1,                  # x 坐标
  y = c(1.2, 1.0,0.8),  # y 坐标手动错开
  color = c("#00FFFF","#5e9ceb","#ff7f27")  # 对应 group 颜色
)


p <- ggplot(group_industry, aes(x = pH, y = shannon, color = group1)) +
  geom_jitter(size = 3, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = " ",
    x = "Soil pH difference\nbetween industry and forest",
    y = "H' difference\nbetween industry and forest"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_line(color = "grey90", size = 0.3)
  ) +
  theme(legend.position = "none") +
  
  # ✅ 对每个 group 绘制趋势线
  geom_smooth(
    aes(linetype = linetype_group, size = abs(R)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  
  # scale_size_continuous(
  #   range = c(1, 4),
  #   guide = "none"
  # ) +
  scale_color_manual(
    name = " ",
    limits = group_order,
    values = c("#FA9A38", "#6B30C5", "#097369")
  ) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(color = "black", size = 16, face = "bold"),
    axis.text.x = element_text(size = 18, color = "black", vjust = 0.5, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm")
  ) +
  geom_text(
    data = label_data_industry,
    aes(x = x, y = y, label = label, color = group1),
    hjust = 0,
    size = 5,
    show.legend = FALSE
  ) +
  guides(color = guide_legend(
    override.aes = list(alpha = 1)
  ))



p

ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/pH_difference_shannon_difference_industry_combine.pdf",p,height=6,width=6,  device = cairo_pdf
)


#urban
group_order<-c("bacteria","arg","high")

label_data_urban <- data.frame(
  group1 = c("bacteria", "arg", "high"),
  label = c(
    "Bacteria: ρ =  0.562, P = 0.003",
    "Resistome: ρ = 0.392, P = 0.044",
    "High risk resistome: ρ = -0.320, P =  0.104"
  ),
  x = 0.1,                  # x 坐标
  y = c(2.0, 1.8,1.6),  # y 坐标手动错开
  color = c("#00FFFF","#5e9ceb","#ff7f27")  # 对应 group 颜色
)


p <- ggplot(group_urban, aes(x = pH, y = shannon, color = group1)) +
  geom_jitter(size = 3, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = " ",
    x = "Soil pH differnece\nbetween urban and forest", y = "H' difference\nbetween urban and forest"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # 主网格线
    panel.grid.minor = element_line(color = "grey90", size = 0.3)   # 次网格线
  ) +
  theme(legend.position = "none")+
  # 对每个 group 绘制趋势线
  geom_smooth(
    aes(x = pH, y = shannon, color = group1, linetype = linetype_group, size = abs(R)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3,
    show.legend = FALSE
  )+
  #scale_size_continuous(
  #range = c(1, 4),
  #guide = "none"
  #) +
  scale_color_manual(
    name = " ",
    limits = group_order,
    values = c("#FA9A38","#6B30C5","#097369")
  ) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(color = "black", size = 16, face = "bold"),
    axis.text.x = element_text(size = 18, color = "black", vjust = 0.5, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm")
  ) +
  geom_text(
    data = label_data_urban,
    aes(x = x, y = y, label = label, color = group1),
    hjust = 0,
    size = 5,
    show.legend = FALSE
  ) +
  guides(color = guide_legend(
    override.aes = list(alpha = 1)  # 图例中的点/线不透明
  ))



p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/pH_difference_shannon_difference_urban_combine.pdf",p,height=6,width=6,  device = cairo_pdf)

#Fig.3j
##MGE copy 和alpha diversity的关系
rm(list=ls())
library(ggplot2)
library(agricolae)
arg<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-4544]
group$shannon <- vegan::diversity(arg, index = "shannon")


mge<-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
mge<-t(mge)
mge <- mge[order(rownames(mge)), ]
identical(rownames(mge), group$sample.id)

group$mge_richness<-rowSums(mge > 0)
group$mge_abundance<-rowSums(mge)

merge<-group


results <- data.frame(Group = character(), 
                      Spearman_Correlation = numeric(), 
                      P_value = numeric(), 
                      stringsAsFactors = FALSE)

# 按 group 列分组计算
for (group_name in unique(merge$group)) {
  # 提取当前组的数据
  group_data <- subset(merge, group == group_name)
  
  # 计算 Spearman 相关系数和 p 值，并抑制警告
  cor_test <- suppressWarnings(
    cor.test(group_data$shannon, group_data$mge_abundance, method = "spearman")
  )
  
  # 将结果添加到 results 数据框
  results <- rbind(results, data.frame(Group = group_name,
                                       Spearman_Correlation = cor_test$estimate,
                                       P_value = cor_test$p.value))
}
results

format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 使用科学计数法
  } else {
    round(p, 3)
  }
}
cropland_P_value <-format_p(results$P_value[results$Group == "Cropland"])
cropland_r_value <-round((results$Spearman_Correlation[results$Group == "Cropland"]),3)
label_Cropland<-paste0("Cropland: ρ = ", cropland_r_value, ", P= ", cropland_P_value)

forest_P_value <-format_p(results$P_value[results$Group == "Forest"])
forest_r_value <-round((results$Spearman_Correlation[results$Group == "Forest"]),3)
label_forest<-paste0("Forest: ρ = ", forest_r_value, ", P = ", forest_P_value)

industry_P_value <-format_p(results$P_value[results$Group == "Industry"])
industry_r_value <-round((results$Spearman_Correlation[results$Group == "Industry"]),3)
label_industry<-paste0("Industry: ρ = ", industry_r_value, ", P = ", industry_P_value)

urban_P_value <-format_p(results$P_value[results$Group == "Urban"])
urban_r_value <-round((results$Spearman_Correlation[results$Group == "Urban"]),3)
label_urban<-paste0("Urban: ρ = ", urban_r_value, ", P = ", urban_P_value)
colnames(results)

results <- results %>%
  mutate(
    label_text = paste0("ρ = ", round(Spearman_Correlation, 2),
                        ", P = ", sapply(P_value, format_p))
  )
results <- results %>%
  mutate(label_text = paste0(Group, ": ", label_text))

ecosystem_order<-c("Forest","Cropland","Industry","Urban")


p_values <- data.frame(
  group = c("Cropland", "Industry", "Forest", "Urban"),
  p_value = c(1.551809e-02, 6.115481e-01, 3.393330e-02, 4.296132e-05)
)

library(dplyr)
p_values <- p_values %>%
  mutate(linetype_group = ifelse(p_value > 0.05, "dashed", "solid"))
merge <- merge %>%
  left_join(p_values[, c("group", "linetype_group")], by = "group")

merge <- merge %>%
  mutate(R_value = case_when(
    group == "Urban" ~ urban_r_value,
    group == "Industry" ~ industry_r_value,
    group == "Forest" ~ forest_r_value,
    group == "Cropland" ~ cropland_r_value
  ))

colnames(results)[1]<-"group"

merge <- merge %>%
  left_join(results %>% select(group, P_value), by = "group")

ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")
results$group <- factor(results$group, levels = ecosystem_order)
results <- results[order(results$group), ]

results$x <- c(0.2, 0.2, 0.2, 0.2)
results$y<-c(6.3, 6.2,6.1, 6.0)

p <- ggplot(merge, aes(x = mge_abundance, y = shannon, color = group)) +
  geom_jitter(size = 3, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = " ",
    x = "Copies of MGEs per cell", y = "Resistome H'"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # 主网格线
    panel.grid.minor = element_line(color = "grey90", size = 0.3)   # 次网格线
  ) +
  # 对每个 group 绘制趋势线
  geom_smooth(
    aes(x = mge_abundance, y = shannon, color = group, linetype = linetype_group, size = abs(R_value)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3,
    show.legend = FALSE
  )+
  #scale_size_continuous(
  #range = c(1, 4),
  #guide = "none"
  #) +
  scale_color_manual(
    name = "Ecosystem",
    limits = ecosystem_order,
    guide="none",
    values = c("#6cb509","#00FFFF","#5e9ceb","#ff7f27")
  ) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  theme(
    plot.title = element_text(size = 18, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "right",
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 18, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
    axis.title.x = element_text(size = 18, color = "black", face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(color = "black", size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", vjust = 0.5, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm")
  ) +
  geom_text(
    data = results,
    aes(x = x, y = y, label = label_text, color = group),
    inherit.aes = FALSE,
    hjust = 0,  # 左对齐
    size = 5,
    show.legend = FALSE 
  ) 

p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/arg_mge_copy_alpha.pdf",p,width = 6,height = 6,device = cairo_pdf)

#Fig.3k
rm(list=ls())
library(ggplot2)
library(agricolae)
arg <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
arg<-arg[,-1]
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-74]
group$shannon <- vegan::diversity(arg, index = "shannon")




mge<-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
mge<-t(mge)
mge <- mge[order(rownames(mge)), ]
identical(rownames(mge), group$sample.id)

group$mge_richness<-rowSums(mge > 0)
group$mge_abundance<-rowSums(mge)

merge<-group


results <- data.frame(Group = character(), 
                      Spearman_Correlation = numeric(), 
                      P_value = numeric(), 
                      stringsAsFactors = FALSE)

# 按 group 列分组计算
for (group_name in unique(merge$group)) {
  # 提取当前组的数据
  group_data <- subset(merge, group == group_name)
  
  # 计算 Spearman 相关系数和 p 值，并抑制警告
  cor_test <- suppressWarnings(
    cor.test(group_data$shannon, group_data$mge_abundance, method = "spearman")
  )
  
  # 将结果添加到 results 数据框
  results <- rbind(results, data.frame(Group = group_name,
                                       Spearman_Correlation = cor_test$estimate,
                                       P_value = cor_test$p.value))
}
results

format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 使用科学计数法
  } else {
    round(p, 3)
  }
}
cropland_P_value <-format_p(results$P_value[results$Group == "Cropland"])
cropland_r_value <-round((results$Spearman_Correlation[results$Group == "Cropland"]),3)
label_Cropland<-paste0("Cropland: ρ = ", cropland_r_value, ", P= ", cropland_P_value)

forest_P_value <-format_p(results$P_value[results$Group == "Forest"])
forest_r_value <-round((results$Spearman_Correlation[results$Group == "Forest"]),3)
label_forest<-paste0("Forest: ρ = ", forest_r_value, ", P = ", forest_P_value)

industry_P_value <-format_p(results$P_value[results$Group == "Industry"])
industry_r_value <-round((results$Spearman_Correlation[results$Group == "Industry"]),3)
label_industry<-paste0("Industry: ρ = ", industry_r_value, ", P = ", industry_P_value)

urban_P_value <-format_p(results$P_value[results$Group == "Urban"])
urban_r_value <-round((results$Spearman_Correlation[results$Group == "Urban"]),3)
label_urban<-paste0("Urban: ρ = ", urban_r_value, ", P = ", urban_P_value)
colnames(results)

results <- results %>%
  mutate(
    label_text = paste0("ρ = ", round(Spearman_Correlation, 2),
                        ", P = ", sapply(P_value, format_p))
  )
results <- results %>%
  mutate(label_text = paste0(Group, ": ", label_text))

ecosystem_order<-c("Forest","Cropland","Industry","Urban")


p_values <- data.frame(
  group = c("Cropland", "Industry", "Forest", "Urban"),
  p_value = c(2.094142e-01, 6.419180e-01, 5.045170e-01, 3.747462e-06)
)

library(dplyr)
p_values <- p_values %>%
  mutate(linetype_group = ifelse(p_value > 0.05, "dashed", "solid"))
merge <- merge %>%
  left_join(p_values[, c("group", "linetype_group")], by = "group")

merge <- merge %>%
  mutate(R_value = case_when(
    group == "Urban" ~ urban_r_value,
    group == "Industry" ~ industry_r_value,
    group == "Forest" ~ forest_r_value,
    group == "Cropland" ~ cropland_r_value
  ))

colnames(results)[1]<-"group"

merge <- merge %>%
  left_join(results %>% select(group, P_value), by = "group")

ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")
results$group <- factor(results$group, levels = ecosystem_order)
results <- results[order(results$group), ]

results$x <- c(0.4, 0.4, 0.4, 0.4)
results$y<-c(1.8, 1.72,1.64, 1.56)

p <- ggplot(merge, aes(x = mge_abundance, y = shannon, color = group)) +
  geom_jitter(size = 3, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = " ",
    x = "Copies of MGEs per cell", y = "High risk resistome H'"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # 主网格线
    panel.grid.minor = element_line(color = "grey90", size = 0.3)   # 次网格线
  ) +
  # 对每个 group 绘制趋势线
  geom_smooth(
    aes(x = mge_abundance, y = shannon, color = group, linetype = linetype_group, size = abs(R_value)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3,
    show.legend = FALSE
  )+
  #scale_size_continuous(
  #range = c(1, 4),
  #guide = "none"
  #) +
  scale_color_manual(
    name = "Ecosystem",
    limits = ecosystem_order,
    guide="none",
    values = c("#6cb509","#00FFFF","#5e9ceb","#ff7f27")
  ) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  theme(
    plot.title = element_text(size = 18, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "right",
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 18, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
    axis.title.x = element_text(size = 18, color = "black", face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(color = "black", size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", vjust = 0.5, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm")
  ) +
  geom_text(
    data = results,
    aes(x = x, y = y, label = label_text, color = group),
    inherit.aes = FALSE,
    hjust = 0,  # 左对齐
    size = 5,
    show.legend = FALSE 
  ) 

p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/alpha_diversity/arg_mge_copy_alpha_high.pdf",p,width=6,height=6,device = cairo_pdf)






