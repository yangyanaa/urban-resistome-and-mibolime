#fig4 effect size on beta diversity
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

arg0<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv",sep=",",header=T,row.names=1)
arg0<-arg0[,-1]
arg0 <- as.data.frame(t(arg0))
#arg0_Flattening = as.data.frame(t(rrarefy(t(arg0), min(colSums(arg0)))))


group0<-read.csv("/Users/yangy/Documents/E//city/metagenome/group.csv",sep=",",header=T,row.names=1)
group0 <- na.omit(group0)
group <- group0[order(rownames(group0)), ]
arg0$sample.id <- rownames(arg0)
group0$sample.id <- rownames(group0)

colnames(arg0)[(ncol(arg0)-9):ncol(arg0)]

arg <- merge(arg0, group0, by.x="sample.id", by.y="sample.id")
colnames(arg)[(ncol(arg)-19):ncol(arg)]

dis <- vegdist(decostand(arg[,2:74],"hellinger"), method = 'bray')
Group <- group$group
Location<-group$Location

mod <- betadisper(dis,Group, type = 'centroid')

distances_to_centroid <- mod$distances

anova_data <- data.frame(
  Distance = distances_to_centroid,
  Location = Location,
  Group = Group
)


kruskal.test(anova_data$Distance,anova_data$Group)
library(broom)
tidy(kruskal.test(anova_data$Distance,anova_data$Group))

ST <- kruskal(anova_data$Distance,anova_data$Group, p.adj = "fdr")
statis <- as.data.frame(rbind(ST$statistics[,1:2]))
rp.value <- data.frame(
  label=sprintf("P = %.3f (χ² = %.3g)", 
                statis$p.chisq,statis$Chisq),
  x = "Tamarix chinensis",
  y=0.55)
statis1 <- as.data.frame(rbind(ST$groups))
statis2 <- as.data.frame(rbind(ST$means))
p1 <- ggplot(anova_data, aes( y=  Distance, x = factor(Group))) +
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
  Group     = statis2$id,
  y     = as.numeric(statis2$ymax))


data<-anova_data
data$Group[data$Group == "Forest"] <- "CK"
colnames(data)[3]<-"TR"
library(tidyverse) 
library(metafor) 

# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(Distance), # 计算生物量均值
            SD = sd(Distance), 
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

otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E//city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)

otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu_sub <- otu_sub[, colSums(otu_sub[,1:6927]) > 0]

colnames(otu_sub)[(ncol(otu_sub)-9):ncol(otu_sub)]
dis <- vegdist(decostand(otu_sub[,1:6641],"hellinger"), method = 'bray')


Location <- group$Location
Group <- group$group
mod <- betadisper(dis,Group, type = 'centroid')

distances_to_centroid <- mod$distances
anova_data <- data.frame(
  Distance = distances_to_centroid,
  Location = Location,
  Group = Group
)

kruskal.test(anova_data$Distance,anova_data$Group)
library(broom)
tidy(kruskal.test(anova_data$Distance,anova_data$Group))

ST <- kruskal(anova_data$Distance,anova_data$Group, p.adj = "fdr")
statis <- as.data.frame(rbind(ST$statistics[,1:2]))
rp.value <- data.frame(
  label=sprintf("P = %.3f (χ² = %.3g)", 
                statis$p.chisq,statis$Chisq),
  x = "Tamarix chinensis",
  y=0.55)
statis1 <- as.data.frame(rbind(ST$groups))
statis2 <- as.data.frame(rbind(ST$means))
p1 <- ggplot(anova_data, aes( y=  Distance, x = factor(Group))) +
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
  Group     = statis2$id,
  y     = as.numeric(statis2$ymax))

multcomp.data$Group <- factor(multcomp.data$Group, levels = c("Forest", "Cropland", "Industry", "Urban"))
library(gghalves)

anova_sub <- subset(anova_data, Group %in% c("Forest", "Urban"))
data<-anova_data
data$Group[data$Group == "Forest"] <- "CK"
colnames(data)[3]<-"TR"
library(tidyverse) 
library(metafor) 

# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(Distance), # 计算生物量均值
            SD = sd(Distance), 
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
# 为显著与非显著设置填充色
plot_data$fill_color <- ifelse(
  plot_data$significant == "Significant",
  custom_colors[plot_data$Subgroup],
  "white"
)
plot_data$point_shape <- ifelse(plot_data$significant == "Significant", 16, 21)
plot_data_bacteria<-plot_data
plot_data_bacteria$group<-"bacteria"
# arg
arg_all0<- read.csv("/Users/yangy/Documents/E//city/metagenome/arg/results/normalized_percell.csv",sep=",",header=T,row.names=1)

arg_all0 <- as.data.frame(t(arg_all0))

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E//city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg_all0$sample.id <- rownames(arg_all0)

colnames(arg_all0)[(ncol(arg_all0)-9):ncol(arg_all0)]
dis <- vegdist(decostand(arg_all0[,1:4543],"hellinger"), method = 'bray')


Location <- group$Location
Group <- group$group
mod <- betadisper(dis,Group, type = 'centroid')

distances_to_centroid <- mod$distances
anova_data <- data.frame(
  Distance = distances_to_centroid,
  Location = Location,
  Group = Group
)

kruskal.test(anova_data$Distance,anova_data$Group)
library(broom)
tidy(kruskal.test(anova_data$Distance,anova_data$Group))

ST <- kruskal(anova_data$Distance,anova_data$Group, p.adj = "fdr")
statis <- as.data.frame(rbind(ST$statistics[,1:2]))
rp.value <- data.frame(
  label=sprintf("P = %.3f (χ² = %.3g)", 
                statis$p.chisq,statis$Chisq),
  x = "Tamarix chinensis",
  y=0.55)
statis1 <- as.data.frame(rbind(ST$groups))
statis2 <- as.data.frame(rbind(ST$means))
p1 <- ggplot(anova_data, aes( y=  Distance, x = factor(Group))) +
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
  Group     = statis2$id,
  y     = as.numeric(statis2$ymax))

multcomp.data$Group <- factor(multcomp.data$Group, levels = c("Forest", "Cropland", "Industry", "Urban"))


data<-anova_data
data$Group[data$Group == "Forest"] <- "CK"
colnames(data)[3]<-"TR"
library(tidyverse) 
library(metafor) 

# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(Distance), # 计算生物量均值
            SD = sd(Distance), 
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
plot_data$group<-"ARG"
plot_data_arg<-plot_data


# vfg
vfg_all0<- read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)

vfg_all0 <- as.data.frame(t(vfg_all0))

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E//city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
vfg_all0$sample.id <- rownames(vfg_all0)

colnames(vfg_all0)[(ncol(vfg_all0)-9):ncol(vfg_all0)]
dis <- vegdist(decostand(vfg_all0[,1:5409],"hellinger"), method = 'bray')


Location <- group$Location
Group <- group$group
mod <- betadisper(dis,Group, type = 'centroid')

distances_to_centroid <- mod$distances
anova_data <- data.frame(
  Distance = distances_to_centroid,
  Location = Location,
  Group = Group
)

kruskal.test(anova_data$Distance,anova_data$Group)
library(broom)
tidy(kruskal.test(anova_data$Distance,anova_data$Group))

ST <- kruskal(anova_data$Distance,anova_data$Group, p.adj = "fdr")
statis <- as.data.frame(rbind(ST$statistics[,1:2]))
rp.value <- data.frame(
  label=sprintf("P = %.3f (χ² = %.3g)", 
                statis$p.chisq,statis$Chisq),
  x = "Tamarix chinensis",
  y=0.55)
statis1 <- as.data.frame(rbind(ST$groups))
statis2 <- as.data.frame(rbind(ST$means))
p1 <- ggplot(anova_data, aes( y=  Distance, x = factor(Group))) +
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
  Group     = statis2$id,
  y     = as.numeric(statis2$ymax))

multcomp.data$Group <- factor(multcomp.data$Group, levels = c("Forest", "Cropland", "Industry", "Urban"))


data<-anova_data
data$Group[data$Group == "Forest"] <- "CK"
colnames(data)[3]<-"TR"
library(tidyverse) 
library(metafor) 

# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(Distance), # 计算生物量均值
            SD = sd(Distance), 
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
plot_data$group<-"vfg"
plot_data_vfg<-plot_data



#MGE
# vfg
mge_all0<- read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)

mge_all0 <- as.data.frame(t(mge_all0))

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E//city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
mge_all0$sample.id <- rownames(mge_all0)

colnames(mge_all0)[(ncol(mge_all0)-9):ncol(mge_all0)]
dis <- vegdist(decostand(mge_all0[,1:1431],"hellinger"), method = 'bray')


Location <- group$Location
Group <- group$group
mod <- betadisper(dis,Group, type = 'centroid')

distances_to_centroid <- mod$distances
anova_data <- data.frame(
  Distance = distances_to_centroid,
  Location = Location,
  Group = Group
)

kruskal.test(anova_data$Distance,anova_data$Group)
library(broom)
tidy(kruskal.test(anova_data$Distance,anova_data$Group))

ST <- kruskal(anova_data$Distance,anova_data$Group, p.adj = "fdr")
statis <- as.data.frame(rbind(ST$statistics[,1:2]))
rp.value <- data.frame(
  label=sprintf("P = %.3f (χ² = %.3g)", 
                statis$p.chisq,statis$Chisq),
  x = "Tamarix chinensis",
  y=0.55)
statis1 <- as.data.frame(rbind(ST$groups))
statis2 <- as.data.frame(rbind(ST$means))
p1 <- ggplot(anova_data, aes( y=  Distance, x = factor(Group))) +
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
  Group     = statis2$id,
  y     = as.numeric(statis2$ymax))

multcomp.data$Group <- factor(multcomp.data$Group, levels = c("Forest", "Cropland", "Industry", "Urban"))


data<-anova_data
data$Group[data$Group == "Forest"] <- "CK"
colnames(data)[3]<-"TR"
library(tidyverse) 
library(metafor) 

# 计算各处理组和对照组的描述统计量
summary_data <- data %>% group_by(TR) %>% # 按处理组分组
  summarize(Mean = mean(Distance), # 计算生物量均值
            SD = sd(Distance), 
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
plot_data$group<-"mge"
plot_data_mge<-plot_data






plot<-rbind(plot_data_bacteria,plot_data_arg,plot_data_high,plot_data_vfg,plot_data_mge)

##urban effect
plot_urban <- plot %>%
  dplyr::filter(Subgroup == "Urban vs Forest")
unique(plot_urban$group)

custom_colors <- c(
  "bacteria" = "#FA9A38",       
  "ARG" = "#6B30C5", 
  "high" = "#097369",
  "vfg"="#A13B46",
  "mge"="#777777"
)
plot_urban$group <- factor(plot_urban$group,
                           levels = rev(c("bacteria", "ARG", "high","vfg","mge")))
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
  
  geom_text(aes(label = sig1), x = -0.2, size = 6, hjust = 0,nudge_y = -0.15 ) +
  
  scale_color_manual(values = custom_colors, name = "") +
  scale_fill_identity() +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Non-significant" = 21)) +
  scale_linetype_manual(name = "Error bars",
                        values = c("95% CI" = "solid", "±SE" = "solid"),
                        guide = guide_legend(override.aes = list(size = c(1.5, 2.5),
                                                                 alpha = c(0.4,0.7),
                                                                 color = "black"))) +
  
  scale_y_discrete(labels = c("bacteria" = "Bacterial community",
                              "ARG" = "Resistome",
                              "high" = "High risk resistome",
                              "vfg"="Virulome",
                              "mge"="Mobilome")) +
  
  scale_x_continuous(
    limits = c(-0.5, 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 0.6) +
  
  labs(x = ' ', y = '', title = "Urban effect on beta diversity") +
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

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/effect_size_β_dispersation_urban_effectt加上vfg和mge版本.pdf",p,width=8,height=6)


#industry
plot_industry <- plot %>%
  dplyr::filter(Subgroup == "Industry vs Forest")
unique(plot_industry$group)

custom_colors <- c(
  "bacteria" = "#FA9A38",       
  "ARG" = "#6B30C5", 
  "high" = "#097369",
  "vfg"="#A13B46",
  "mge"="#777777"
)
plot_industry$group <- factor(plot_industry$group,
                              levels = rev(c("bacteria", "ARG", "high","vfg","mge")))


format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 科学计数法
  } else {
    sprintf("%.3f", round(p, 3))          # 三位小数
  }
}

plot_industry$sig1 <- sprintf(
  "estimate = %.3f, P = %s",
  round(plot_industry$estimate, 3),
  sapply(plot_industry$pval, format_p)
)

p <- ggplot(plot_industry, aes(x = estimate, y = group)) +
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
  
  # 点估计
  geom_point(aes(
    fill = fill_color,
    color = group,
    shape = significant
  ),
  size = 6, stroke = 1) +
  
  geom_text(aes(label = sig1), x = -0.2, size = 6, hjust = 0,nudge_y = -0.15 ) +
  
  scale_color_manual(values = custom_colors, name = "") +
  scale_fill_identity() +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Non-significant" = 21)) +
  # ✅ 自定义 y 轴显示名称
  scale_y_discrete(labels = c("bacteria" = "Bacterial community",
                              "ARG" = "Resistome",
                              "high" = "High risk resistome",
                              "vfg"="Virulome",
                              "mge"="Mobilome")) +
  
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 0.6) +
  
  labs(x = ' ', y = '', title = "Industry effect on beta diversity") +
  theme_few() +
  # 线型图例
  scale_linetype_manual(
    name = "Error bars",
    values = c("95% CI" = "solid", "±SE" = "solid"),
    guide = guide_legend(
      override.aes = list(
        size = c(1.5, 2.5),
        alpha = c(0.4, 0.7),
        color = "black"
      )
    )
  ) +
  scale_x_continuous(
    limits = c(-0.5, 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",  # 隐藏图例
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 18, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    #panel.grid.major = element_line(color = "grey80", size = 0.1),
    axis.title.x = element_blank(),  # 隐藏x轴标题
    #panel.grid.minor = element_line(color = "grey90", size = 0.1)
  )

p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/effect_size_β_dispersation_industry_effect加上vfg和mge的版本.pdf",p,width=8,height=6)

#cropland effect on beta diversity

plot_cropland <- plot %>%
  dplyr::filter(Subgroup == "Cropland vs Forest")
unique(plot_cropland$group)


plot_cropland$group <- factor(plot_cropland$group,
                              levels = rev(c("bacteria", "ARG", "high","vfg","mge")))

format_p <- function(p) {
  if (p < 0.01) {
    formatC(p, format = "e", digits = 2)  # 科学计数法
  } else {
    sprintf("%.3f", round(p, 3))          # 三位小数
  }
}

plot_cropland$sig1 <- sprintf(
  "estimate = %.3f, P = %s",
  round(plot_cropland$estimate, 3),
  sapply(plot_cropland$pval, format_p)
)
p <- ggplot(plot_cropland, aes(x = estimate, y = group)) +
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
  
  # 点估计
  geom_point(aes(
    fill = fill_color,
    color = group,
    shape = significant
  ),
  size = 6, stroke = 1) +
  
  geom_text(aes(label = sig1), x = -0.2, size = 6, hjust = 0,nudge_y = -0.15 ) +
  
  scale_color_manual(values = custom_colors, name = "") +
  scale_fill_identity() +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Non-significant" = 21)) +
  
  # ✅ 自定义 y 轴显示名称
  scale_y_discrete(labels = c("bacteria" = "Bacterial community",
                              "ARG" = "Resistome",
                              "high" = "High risk resistome",
                              "vfg"="Virulome",
                              "mge"="Mobilome")) +
  
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 0.6) +
  
  labs(x = ' ', y = '', title = "Cropland effect on beta diversity") +
  theme_few() +
  # 线型图例
  scale_linetype_manual(
    name = "Error bars",
    values = c("95% CI" = "solid", "±SE" = "solid"),
    guide = guide_legend(
      override.aes = list(
        size = c(1.5, 2.5),
        alpha = c(0.4, 0.7),
        color = "black"
      )
    )
  ) +
  scale_x_continuous(
    limits = c(-0.5, 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",  # 隐藏图例
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 18, color = "black", face = "bold", vjust = 1.9, hjust = 0.5, angle = 90),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    #panel.grid.major = element_line(color = "grey80", size = 0.1),
    axis.title.x = element_blank(),  # 隐藏x轴标题
    #panel.grid.minor = element_line(color = "grey90", size = 0.1)
  )
p




#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/effect_size_β_dispersation_cropland_effect加上vfg和mge版本.pdf",p,width=8,height=6)



#d-e DDR
library(agricolae)
library(ggplot2)
library(geosphere)
library(fdrtool)
library(vegan)
library(maps)
library(dplyr)
rm(list=ls())
arg<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-4544]

dis_mat <- as.matrix(vegdist(decostand(arg, "hellinger"), method="bray"))
dist_otu_long <- as.data.frame(as.table(dis_mat))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")

geo <- group[,c("sample.id","Longitude","Latitude","group")]  # 假设 lon, lat 列名
rownames(geo) <- geo$sample.id
geo_mat <- distm(geo[,c("Longitude","Latitude")], fun = distHaversine)/1000  # km
rownames(geo_mat) <- geo$sample.id
colnames(geo_mat) <- geo$sample.id
dist_geo_long <- as.data.frame(as.table(geo_mat))

colnames(dist_geo_long) <- c("Sample1","Sample2","dist_geo")

data <- merge(dist_otu_long, dist_geo_long, by=c("Sample1","Sample2"))

data <- merge(data, group[,c("sample.id","group")], by.x="Sample1", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group1"
data <- merge(data, group[,c("sample.id","group")], by.x="Sample2", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group2"

data_same_group <- data[data$Group1 == data$Group2, ]
data_same_group <- data_same_group[data_same_group$Sample1 != data_same_group$Sample2, ]
ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")

library(ggrepel)

library(ggplot2)
library(grid)  # 用于 unit()

env<-group[,-1]
row.names(env)<-env[,1]
env$group <- factor(env$group, labels = c("Forest", "Cropland", "Industry", "Urban"))
#计算每个组的Mantel test
mantel_results_arg <- data.frame(
  group = character(),
  R = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# 循环每个组
for (treat in levels(env$group)) {
  
  arg1 <- arg[env$group == treat, ]
  env1 <- env[env$group == treat, ]
  
  arg.dist <- vegdist(decostand(arg1, "hellinger"), "bray")
  geo.dist <- as.dist(geo_mat[rownames(env1), rownames(env1)])
  
  # Mantel test
  ma <- mantel(arg.dist, geo.dist, permutations = 999)
  
  # 保存结果
  mantel_results_arg <- rbind(
    mantel_results_arg,
    data.frame(
      group = treat,
      R = ma$statistic,
      P = ma$signif
    )
  )
}

# 查看结果
print(mantel_results_arg)
data_same_group$group3<-"arg"
data_arg<-data_same_group

#108bac
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
otu_clean <- otu[, colSums(otu) != 0]
otu<-otu_clean 


dis_mat <- as.matrix(vegdist(decostand(otu, "hellinger"), method="bray"))
dist_otu_long <- as.data.frame(as.table(dis_mat))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")

geo <- group[,c("sample.id","Longitude","Latitude","group")]  # 假设 lon, lat 列名
rownames(geo) <- geo$sample.id
geo_mat <- distm(geo[,c("Longitude","Latitude")], fun = distHaversine)/1000  # km
rownames(geo_mat) <- geo$sample.id
colnames(geo_mat) <- geo$sample.id
dist_geo_long <- as.data.frame(as.table(geo_mat))

colnames(dist_geo_long) <- c("Sample1","Sample2","dist_geo")

data <- merge(dist_otu_long, dist_geo_long, by=c("Sample1","Sample2"))

data <- merge(data, group[,c("sample.id","group")], by.x="Sample1", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group1"
data <- merge(data, group[,c("sample.id","group")], by.x="Sample2", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group2"

data_same_group <- data[data$Group1 == data$Group2, ]
data_same_group <- data_same_group[data_same_group$Sample1 != data_same_group$Sample2, ]
ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")

library(ggplot2)
library(grid)  # 用于 unit()

env<-group[,-1]
row.names(env)<-env[,1]
env$group <- factor(env$group, labels = c("Forest", "Cropland", "Industry", "Urban"))
#计算每个组的Mantel test
mantel_results_bacteria <- data.frame(
  group = character(),
  R = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# 循环每个组
for (treat in levels(env$group)) {
  
  otu1 <- otu[env$group == treat, ]
  env1 <- env[env$group == treat, ]
  
  arg.dist <- vegdist(decostand(otu1, "hellinger"), "bray")
  geo.dist <- as.dist(geo_mat[rownames(env1), rownames(env1)])
  
  # Mantel test
  ma <- mantel(arg.dist, geo.dist, permutations = 999)
  
  # 保存结果
  mantel_results_bacteria <- rbind(
    mantel_results_bacteria,
    data.frame(
      group = treat,
      R = ma$statistic,
      P = ma$signif
    )
  )
}

# 查看结果
print(mantel_results_bacteria)

data_same_group$group3<-"bacteria"
data_bacteria<-data_same_group

##high risk arg
arg<-arg <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
arg<-arg[,-1]
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-74]
arg_clean <- arg[, colSums(arg) != 0]
arg<-arg_clean 

dis_mat <- as.matrix(vegdist(decostand(arg, "hellinger"), method="bray"))
dist_otu_long <- as.data.frame(as.table(dis_mat))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")

geo <- group[,c("sample.id","Longitude","Latitude","group")]  # 假设 lon, lat 列名
rownames(geo) <- geo$sample.id
geo_mat <- distm(geo[,c("Longitude","Latitude")], fun = distHaversine)/1000  # km
rownames(geo_mat) <- geo$sample.id
colnames(geo_mat) <- geo$sample.id
dist_geo_long <- as.data.frame(as.table(geo_mat))

colnames(dist_geo_long) <- c("Sample1","Sample2","dist_geo")

data <- merge(dist_otu_long, dist_geo_long, by=c("Sample1","Sample2"))

data <- merge(data, group[,c("sample.id","group")], by.x="Sample1", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group1"
data <- merge(data, group[,c("sample.id","group")], by.x="Sample2", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group2"

data_same_group <- data[data$Group1 == data$Group2, ]
data_same_group <- data_same_group[data_same_group$Sample1 != data_same_group$Sample2, ]
ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")

library(ggplot2)
library(grid)  # 用于 unit()

env<-group[,-1]
row.names(env)<-env[,1]
env$group <- factor(env$group, labels = c("Forest", "Cropland", "Industry", "Urban"))
#计算每个组的Mantel test
mantel_results_high <- data.frame(
  group = character(),
  R = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# 循环每个组
for (treat in levels(env$group)) {
  
  arg1 <- arg[env$group == treat, ]
  env1 <- env[env$group == treat, ]
  
  arg.dist <- vegdist(decostand(arg1, "hellinger"), "bray")
  geo.dist <- as.dist(geo_mat[rownames(env1), rownames(env1)])
  
  # Mantel test
  ma <- mantel(arg.dist, geo.dist, permutations = 999)
  
  # 保存结果
  mantel_results_high <- rbind(
    mantel_results_high,
    data.frame(
      group = treat,
      R = ma$statistic,
      P = ma$signif
    )
  )
}

# 查看结果
print(mantel_results_high)
data_same_group$group3<-"high"
data_high<-data_same_group

data<-rbind(data_bacteria,data_arg,data_high)

library(ggplot2)
library(grid)  # 用于 unit()

mantel_results_high$group3<-"high"
mantel_results_bacteria$group3<-"bacteria"
mantel_results_arg$group3<-"arg"
mantel<-rbind(mantel_results_bacteria,mantel_results_arg,mantel_results_high)
mantel$combo <- paste(mantel$group3, mantel$group, sep = "_")
mantel$xpos <- max(data$dist_geo) * 0.8
mantel$ypos <- seq(max(data$dist_otu), min(data$dist_otu), length.out = nrow(mantel))

# 先生成组合变量
data$combo <- paste(data$group3, data$Group1, sep = "_")

library(ggplot2)

# 确保 combo 顺序正确
data$combo <- factor(
  data$combo,
  levels = c(
    "bacteria_Forest", "bacteria_Cropland", "bacteria_Industry", "bacteria_Urban",
    "arg_Forest", "arg_Cropland", "arg_Industry", "arg_Urban",
    "high_Forest", "high_Cropland", "high_Industry", "high_Urban"
  )
)

# 自定义颜色
combo_colors <- c(
  "#6cb509","#00FFFF","#5e9ceb","#ff7f27",  # bacteria 内的 4 个 group
  "#6cb509","#00FFFF","#5e9ceb","#ff7f27",  # ARG 内的 4 个 group
  "#6cb509","#00FFFF","#5e9ceb","#ff7f27"   # high-risk ARG 内的 4 个 group
)

# 创建图例标签，把 R 和 P 加进去
mantel <- mantel %>%
  mutate(label = paste0(combo, " (R=", round(R, 3), ", P=", P, ")"))

# 把 label 对应回 data$combo
data$combo_label <- mantel$label[match(data$combo, mantel$combo)]
# 替换文本
mantel$label <- mantel$label %>%
  gsub("^bacteria", "Bacteria", .) %>%
  gsub("^arg", "ARG", .) %>%
  gsub("^high", "High risk arg", .)

# 把 label 对应回 data$combo
data$combo_label <- mantel$label[match(data$combo, mantel$combo)]

#筛选bacteria和ARG
data1 <- data %>%
  filter(group3 %in% c("bacteria", "arg"))
p_values <- data.frame(
  combo = c("bacteria_Forest", "bacteria_Cropland", "bacteria_Industry", "bacteria_Urban","arg_Forest","arg_Cropland","arg_Industry","arg_Urban"),
  p_value = c(0.002, 0.004, 0.151, 0.041,0.002,0.012,0.352,0.460)
)

library(dplyr)
p_values <- p_values %>%
  mutate(linetype_group = ifelse(p_value >= 0.01, "dashed", "solid"))

data1 <- data1 %>%
  left_join(p_values[, c("combo", "linetype_group")], by = "combo")

data1 <- data1 %>%
  left_join(mantel[, c("combo", "R")], by = "combo")


label_positions <- data1 %>%
  group_by(combo_label) %>%
  summarize(x = mean(dist_geo, na.rm = TRUE),
            y = max(dist_otu, na.rm = TRUE) * 0.6) %>%
  left_join(mantel, by = c("combo_label" = "label")) %>%
  mutate(label_text = sprintf("R = %.2f, P = %.3f", R, P))

label_positions <- label_positions %>%
  mutate(label_text = paste0(group, ": ", label_text))

ecosystem_order <- c("bacteria_Forest", "bacteria_Cropland", "bacteria_Industry", "bacteria_Urban","arg_Forest","arg_Cropland","arg_Industry","arg_Urban")
label_positions$combo <- factor(label_positions$combo, levels = ecosystem_order)
label_positions <- label_positions[order(label_positions$combo), ]
label_positions$x <- 200
label_positions$y <- c(0.8, 0.785, 0.77, 0.755, 0.48, 0.465, 0.45, 0.435)




p <- ggplot(data1, aes(x = dist_geo, y = dist_otu, color = combo_label)) +
  # 背景色（可保留或删除）
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 0.5, ymax = Inf,
    fill = "#FA9A38", alpha = 0.1
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = 0.5,
    fill = "#6B30C5", alpha = 0.1
  ) +
  
  # 主拟合线
  geom_smooth(
    aes(linetype = linetype_group, size = abs(R)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3
  ) +
  
  # 颜色映射
  scale_color_manual(
    values = combo_colors,
    breaks = mantel$label
  ) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  #scale_size_continuous(
  #range = c(1, 4),
  #guide = "none"
  #) +
  
  # ✅ 在图上标出 R 和 P 值
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label_text, color = combo_label),
    size = 5,
    fontface = "bold",
    show.legend = FALSE,
    hjust = 0
    
  ) +
  
  labs(
    x = "Geographic distance (km)",
    y = "Compositional dissimilarity"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text = element_text(size = 18, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black", size = 0.8),
    legend.position = "none" ,  # ✅ 完全隐藏图例
    panel.grid = element_blank()#不显示网格线
  )
p


#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/total_整合距离衰减_bacteria_arg.pdf",height=12,width=6)

#FigS11
#high risk arg 的距离衰减
data2 <- data %>%
  filter(group3 %in% c("high"))
data2<-data2[,-10]
p_values <- data.frame(
  combo = c("high_Forest", "high_Cropland", "high_Industry", "high_Urban"),
  p_value = c(0.014, 0.867, 0.277, 0.638)
)

library(dplyr)
p_values <- p_values %>%
  mutate(linetype_group = ifelse(p_value >= 0.01, "dashed", "solid"))
data2 <- data2 %>%
  left_join(p_values[, c("combo", "linetype_group")], by = "combo")
data2 <- data2 %>%
  left_join(mantel[, c("combo", "R")], by = "combo")

label_positions <- data2 %>%
  group_by(combo_label) %>%
  summarize(x = mean(dist_geo, na.rm = TRUE),
            y = max(dist_otu, na.rm = TRUE) * 0.6) %>%
  left_join(mantel, by = c("combo_label" = "label")) %>%
  mutate(label_text = sprintf("R = %.2f, P = %.3f", R, P))

label_positions

label_positions <- label_positions %>%
  mutate(label_text = paste0(group, ": ", label_text))


ecosystem_order <- c("high_Forest", "high_Cropland", "high_Industry", "high_Urban")
label_positions$combo <- factor(label_positions$combo, levels = ecosystem_order)
label_positions <- label_positions[order(label_positions$combo), ]
label_positions$x <- 200
label_positions$y <- c(0.34, 0.33, 0.32, 0.31)




p <- ggplot(data2, aes(x = dist_geo, y = dist_otu, color = combo_label)) +
  geom_smooth(
    aes(linetype = linetype_group, size = abs(R)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3
  )+
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(1, 4),
    guide = "none"
  ) +
  scale_color_manual(
    values = combo_colors,
    breaks = mantel$label,
    guide = guide_legend(nrow = 2, byrow = TRUE),
    name = NULL  # 保持颜色图例无标题
  ) +
  
  # ✅ 在图上标出 R 和 P 值
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label_text, color = combo_label),
    size = 5,
    fontface = "bold",
    show.legend = FALSE,
    hjust=0
  ) +
  labs(
    x = "Geographic distance (km)",
    y = "Community dissimilarity"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text = element_text(size = 18, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black", size = 0.8),
    legend.position = "none" ,  # ✅ 完全隐藏图例
    panel.grid = element_blank()#不显示网格线
  )+
  ggtitle("High risk resistome")



p



#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/high_整合距离衰减.pdf",height=6,width=6)

#环境距离
library(agricolae)
library(ggplot2)
library(geosphere)
library(fdrtool)
library(vegan)
library(maps)
rm(list=ls())
arg<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T,row.names=1,sep=",")
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-4544]
row.names(group)<-group$sample.id

dis_mat <- as.matrix(vegdist(decostand(arg, "hellinger"), method="bray"))
dist_otu_long <- as.data.frame(as.table(dis_mat))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")

env<-group[,4:18]
row.names(env)<-group$sample.id
scale.env<-scale(env,center=TRUE,scale=TRUE)
dist.env<-as.matrix(dist(scale.env,method='euclidean'))
dist_env_long <- as.data.frame(as.table(dist.env))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")
colnames(dist_env_long) <- c("Sample1","Sample2","dist_env")


data <- merge(dist_otu_long, dist_env_long, by=c("Sample1","Sample2"))


data <- merge(data, group[,c("sample.id","group")], by.x="Sample1", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group1"
data <- merge(data, group[,c("sample.id","group")], by.x="Sample2", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group2"

data_same_group <- data[data$Group1 == data$Group2, ]
data_same_group <- data_same_group[data_same_group$Sample1 != data_same_group$Sample2, ]
ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")

library(ggrepel)

library(ggplot2)

library(grid)  # 用于 unit()

identical(rownames(group), rownames(env))
env$group<-group$group
env$group <- factor(env$group, labels = c("Forest", "Cropland", "Industry", "Urban"))
#计算每个组的Mantel test
mantel_results_arg <- data.frame(
  group = character(),
  R = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# 循环每个组
for (treat in levels(env$group)) {
  arg1 <- arg[env$group == treat, ]
  env1 <- env[env$group == treat, ]
  
  arg.dist <- vegdist(decostand(arg1, "hellinger"), "bray")
  scale.env<-scale(env1[,1:15],center=TRUE,scale=TRUE)
  dist.env<-dist(scale.env,method='euclidean')
  
  # Mantel test
  ma <- mantel(arg.dist, dist.env, permutations = 999)
  
  # 保存结果
  mantel_results_arg <- rbind(
    mantel_results_arg,
    data.frame(
      group = treat,
      R = ma$statistic,
      P = ma$signif
    )
  )
}

# 查看结果
print(mantel_results_arg)
data_same_group$group3<-"arg"
data_arg<-data_same_group

#108bac
otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T,row.names=1)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
row.names(group)<-group$sample.id

otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]
identical(rownames(otu), group$sample.id)
otu_clean <- otu[, colSums(otu) != 0]
otu<-otu_clean 


dis_mat <- as.matrix(vegdist(decostand(otu, "hellinger"), method="bray"))
dist_otu_long <- as.data.frame(as.table(dis_mat))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")

env<-group[,4:18]
row.names(env)<-group$sample.id
scale.env<-scale(env,center=TRUE,scale=TRUE)
dist.env<-as.matrix(dist(scale.env,method='euclidean'))
dist_env_long <- as.data.frame(as.table(dist.env))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")
colnames(dist_env_long) <- c("Sample1","Sample2","dist_env")


data <- merge(dist_otu_long, dist_env_long, by=c("Sample1","Sample2"))


data <- merge(data, group[,c("sample.id","group")], by.x="Sample1", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group1"
data <- merge(data, group[,c("sample.id","group")], by.x="Sample2", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group2"

data_same_group <- data[data$Group1 == data$Group2, ]
data_same_group <- data_same_group[data_same_group$Sample1 != data_same_group$Sample2, ]
ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")

library(ggrepel)

library(ggplot2)

library(grid)  # 用于 unit()

identical(rownames(group), rownames(env))
env$group<-group$group
env$group <- factor(env$group, labels = c("Forest", "Cropland", "Industry", "Urban"))
#计算每个组的Mantel test
mantel_results_bacteria <- data.frame(
  group = character(),
  R = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# 循环每个组
for (treat in levels(env$group)) {
  otu1 <- otu[env$group == treat, ]
  env1 <- env[env$group == treat, ]
  
  arg.dist <- vegdist(decostand(otu1, "hellinger"), "bray")
  scale.env<-scale(env1[,1:15],center=TRUE,scale=TRUE)
  dist.env<-dist(scale.env,method='euclidean')
  
  # Mantel test
  ma <- mantel(arg.dist, dist.env, permutations = 999)
  
  # 保存结果
  mantel_results_bacteria <- rbind(
    mantel_results_bacteria,
    data.frame(
      group = treat,
      R = ma$statistic,
      P = ma$signif
    )
  )
}

# 查看结果
print(mantel_results_bacteria)
data_same_group$group3<-"bacteria"
data_bacteria<-data_same_group

#high risk arg
arg<-arg <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
arg<-arg[,-1]
arg<-as.data.frame(t(arg))
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
group<-group[,-1]
row.names(group)<-group$sample.id
arg$sample.id <- rownames(arg)
arg <- arg[order(arg$sample.id),]
identical(rownames(arg), group$sample.id)
arg<-arg[,-74]
arg_clean <- arg[, colSums(arg) != 0]
arg<-arg_clean 
identical(rownames(arg), group$sample.id)

dis_mat <- as.matrix(vegdist(decostand(arg, "hellinger"), method="bray"))
dist_otu_long <- as.data.frame(as.table(dis_mat))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")

env<-group[,4:18]
row.names(env)<-group$sample.id
scale.env<-scale(env,center=TRUE,scale=TRUE)
dist.env<-as.matrix(dist(scale.env,method='euclidean'))
dist_env_long <- as.data.frame(as.table(dist.env))
colnames(dist_otu_long) <- c("Sample1","Sample2","dist_otu")
colnames(dist_env_long) <- c("Sample1","Sample2","dist_env")


data <- merge(dist_otu_long, dist_env_long, by=c("Sample1","Sample2"))


data <- merge(data, group[,c("sample.id","group")], by.x="Sample1", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group1"
data <- merge(data, group[,c("sample.id","group")], by.x="Sample2", by.y="sample.id")
colnames(data)[ncol(data)] <- "Group2"

data_same_group <- data[data$Group1 == data$Group2, ]
data_same_group <- data_same_group[data_same_group$Sample1 != data_same_group$Sample2, ]
ecosystem_order <- c("Forest", "Cropland", "Industry", "Urban")

library(ggrepel)

library(ggplot2)

library(grid)  # 用于 unit()

identical(rownames(group), rownames(env))
env$group<-group$group
env$group <- factor(env$group, labels = c("Forest", "Cropland", "Industry", "Urban"))
#计算每个组的Mantel test
mantel_results_high <- data.frame(
  group = character(),
  R = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# 循环每个组
for (treat in levels(env$group)) {
  arg1 <- arg[env$group == treat, ]
  env1 <- env[env$group == treat, ]
  
  arg.dist <- vegdist(decostand(arg1, "hellinger"), "bray")
  scale.env<-scale(env1[,1:15],center=TRUE,scale=TRUE)
  dist.env<-dist(scale.env,method='euclidean')
  
  # Mantel test
  ma <- mantel(arg.dist, dist.env, permutations = 999)
  
  # 保存结果
  mantel_results_high <- rbind(
    mantel_results_high,
    data.frame(
      group = treat,
      R = ma$statistic,
      P = ma$signif
    )
  )
}

# 查看结果
print(mantel_results_high)
data_same_group$group3<-"high"
data_high<-data_same_group


data<-rbind(data_bacteria,data_arg,data_high)

library(ggplot2)
library(grid)  # 用于 unit()

mantel_results_high$group3<-"high"
mantel_results_bacteria$group3<-"bacteria"
mantel_results_arg$group3<-"arg"
mantel<-rbind(mantel_results_bacteria,mantel_results_arg,mantel_results_high)
mantel$combo <- paste(mantel$group3, mantel$group, sep = "_")
mantel$xpos <- max(data$dist_env) * 0.8
mantel$ypos <- seq(max(data$dist_otu), min(data$dist_otu), length.out = nrow(mantel))

# 先生成组合变量
data$combo <- paste(data$group3, data$Group1, sep = "_")

library(ggplot2)

# 确保 combo 顺序正确
data$combo <- factor(
  data$combo,
  levels = c(
    "bacteria_Forest", "bacteria_Cropland", "bacteria_Industry", "bacteria_Urban",
    "arg_Forest", "arg_Cropland", "arg_Industry", "arg_Urban",
    "high_Forest", "high_Cropland", "high_Industry", "high_Urban"
  )
)

# 自定义颜色
combo_colors <- c(
  "#6cb509","#00FFFF","#5e9ceb","#ff7f27",  # bacteria 内的 4 个 group
  "#6cb509","#00FFFF","#5e9ceb","#ff7f27",  # ARG 内的 4 个 group
  "#6cb509","#00FFFF","#5e9ceb","#ff7f27"   # high-risk ARG 内的 4 个 group
)

library(dplyr)
# 创建图例标签，把 R 和 P 加进去
mantel <- mantel %>%
  mutate(label = paste0(combo, " (R=", round(R, 3), ", P=", P, ")"))

# 把 label 对应回 data$combo
data$combo_label <- mantel$label[match(data$combo, mantel$combo)]
# 替换文本
mantel$label <- mantel$label %>%
  gsub("^bacteria", "Bacteria", .) %>%
  gsub("^arg", "ARG", .) %>%
  gsub("^high", "High risk arg", .)

# 把 label 对应回 data$combo
data$combo_label <- mantel$label[match(data$combo, mantel$combo)]
#筛选bacteria和arg的的数据
data1 <- data %>%
  filter(group3 %in% c("bacteria", "arg"))
p_values <- data.frame(
  combo = c("bacteria_Forest", "bacteria_Cropland", "bacteria_Industry", "bacteria_Urban","arg_Forest","arg_Cropland","arg_Industry","arg_Urban"),
  p_value = c(0.001, 0.001, 0.004, 0.021,0.001,0.001,0.206,0.625)
)

library(dplyr)
p_values <- p_values %>%
  mutate(linetype_group = ifelse(p_value > 0.01, "dashed", "solid"))
data1 <- data1 %>%
  left_join(p_values[, c("combo", "linetype_group")], by = "combo")

data1 <- data1 %>%
  left_join(mantel[, c("combo", "R")], by = "combo")

label_positions <- data1 %>%
  group_by(combo_label) %>%
  summarize(x = mean(dist_env, na.rm = TRUE),
            y = max(dist_otu, na.rm = TRUE) * 0.6) %>%
  left_join(mantel, by = c("combo_label" = "label")) %>%
  mutate(label_text = sprintf("R = %.2f, P = %.3f", R, P))
label_positions <- label_positions %>%
  mutate(label_text = paste0(group, ": ", label_text))

label_positions
ecosystem_order <- c("bacteria_Forest", "bacteria_Cropland", "bacteria_Industry", "bacteria_Urban","arg_Forest","arg_Cropland","arg_Industry","arg_Urban")
label_positions$combo <- factor(label_positions$combo, levels = ecosystem_order)
label_positions <- label_positions[order(label_positions$combo), ]
label_positions$x <- 1
label_positions$y <- c(0.86, 0.84, 0.82, 0.8, 0.48, 0.46, 0.44, 0.42)


# 绘图
p <- ggplot(data1, aes(x = dist_env, y = dist_otu, color = combo_label)) +
  # 主拟合线
  geom_smooth(
    aes(linetype = linetype_group, size = abs(R)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 0.5, ymax = Inf,
    fill = "#FA9A38", alpha = 0.1
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = 0.5,
    fill = "#6B30C5", alpha = 0.1
  ) +
  #scale_size_continuous(
  #range = c(1, 4),
  #guide = "none"
  #) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  # ✅ 在图上标出 R 和 P 值
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label_text, color = combo_label),
    size = 5,
    fontface = "bold",
    show.legend = FALSE,
    hjust=0
  ) +
  scale_color_manual(
    values = combo_colors,
    breaks = mantel$label,
    guide = guide_legend(nrow = 4, byrow = TRUE),
    name = NULL  # 保持颜色图例无标题
  ) +
  labs(x = "Environmental difference", y = "Compositional dissimilarity", color = NULL) +  # 不显示图例标题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text = element_text(size = 18, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black", size = 0.8),
    legend.position = "none" ,  # ✅ 完全隐藏图例
    panel.grid = element_blank()#不显示网格线
  )

p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/bacteria_arg_整合环境距离.pdf",height=12,width=6)

#Fig.S11
##high risk arg 的环境衰减

data2 <- data %>%
  filter(group3 %in% c("high"))

p_values <- data.frame(
  combo = c("high_Forest", "high_Cropland", "high_Industry", "high_Urban"),
  p_value = c(0.010, 0.315,0.024,  0.795)
)

library(dplyr)
p_values <- p_values %>%
  mutate(linetype_group = ifelse(p_value >= 0.01, "dashed", "solid"))
data2 <- data2 %>%
  left_join(p_values[, c("combo", "linetype_group")], by = "combo")

data2 <- data2 %>%
  left_join(mantel[, c("combo", "R")], by = "combo")
label_positions <- data2 %>%
  group_by(combo_label) %>%
  summarize(x = mean(dist_env, na.rm = TRUE),
            y = max(dist_otu, na.rm = TRUE) * 0.6) %>%
  left_join(mantel, by = c("combo_label" = "label")) %>%
  mutate(label_text = sprintf("R = %.2f, P = %.3f", R, P))

label_positions
label_positions <- label_positions %>%
  mutate(label_text = paste0(group, ": ", label_text))
ecosystem_order <- c("high_Forest", "high_Cropland", "high_Industry", "high_Urban")
label_positions$combo <- factor(label_positions$combo, levels = ecosystem_order)
label_positions <- label_positions[order(label_positions$combo), ]
label_positions$x <- 1
label_positions$y <- c(0.32, 0.312, 0.304, 0.296)

p <- ggplot(data2, aes(x = dist_env, y = dist_otu, color = combo_label)) +
  # 主拟合线
  geom_smooth(
    aes(linetype = linetype_group, size = abs(R)),
    method = "lm",
    se = TRUE,
    fill = "grey90",
    alpha = 0.3
  ) +
  scale_size_continuous(
    range = c(1, 4),
    guide = "none"
  ) +
  # ✅ 在图上标出 R 和 P 值
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label_text, color = combo_label),
    size = 5,
    fontface = "bold",
    show.legend = FALSE,
    hjust=0
  ) +
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "42"),
    guide = "none"
  ) +
  scale_color_manual(
    values = combo_colors,
    breaks = mantel$label,
    guide = guide_legend(nrow = 2, byrow = TRUE),
    name = NULL  # 保持颜色图例无标题
  ) +
  labs(x = "Environmental difference", y = "Compositional dissimilarity", color = NULL) +  # 不显示图例标题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text = element_text(size = 18, color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black", size = 0.8),
    legend.position = "none" ,  # ✅ 完全隐藏图例
    panel.grid = element_blank()#不显示网格线
  )+
  ggtitle("High risk resistome")

p

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/high_整合环境距离.pdf",p,height=6,width=6)


#Fig.4f
##HGT 事件与β多样性
library(betapart)
library(agricolae)
library(ggplot2)
library(patchwork)
library(vegan)
library(reshape2)
library(Cairo)
rm(list=ls())
arg_all0<- read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",sep=",",header=T,row.names=1)

arg_all0 <- as.data.frame(t(arg_all0))

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg_all0$sample.id <- rownames(arg_all0)
arg_all0 <- arg_all0[order(arg_all0$sample.id),]

colnames(arg_all0)[(ncol(arg_all0)-9):ncol(arg_all0)]
dis <- vegdist(decostand(arg_all0[,1:4543],"hellinger"), method = 'bray')


Location <- group$Location
Group <- group$group
mod <- betadisper(dis,Group, type = 'centroid')

distances_to_centroid <- mod$distances
anova_data <- data.frame(
  Distance = distances_to_centroid,
  Location = Location,
  Group = Group
)
anova_data$sample.id<-row.names(anova_data)


#hgt
hgt<-read.csv("/Users/yangy/Documents/E/city/metagenome/HGT/results/hgt_arg_gene_count.csv",sep=",",header=T,row.names=1)
merge<-merge(anova_data,hgt,by="sample.id")
colnames(merge)[4]<-"group"

hgt1<-read.csv("/Users/yangy/Documents/E/city/metagenome/HGT/results/hgt_arg_gene_richness.csv",sep=",",header=T,row.names=1)
merge<-merge(merge,hgt1,by="sample.id")


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
    cor.test(group_data$Distance, group_data$count, method = "spearman")
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
label_Cropland<-paste0("Cropland: ρ = ",cropland_r_value, ", P= ", cropland_P_value)

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
                        " P = ", sapply(P_value, format_p))
  )

ecosystem_order<-c("Forest","Cropland","Industry","Urban")




p_values <- data.frame(
  group = c("Cropland", "Industry", "Forest", "Urban"),
  p_value = c(0.8341, 0.1061, 0.9343, 0.0269)
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
results$x <- c(8, 8, 8, 8)
results$y<-c(0.42, 0.405,0.39, 0.375)
results <- results %>%
  mutate(label_text = paste0(group, ": ", label_text))

p <- ggplot(merge, aes(x = count, y = Distance, color = group)) +
  geom_jitter(size = 3, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = " ",
    x = "HGT events-ARG", y = "Beta diversity of resistome"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # 主网格线
    panel.grid.minor = element_line(color = "grey90", size = 0.3)   # 次网格线
  ) +
  # 对每个 group 绘制趋势线
  geom_smooth(
    aes(x = count, y = Distance, color = group, linetype = linetype_group, size = abs(R_value)),
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


#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/arg_hgt_count_beta.pdf",p,device = cairo_pdf,width=6,height=6)


#Fig.4g
library(betapart)
library(agricolae)
library(ggplot2)
library(patchwork)
library(vegan)
library(reshape2)
library(Cairo)
rm(list=ls())
arg_all0<- read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",sep=",",header=T,row.names=1)

arg_all0 <- as.data.frame(t(arg_all0))

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
arg_all0$sample.id <- rownames(arg_all0)
arg_all0 <- arg_all0[order(arg_all0$sample.id),]

colnames(arg_all0)[(ncol(arg_all0)-9):ncol(arg_all0)]
dis <- vegdist(decostand(arg_all0[,1:4543],"hellinger"), method = 'bray')


Location <- group$Location
Group <- group$group
mod <- betadisper(dis,Group, type = 'centroid')

distances_to_centroid <- mod$distances
anova_data <- data.frame(
  Distance = distances_to_centroid,
  Location = Location,
  Group = Group
)
anova_data$sample.id<-row.names(anova_data)

group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
mge<-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
mge<-t(mge)
mge <- mge[order(rownames(mge)), ]
identical(rownames(mge), group$sample.id)
group$mge_abundance<-rowSums(mge)


merge<-merge(anova_data,group,by="sample.id")



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
    cor.test(group_data$Distance, group_data$mge_abundance, method = "spearman")
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
library(dplyr)
results <- results %>%
  mutate(
    label_text = paste0("ρ = ", round(Spearman_Correlation, 2),
                        " P = ", sapply(P_value, format_p))
  )

ecosystem_order<-c("Forest","Cropland","Industry","Urban")


p_values <- data.frame(
  group = c("Cropland", "Industry", "Forest", "Urban"),
  p_value = c(0.22466415, 0.23318172, 0.83475111, 0.01773069)
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

results$x <- c(0.3, 0.3, 0.3, 0.3)
results$y<-c(0.44, 0.425,0.41, 0.395)
results <- results %>%
  mutate(label_text = paste0(group, ": ", label_text))




p <- ggplot(merge, aes(x = mge_abundance, y = Distance, color = group)) +
  geom_jitter(size = 3, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = " ",
    x = "Copies of MGEs per cell", y = "Beta diversity of resistome"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # 主网格线
    panel.grid.minor = element_line(color = "grey90", size = 0.3)   # 次网格线
  ) +
  # 对每个 group 绘制趋势线
  geom_smooth(
    aes(x = mge_abundance, y = Distance, color = group, linetype = linetype_group, size = abs(R_value)),
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
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/beta_diversity/arg_mge_copy_beta_diversity.pdf",p,width=6,height=6,device = cairo_pdf)
