#Fig.5
#vfg network
library(GO.db)
library(vegan)
library(permute)        
library(lattice)
library(WGCNA)
library(igraph)
library(WGCNA)
library(psych)
library(reshape2)
library(igraph)
library(Hmisc)
library(GUniFrac)
library(sciplot)
library(ggpmisc)
library(indicspecies)
library(BiocManager)
library(patchwork)
library(agricolae)
library(ggplot2)
library(ggraph)
library(colorRamps)
rm(list=ls())

vfg<-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
vfg<-t(as.data.frame(vfg))
rs.d<-vfg[,colSums(vfg)>0]
group<-read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv", header=TRUE)
vfgid<-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE)


rs_forest<-rs.d[group$group%in%c("Forest"),]
rs_forest1<-rs_forest[,colSums(rs_forest)>0]
rs_forest1<-as.data.frame(t(rs_forest1))
rs_forest1 <- t(rs_forest1)[, (colSums(t(rs_forest1)) / sum(colSums(t(rs_forest1)))) >= 0.0001]
trs_forest1<-t(rs_forest1)
trs_forest1<-data.frame(row.names(trs_forest1),trs_forest1)
names(trs_forest1)[1]<-c("gene")
forest<-merge(trs_forest1,vfgid,by="gene")
forest_order<-forest[order(forest$Type),]
forest_order1<-forest_order[-29]
forest_order2<-data.frame(forest_order1)
#write.csv(forest_order2,file = "forest_order3.csv")
rownames(forest_order2)<-forest_order2[,1]
forest_order3<-forest_order2[,-1]



rs_cropland<-rs.d[group$group%in%c("Cropland"),]
rs_cropland1<-rs_cropland[,colSums(rs_cropland)>0]
rs_cropland1<-as.data.frame(t(rs_cropland1))
rs_cropland1 <- t(rs_cropland1)[, (colSums(t(rs_cropland1)) / sum(colSums(t(rs_cropland1)))) >= 0.0001]
trs_cropland1<-t(rs_cropland1)
trs_cropland1<-data.frame(row.names(trs_cropland1),trs_cropland1)
names(trs_cropland1)[1]<-c("gene")
cropland<-merge(trs_cropland1,vfgid,by="gene")
cropland_order<-cropland[order(cropland$Type),]
cropland_order1<-cropland_order[-29]
cropland_order2<-data.frame(cropland_order1)
#write.csv(cropland_order2,file = "cropland_order3.csv")
rownames(cropland_order2)<-cropland_order2[,1]
cropland_order3<-cropland_order2[,-1]


rs_industry<-rs.d[group$group%in%c("Industry"),]
rs_industry1<-rs_industry[,colSums(rs_industry)>0]
rs_industry1<-as.data.frame(t(rs_industry1))
rs_industry1 <- t(rs_industry1)[, (colSums(t(rs_industry1)) / sum(colSums(t(rs_industry1)))) >= 0.0001]
trs_industry1<-t(rs_industry1)
trs_industry1<-data.frame(row.names(trs_industry1),trs_industry1)
names(trs_industry1)[1]<-c("gene")
industry<-merge(trs_industry1,vfgid,by="gene")
industry_order<-industry[order(industry$Type),]
industry_order1<-industry_order[-29]
industry_order2<-data.frame(industry_order1)
#write.csv(industry_order2,file = "industry_order3.csv")
rownames(industry_order2)<-industry_order2[,1]
industry_order3<-industry_order2[,-1]

rs_urban<-rs.d[group$group%in%c("Urban"),]
rs_urban1<-rs_urban[,colSums(rs_urban)>0]
rs_urban1<-as.data.frame(t(rs_urban1))
rs_urban1 <- t(rs_urban1)[, (colSums(t(rs_urban1)) / sum(colSums(t(rs_urban1)))) >= 0.0001]
trs_urban1<-t(rs_urban1)
trs_urban1<-data.frame(row.names(trs_urban1),trs_urban1)
names(trs_urban1)[1]<-c("gene")
urban<-merge(trs_urban1,vfgid,by="gene")
urban_order<-urban[order(urban$Type),]
urban_order1<-urban_order[-29]
urban_order2<-data.frame(urban_order1)
#write.csv(urban_order2,file = "urban_order3.csv")
rownames(urban_order2)<-urban_order2[,1]
urban_order3<-urban_order2[,-1]



############### forest
df1<-forest_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_forest<-igraph
#igraph<-igraph_forest
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))


E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####
par(mfrow=c(1,4))

data2=read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE,row.names=1)
igraph_col = data2[V(igraph)$name,]

igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","darkgreen","#FFFF05", "grey","pink",   "black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_forest<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_forest<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                    diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)
set.seed(123)

forest <- plot(
  igraph,
  main = "Forest\n Nodes: 879\n Edges: 18687\n Average degree: 42.5",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)





##cropland
df1<-cropland_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_cropland<-igraph
#igraph<-igraph_cropland
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight


#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0,"grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE,row.names=1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","darkgreen","#FFFF05", "grey","pink",   "black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_cropland<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_cropland<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)

cropland <- plot(
  igraph,
  main = "Cropland\n Nodes: 850\n Edges: 11832\n Average degree: 27.8",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)




##industry

df1<-industry_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_industry<-igraph 

tax = read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))
#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE,row.names=1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

color = c("red","orange" ,"green","cyan", "#0000FF","#7506F7","darkgreen","#FFFF05", "grey","pink",   "black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)


num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_industry<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)

prold_industry<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)
set.seed(123)

industry<- plot(
  igraph,
  main = "Industry\n Nodes: 731\n Edges: 12153\n Average degree: 33.3",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)


##Urban
df1<-urban_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_urban<-igraph
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0,"grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#7506F7","darkgreen","#FFFF05", "grey","pink",   "black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)



num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_urban <- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_urban<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                   diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)


urban<- plot(
  igraph,
  main = "Urban\n Nodes: 877\n Edges: 38759\n Average degree: 88.4",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)
###########################################################################property_network

proall<-rbind(prold_forest,prold_cropland,prold_industry,prold_urban)
proall_network<-data.frame(proall)
proall_network$group<-c('Forest','Cropland','Industry','Urban')

##################################################################################Contribution to Degree

node_degree_df_forest$group<-rep("Forest",nrow(node_degree_df_forest))
node_degree_df_cropland$group<-rep("Cropland",nrow(node_degree_df_cropland))
node_degree_df_industry$group<-rep("Industry",nrow(node_degree_df_industry))
node_degree_df_urban$group<-rep("Urban",nrow(node_degree_df_urban))

degree_all<-rbind(node_degree_df_forest,node_degree_df_cropland,node_degree_df_industry,node_degree_df_urban)
degree_all$group<-factor(degree_all$group,levels = c('Forest','Cropland','Industry','Urban'))
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE,row.names = 1)
degree_all$Type <- data2[degree_all$node, "Type"]

col_map <- c(
  "Adherence" = "red",
  "Antimicrobial activity/Competitive advantage"              = "orange",
  "Biofilm"                  = "green",
  "Effector delivery system"                   = "cyan",
  "Exoenzyme"                  = "#0000FF",
  "Exotoxin"                   = "#F706DB",
  "Immune modulation"="#7506F7",
  "Motility"                  = "darkgreen",
  "Nutritional/Metabolic factor"                 = "#FFFF05",      # Others 固定为灰色
  "Others"                   = "grey",
  "Regulation"                   = "pink",
  "Stress survival"            = "black"
)
degree_all$Type <- factor(
  degree_all$Type,
  levels = c("Others", setdiff(unique(degree_all$Type), "Others"))
)


f_degree_contri <- ggplot(degree_all, aes(x = group, y = degree, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = " ", y = "Contribution to Degree") +
  labs(fill = "  ") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    panel.spacing = unit(0, "lines"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text  = element_text(colour = "black", size = 15, face = "bold"),
    axis.text = element_text(size = 15, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 15, face = "bold")
  ) +
  scale_fill_manual(values = col_map)

# 输出图
f_degree_contri

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/contribution_to_degree_vfg_all.pdf",f_degree_contri,heigh=4,width=8)

##正负连接图
#正负连接饼图
#forest
library(ggplot2)
library(ggsci)

# 示例数据
pie_forest <- data.frame(
  links = c("positive", "negative"),
  value = c(17102, 1585)
)
library(dplyr)

pie_forest <- pie_forest %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_forest <- ggplot(pie_forest, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_forest

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_forest_vfg.pdf",p1_forest,heigh=5,width=6)

library(ggplot2)
library(ggsci)

# 示例数据
pie_cropland <- data.frame(
  links = c("positive", "negative"),
  value = c(11739, 93)
)
library(dplyr)

pie_cropland <- pie_cropland %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_cropland <- ggplot(pie_cropland, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_cropland
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_cropland_vfg.pdf",p1_cropland,heigh=5,width=6)

#industry
library(ggplot2)
library(ggsci)

# 示例数据
pie_industry <- data.frame(
  links = c("positive", "negative"),
  value = c(11707, 446)
)
library(dplyr)

pie_industry <- pie_industry %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_industry <- ggplot(pie_industry, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_industry
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_industry_vfg.pdf",p1_industry,heigh=5,width=6)


#urban
library(ggplot2)
library(ggsci)

# 示例数据
pie_urban <- data.frame(
  links = c("positive", "negative"),
  value = c(38012, 747)
)
library(dplyr)

pie_urban <- pie_urban %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_urban <- ggplot(pie_urban, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_urban
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_urban_vfg.pdf",p1_urban,heigh=5,width=6)
#VFG鲁棒性
library(vegan)
library(psych)
library(igraph)
library(ggplot2)
library(colorRamps)
library(ape)
library(stringr)
library(reshape2)
library(splitstackshape)
library(emmeans)
library(graphics)
library(phyloseq)
library(dplyr)
library(ggraph)
library(pals)
library(RColorBrewer)
library(patchwork)
library(scales) 
rm(list=ls())

bacteria.network<-function(k){
  otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
  otu_all <- otu_all[, order(colnames(otu_all))]
  group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
  env0 <- group2[order(group2$sample.id), ]
  identical(group2$sample.id, colnames(otu_all))
  
  ID0 <-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/all_vfg_ID_network.csv", header=TRUE)
  # Set node color
  test <- as.data.frame(rowSums(otu_all))
  test$id <- rownames(test)
  test2 <- test[order(test$`rowSums(otu_all)`,decreasing = T),]
  otu_all <- otu_all[test2[test2$`rowSums(otu_all)` >= sum(test2$`rowSums(otu_all)`)*0.000,]$id,]
  # Subset data
  env.tmp <- env0[env0$group == k, ]
  otu_all2 <- otu_all[,colnames(otu_all) %in% env.tmp$sample.id]
  
  ID.tmp0<-ID0[ID0$subtype %in% rownames(otu_all),]
  
  otu <- t(otu_all2)[, (colSums(t(otu_all2)) / sum(colSums(t(otu_all2)))) >= 0.0001]
  otu2 <- otu_all2[colnames(otu),]
  otu2$otusum <- rowSums(otu2)  
  otu2$abu <- otu2$otusum/sum(otu2$otusum)
  spman.r0 = corr.test(otu, use="pairwise",method="spearman",adjust="fdr", alpha = .05, ci=FALSE)
  
  r_value<-spman.r0$r #提取相关系数
  p_value<-spman.r0$p
  r_value[p_value>0.01|abs(r_value)<0.6] = 0
  r_value1<<-r_value
  Cor<-as.matrix(spman.r0$r)  
  Cor[lower.tri(Cor)] <- NA
  diag(Cor) <- NA
  Cor.df <- reshape2::melt(Cor, na.rm = TRUE) 
  
  P0<-as.matrix(spman.r0$p)
  P0[lower.tri(P0)] <- NA
  diag(P0) <- NA
  P.df <- reshape2::melt(P0, na.rm = TRUE)
  
  df <- merge(Cor.df, P.df, by=c("Var1","Var2"))
  names(df) <- c("row","col","Cor","p")
  
  da.tmp <- df.sig <- filter(df, abs(Cor) >= 0.6, p < 0.01)
  
  da.g1 <- data.frame(da.tmp)
  #da.g1 <- na.omit(da.g1)
  da.g2 <<- da.g1
  g <- graph_from_data_frame(da.g1, 
                             directed=FALSE,
                             vertices = NULL)
  g1 <<- g
}


bacteria.network("Forest")
da.g1_forest<-da.g2
r_forest<-r_value1
r_upper <- r_forest
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_forest<-g1

bacteria.network("Cropland")
da.g1_cropland<-da.g2
r_cropland<-r_value1
r_upper <- r_cropland
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_cropland<-g1

bacteria.network("Industry")
da.g1_industry<-da.g2
r_industry<-r_value1
r_upper <- r_industry
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_industry<-g1

bacteria.network("Urban")
da.g1_urban<-da.g2
r_urban<-r_value1
r_upper <- r_urban
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_urban<-g1


source("/Users/yangy/Documents/E/city/bacteria/bacteria1/code/info.centrality.R")#自定义函数
cormatrix2<- r_forest
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_forest <- env0  %>% filter(group == "Forest")

otu1_forest <- otu1[group2_forest$sample.id,]
otu.rel <- sweep(otu1_forest, 1, rowSums(otu1_forest), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_forest) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_forest<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                          rbind(Weighted.simu, Unweighted.simu),
                          weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("Forest", 40)#根据自己的处理修改treat名称
)



###计算Cropland
cormatrix2<- r_cropland
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_cropland <- env0  %>% filter(group == "Cropland")

otu1_cropland <- otu1[group2_cropland$sample.id,]
otu.rel <- sweep(otu1_cropland, 1, rowSums(otu1_cropland), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_cropland) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_cropland<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("cropland", 40)#根据自己的处理修改treat名称
)


cormatrix2<- r_industry
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_industry <- env0  %>% filter(group == "Industry")

otu1_industry <- otu1[group2_industry$sample.id,]
otu.rel <- sweep(otu1_industry, 1, rowSums(otu1_industry), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_industry) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_industry<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("industry", 40)#根据自己的处理修改treat名称
)



cormatrix2<- r_urban
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/vfdb/results/normalized_level1.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_urban <- env0  %>% filter(group == "Urban")

otu1_urban <- otu1[group2_urban$sample.id,]
otu.rel <- sweep(otu1_urban, 1, rowSums(otu1_urban), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_urban) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_urban<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                         rbind(Weighted.simu, Unweighted.simu),
                         weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("urban", 40)#根据自己的处理修改treat名称
)
df_all <- rbind(dat1_forest, dat1_cropland, dat1_industry, dat1_urban)

df_all$treat <- dplyr::recode(df_all$treat,
                              cropland = "Cropland",
                              industry = "Industry",
                              urban    = "Urban")


#生成非加权网络的结果图
p_Robustness_unweighted <- ggplot(df_all[df_all$weighted == "unweighted",], 
                                  aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +  # 注意 theme_light 要放前面，避免覆盖自定义 theme
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),   # 加粗 x 轴刻度
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),  # 加粗 x 轴标题
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "Virulome-Robustness-Robustness(Unweighted)")
p_Robustness_unweighted 

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/VFG_network_鲁棒性unweighted.pdf",p_Robustness_unweighted,width = 7, height = 6)


##生成加权网络的结果图
p_Robustness_weighted <- ggplot(df_all[df_all$weighted == "weighted",], 
                                aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "Virulome-Robustness(Weighted)")


p_Robustness_weighted
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/VFGs_network_鲁棒性weighted.pdf",p_Robustness_unweighted,width = 7, height = 6)






#MGE
##MGE network
rm(list=ls())

mge<-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
mge<-t(as.data.frame(mge))
rs.d<-mge[,colSums(mge)>0]
group<-read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv", header=TRUE)
mgeid<-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE)


rs_forest<-rs.d[group$group%in%c("Forest"),]
rs_forest1<-rs_forest[,colSums(rs_forest)>0]
rs_forest1<-as.data.frame(t(rs_forest1))
rs_forest1 <- t(rs_forest1)[, (colSums(t(rs_forest1)) / sum(colSums(t(rs_forest1)))) >= 0.0001]
trs_forest1<-t(rs_forest1)
trs_forest1<-data.frame(row.names(trs_forest1),trs_forest1)
names(trs_forest1)[1]<-c("gene")
forest<-merge(trs_forest1,mgeid,by="gene")
forest_order<-forest[order(forest$Type),]
forest_order1<-forest_order[-29]
forest_order2<-data.frame(forest_order1)
#write.csv(forest_order2,file = "forest_order3.csv")
rownames(forest_order2)<-forest_order2[,1]
forest_order3<-forest_order2[,-1]



rs_cropland<-rs.d[group$group%in%c("Cropland"),]
rs_cropland1<-rs_cropland[,colSums(rs_cropland)>0]
rs_cropland1<-as.data.frame(t(rs_cropland1))
rs_cropland1 <- t(rs_cropland1)[, (colSums(t(rs_cropland1)) / sum(colSums(t(rs_cropland1)))) >= 0.0001]
trs_cropland1<-t(rs_cropland1)
trs_cropland1<-data.frame(row.names(trs_cropland1),trs_cropland1)
names(trs_cropland1)[1]<-c("gene")
cropland<-merge(trs_cropland1,mgeid,by="gene")
cropland_order<-cropland[order(cropland$Type),]
cropland_order1<-cropland_order[-29]
cropland_order2<-data.frame(cropland_order1)
#write.csv(cropland_order2,file = "cropland_order3.csv")
rownames(cropland_order2)<-cropland_order2[,1]
cropland_order3<-cropland_order2[,-1]


rs_industry<-rs.d[group$group%in%c("Industry"),]
rs_industry1<-rs_industry[,colSums(rs_industry)>0]
rs_industry1<-as.data.frame(t(rs_industry1))
rs_industry1 <- t(rs_industry1)[, (colSums(t(rs_industry1)) / sum(colSums(t(rs_industry1)))) >= 0.0001]
trs_industry1<-t(rs_industry1)
trs_industry1<-data.frame(row.names(trs_industry1),trs_industry1)
names(trs_industry1)[1]<-c("gene")
industry<-merge(trs_industry1,mgeid,by="gene")
industry_order<-industry[order(industry$Type),]
industry_order1<-industry_order[-29]
industry_order2<-data.frame(industry_order1)
#write.csv(industry_order2,file = "industry_order3.csv")
rownames(industry_order2)<-industry_order2[,1]
industry_order3<-industry_order2[,-1]

rs_urban<-rs.d[group$group%in%c("Urban"),]
rs_urban1<-rs_urban[,colSums(rs_urban)>0]
rs_urban1<-as.data.frame(t(rs_urban1))
rs_urban1 <- t(rs_urban1)[, (colSums(t(rs_urban1)) / sum(colSums(t(rs_urban1)))) >= 0.0001]
trs_urban1<-t(rs_urban1)
trs_urban1<-data.frame(row.names(trs_urban1),trs_urban1)
names(trs_urban1)[1]<-c("gene")
urban<-merge(trs_urban1,mgeid,by="gene")
urban_order<-urban[order(urban$Type),]
urban_order1<-urban_order[-29]
urban_order2<-data.frame(urban_order1)
#write.csv(urban_order2,file = "urban_order3.csv")
rownames(urban_order2)<-urban_order2[,1]
urban_order3<-urban_order2[,-1]



############### forest
df1<-forest_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_forest<-igraph
#igraph<-igraph_forest
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))


E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####
par(mfrow=c(1,4))

data2=read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE,row.names=1)
igraph_col = data2[V(igraph)$name,]

igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color <- c(
  "red",
  "orange",
  "green",
  "cyan",
  "#0000FF",
  "#F706DB",
  "#7506F7",
  "grey",
  "darkgreen",
  "#FFFF05",
  "pink"
)

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_forest<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_forest<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                    diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)
set.seed(123)

forest <- plot(
  igraph,
  main = "Forest\n Nodes: 280\n Edges: 823\n Average degree: 5.9",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)





##cropland
df1<-cropland_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_cropland<-igraph
#igraph<-igraph_cropland
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight


#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0,"grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color <- c(
  "red",
  "orange",
  "green",
  "cyan",
  "#0000FF",
  "#F706DB",
  "#7506F7",
  "grey",
  "darkgreen",
  "#FFFF05",
  "pink"
)
levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_cropland<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_cropland<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)

cropland <- plot(
  igraph,
  main = "Cropland\n Nodes: 460\n Edges: 1489\n Average degree: 6.5",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)




##industry

df1<-industry_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_industry<-igraph 

tax = read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))
#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE,row.names=1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

color <- c(
  "red",
  "orange",
  "green",
  "cyan",
  "#0000FF",
  "#F706DB",
  "#7506F7",
  "grey",
  "darkgreen",
  "#FFFF05",
  "pink"
)
levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)


num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_industry<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)

prold_industry<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)
set.seed(123)

industry<- plot(
  igraph,
  main = "Industry\n Nodes: 380\n Edges: 1173\n Average degree: 6.2",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)


##Urban
df1<-urban_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_urban<-igraph
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0,"grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE,row.names=1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color <- c(
  "red",
  "orange",
  "cyan",
  "#0000FF",
  "#F706DB",
  "#7506F7",
  "grey",
  "darkgreen",
  "#FFFF05",
  "pink"
)

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)



num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_urban <- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_urban<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                   diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)


urban<- plot(
  igraph,
  main = "Urban\n Nodes: 565\n Edges: 22424\n Average degree: 79.4",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.1,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)
)





###########################################################################property_network

proall<-rbind(prold_forest,prold_cropland,prold_industry,prold_urban)
proall_network<-data.frame(proall)
proall_network$group<-c('Forest','Cropland','Industry','Urban')

##################################################################################Contribution to Degree

node_degree_df_forest$group<-rep("Forest",nrow(node_degree_df_forest))
node_degree_df_cropland$group<-rep("Cropland",nrow(node_degree_df_cropland))
node_degree_df_industry$group<-rep("Industry",nrow(node_degree_df_industry))
node_degree_df_urban$group<-rep("Urban",nrow(node_degree_df_urban))

degree_all<-rbind(node_degree_df_forest,node_degree_df_cropland,node_degree_df_industry,node_degree_df_urban)
degree_all$group<-factor(degree_all$group,levels = c('Forest','Cropland','Industry','Urban'))
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE,row.names=1)
degree_all$Type <- data2[degree_all$node, "Type"]


# 固定颜色映射
col_map <- c(
  "insertion_element_IS91" = "red",
  "integrase"              = "orange",
  "ISRj1"                  = "green",
  "istA"                   = "cyan",
  "istA3"                  = "#0000FF",
  "istB"                   = "#F706DB",
  "istB1"                  = "#7506F7",
  "Others"                 = "grey",      # Others 固定为灰色
  "tniA"                   = "darkgreen",
  "tniB"                   = "#FFFF05",
  "transposase"            = "pink"
)

# ------- 图例排序：将 Others 放在第一个 ----------
degree_all$Type <- factor(
  degree_all$Type,
  levels = c("Others", setdiff(unique(degree_all$Type), "Others"))
)

# -------------------- 绘图 ---------------------
f_degree_contri <- ggplot(degree_all, aes(x = group, y = degree, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = " ", y = "Contribution to Degree") +
  labs(fill = "  ") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    panel.spacing = unit(0, "lines"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text  = element_text(colour = "black", size = 15, face = "bold"),
    axis.text = element_text(size = 15, face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 15, face = "bold")
  ) +
  scale_fill_manual(values = col_map)

# 输出图
f_degree_contri



#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/contribution_to_degree_mge_all.pdf",f_degree_contri,heigh=4,width=8)




##正负连接图
#正负连接饼图
#forest
library(ggplot2)
library(ggsci)

# 示例数据
pie_forest <- data.frame(
  links = c("positive", "negative"),
  value = c(745, 78)
)
library(dplyr)

pie_forest <- pie_forest %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_forest <- ggplot(pie_forest, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_forest

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_forest_mge.pdf",p1_forest,heigh=5,width=6)

library(ggplot2)
library(ggsci)

# 示例数据
pie_cropland <- data.frame(
  links = c("positive", "negative"),
  value = c(1482, 7)
)
library(dplyr)

pie_cropland <- pie_cropland %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_cropland <- ggplot(pie_cropland, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_cropland
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_cropland_mge.pdf",p1_cropland,heigh=5,width=6)

#industry
library(ggplot2)
library(ggsci)

# 示例数据
pie_industry <- data.frame(
  links = c("positive", "negative"),
  value = c(1172, 1)
)
library(dplyr)

pie_industry <- pie_industry %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_industry <- ggplot(pie_industry, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_industry
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_industry_mge.pdf",p1_industry,heigh=5,width=6)


library(ggplot2)
library(ggsci)

# 示例数据
pie_urban <- data.frame(
  links = c("positive", "negative"),
  value = c(22409, 15)
)
library(dplyr)

pie_urban <- pie_urban %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_urban <- ggplot(pie_urban, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_urban
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/positive_negative_link_urban_mge.pdf",p1_urban,heigh=5,width=6)

#MGE鲁棒性


library(vegan)
library(psych)
library(igraph)
library(ggplot2)
library(colorRamps)
library(ape)
library(stringr)
library(reshape2)
library(splitstackshape)
library(emmeans)
library(graphics)
library(phyloseq)
library(dplyr)
library(ggraph)
library(pals)
library(RColorBrewer)
library(patchwork)
library(scales) 
rm(list=ls())

bacteria.network<-function(k){
  otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
  otu_all <- otu_all[, order(colnames(otu_all))]
  group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
  env0 <- group2[order(group2$sample.id), ]
  identical(group2$sample.id, colnames(otu_all))
  
  ID0 <-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/all_mge_ID_network.csv", header=TRUE)
  # Set node color
  test <- as.data.frame(rowSums(otu_all))
  test$id <- rownames(test)
  test2 <- test[order(test$`rowSums(otu_all)`,decreasing = T),]
  otu_all <- otu_all[test2[test2$`rowSums(otu_all)` >= sum(test2$`rowSums(otu_all)`)*0.000,]$id,]
  # Subset data
  env.tmp <- env0[env0$group == k, ]
  otu_all2 <- otu_all[,colnames(otu_all) %in% env.tmp$sample.id]
  
  ID.tmp0<-ID0[ID0$subtype %in% rownames(otu_all),]
  
  otu <- t(otu_all2)[, (colSums(t(otu_all2)) / sum(colSums(t(otu_all2)))) >= 0.0001]
  otu2 <- otu_all2[colnames(otu),]
  otu2$otusum <- rowSums(otu2)  
  otu2$abu <- otu2$otusum/sum(otu2$otusum)
  spman.r0 = corr.test(otu, use="pairwise",method="spearman",adjust="fdr", alpha = .05, ci=FALSE)
  
  r_value<-spman.r0$r #提取相关系数
  p_value<-spman.r0$p
  r_value[p_value>0.01|abs(r_value)<0.6] = 0
  r_value1<<-r_value
  Cor<-as.matrix(spman.r0$r)  
  Cor[lower.tri(Cor)] <- NA
  diag(Cor) <- NA
  Cor.df <- reshape2::melt(Cor, na.rm = TRUE) 
  
  P0<-as.matrix(spman.r0$p)
  P0[lower.tri(P0)] <- NA
  diag(P0) <- NA
  P.df <- reshape2::melt(P0, na.rm = TRUE)
  
  df <- merge(Cor.df, P.df, by=c("Var1","Var2"))
  names(df) <- c("row","col","Cor","p")
  
  da.tmp <- df.sig <- filter(df, abs(Cor) >= 0.6, p < 0.01)
  
  da.g1 <- data.frame(da.tmp)
  #da.g1 <- na.omit(da.g1)
  da.g2 <<- da.g1
  g <- graph_from_data_frame(da.g1, 
                             directed=FALSE,
                             vertices = NULL)
  g1 <<- g
}


bacteria.network("Forest")
da.g1_forest<-da.g2
r_forest<-r_value1
r_upper <- r_forest
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_forest<-g1

bacteria.network("Cropland")
da.g1_cropland<-da.g2
r_cropland<-r_value1
r_upper <- r_cropland
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_cropland<-g1

bacteria.network("Industry")
da.g1_industry<-da.g2
r_industry<-r_value1
r_upper <- r_industry
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_industry<-g1

bacteria.network("Urban")
da.g1_urban<-da.g2
r_urban<-r_value1
r_upper <- r_urban
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_urban<-g1


source("/Users/yangy/Documents/E/city/bacteria/bacteria1/code/info.centrality.R")#自定义函数
cormatrix2<- r_forest
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_forest <- env0  %>% filter(group == "Forest")

otu1_forest <- otu1[group2_forest$sample.id,]
otu.rel <- sweep(otu1_forest, 1, rowSums(otu1_forest), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_forest) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_forest<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                          rbind(Weighted.simu, Unweighted.simu),
                          weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("Forest", 40)#根据自己的处理修改treat名称
)



###计算Cropland
cormatrix2<- r_cropland
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_cropland <- env0  %>% filter(group == "Cropland")

otu1_cropland <- otu1[group2_cropland$sample.id,]
otu.rel <- sweep(otu1_cropland, 1, rowSums(otu1_cropland), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_cropland) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_cropland<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("cropland", 40)#根据自己的处理修改treat名称
)


cormatrix2<- r_industry
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_industry <- env0  %>% filter(group == "Industry")

otu1_industry <- otu1[group2_industry$sample.id,]
otu.rel <- sweep(otu1_industry, 1, rowSums(otu1_industry), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_industry) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_industry<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("industry", 40)#根据自己的处理修改treat名称
)



cormatrix2<- r_urban
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/mge/results/normalized_level1_percell.csv",sep=",",header=T,row.names=1)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_urban <- env0  %>% filter(group == "Urban")

otu1_urban <- otu1[group2_urban$sample.id,]
otu.rel <- sweep(otu1_urban, 1, rowSums(otu1_urban), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_urban) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_urban<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                         rbind(Weighted.simu, Unweighted.simu),
                         weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("urban", 40)#根据自己的处理修改treat名称
)
df_all <- rbind(dat1_forest, dat1_cropland, dat1_industry, dat1_urban)

df_all$treat <- dplyr::recode(df_all$treat,
                              cropland = "Cropland",
                              industry = "Industry",
                              urban    = "Urban")


#生成非加权网络的结果图
p_Robustness_unweighted <- ggplot(df_all[df_all$weighted == "unweighted",], 
                                  aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +  # 注意 theme_light 要放前面，避免覆盖自定义 theme
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),   # 加粗 x 轴刻度
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),  # 加粗 x 轴标题
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "Mobilome-Robustness(Unweighted)")
p_Robustness_unweighted 

#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/mge_network_鲁棒性unweighted.pdf",p_Robustness_unweighted,width = 7, height = 6)


##生成加权网络的结果图
p_Robustness_weighted <- ggplot(df_all[df_all$weighted == "weighted",], 
                                aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "Mobilome-Robustness(Weighted)")


p_Robustness_weighted
#ggsave("/Users/yangy/Documents/D/city/ARG/figure/network/mge_network_鲁棒性weighted.pdf",p_Robustness_weighted,width = 7, height = 6)


##ARG network
library(GO.db)
library(vegan)
library(permute)        
library(lattice)
library(WGCNA)
library(igraph)
library(WGCNA)
library(psych)
library(reshape2)
library(igraph)
library(Hmisc)
library(GUniFrac)
library(sciplot)
library(ggpmisc)
library(indicspecies)
library(BiocManager)
library(patchwork)
library(agricolae)
library(ggplot2)
library(ggraph)
library(colorRamps)
rm(list=ls())

arg<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",header=T,row.names = 1) 
arg<-t(as.data.frame(arg))
rs.d<-arg[,colSums(arg)>0]
group<-read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv", header=TRUE)
argid<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv", header=TRUE)


rs_forest<-rs.d[group$group%in%c("Forest"),]
rs_forest1<-rs_forest[,colSums(rs_forest)>0]
rs_forest1<-as.data.frame(t(rs_forest1))
rs_forest1 <- t(rs_forest1)[, (colSums(t(rs_forest1)) / sum(colSums(t(rs_forest1)))) >= 0.0001]
trs_forest1<-t(rs_forest1)
trs_forest1<-data.frame(row.names(trs_forest1),trs_forest1)
names(trs_forest1)[1]<-c("gene")
forest<-merge(trs_forest1,argid,by="gene")
forest_order<-forest[order(forest$Type),]
forest_order1<-forest_order[-29]
forest_order2<-data.frame(forest_order1)
#write.csv(forest_order2,file = "forest_order3.csv")
rownames(forest_order2)<-forest_order2[,1]
forest_order3<-forest_order2[,-1]



rs_cropland<-rs.d[group$group%in%c("Cropland"),]
rs_cropland1<-rs_cropland[,colSums(rs_cropland)>0]
rs_cropland1<-as.data.frame(t(rs_cropland1))
rs_cropland1 <- t(rs_cropland1)[, (colSums(t(rs_cropland1)) / sum(colSums(t(rs_cropland1)))) >= 0.0001]
trs_cropland1<-t(rs_cropland1)
trs_cropland1<-data.frame(row.names(trs_cropland1),trs_cropland1)
names(trs_cropland1)[1]<-c("gene")
cropland<-merge(trs_cropland1,argid,by="gene")
cropland_order<-cropland[order(cropland$Type),]
cropland_order1<-cropland_order[-29]
cropland_order2<-data.frame(cropland_order1)
#write.csv(cropland_order2,file = "cropland_order3.csv")
rownames(cropland_order2)<-cropland_order2[,1]
cropland_order3<-cropland_order2[,-1]


rs_industry<-rs.d[group$group%in%c("Industry"),]
rs_industry1<-rs_industry[,colSums(rs_industry)>0]
rs_industry1<-as.data.frame(t(rs_industry1))
rs_industry1 <- t(rs_industry1)[, (colSums(t(rs_industry1)) / sum(colSums(t(rs_industry1)))) >= 0.0001]
trs_industry1<-t(rs_industry1)
trs_industry1<-data.frame(row.names(trs_industry1),trs_industry1)
names(trs_industry1)[1]<-c("gene")
industry<-merge(trs_industry1,argid,by="gene")
industry_order<-industry[order(industry$Type),]
industry_order1<-industry_order[-29]
industry_order2<-data.frame(industry_order1)
#write.csv(industry_order2,file = "industry_order3.csv")
rownames(industry_order2)<-industry_order2[,1]
industry_order3<-industry_order2[,-1]

rs_urban<-rs.d[group$group%in%c("Urban"),]
rs_urban1<-rs_urban[,colSums(rs_urban)>0]
rs_urban1<-as.data.frame(t(rs_urban1))
rs_urban1 <- t(rs_urban1)[, (colSums(t(rs_urban1)) / sum(colSums(t(rs_urban1)))) >= 0.0001]
trs_urban1<-t(rs_urban1)
trs_urban1<-data.frame(row.names(trs_urban1),trs_urban1)
names(trs_urban1)[1]<-c("gene")
urban<-merge(trs_urban1,argid,by="gene")
urban_order<-urban[order(urban$Type),]
urban_order1<-urban_order[-29]
urban_order2<-data.frame(urban_order1)
#write.csv(urban_order2,file = "urban_order3.csv")
rownames(urban_order2)<-urban_order2[,1]
urban_order3<-urban_order2[,-1]



############### forest
df1<-forest_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_forest<-igraph
#igraph<-igraph_forest
tax = read.csv('/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv',  header=T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))


E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####
par(mfrow=c(1,4))

data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv",header=T,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","grey","#7506F7","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_forest<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_forest<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                    diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)



forest <- plot(
  igraph,
  main = "Forest\n Nodes:  522\n Edges: 9830\n Average degree:  37.7",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)


##cropland
df1<-cropland_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_cropland<-igraph
#igraph<-igraph_cropland
tax = read.csv('/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv',  header=T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight


#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0,"grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv",header=T,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","grey","#7506F7","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)


num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_cropland<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_cropland<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)



cropland <- plot(
  igraph,
  main = "Cropland\n Nodes:  664\n Edges: 8369\n Average degree:  25.2",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)


##industry

df1<-industry_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_industry<-igraph 

tax = read.csv('/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv',  header=T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))
#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv",header=T,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","grey","#7506F7","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_industry<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_industry<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)


industry <- plot(
  igraph,
  main = "Industry\n Nodes:  575\n Edges: 2632\n Average degree:  9.2",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)


##Urban
df1<-urban_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_urban<-igraph
tax = read.csv('/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv',  header=T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0,"grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)



####
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv",header=T,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","grey","#7506F7","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_urban <- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_urban<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                   diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)



urban <- plot(
  igraph,
  main = "Urban\n Nodes:  755\n Edges: 19749\n Average degree:  52.3",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)



#save.image(file = "/Users/yangy/Documents/D/test/arg0.01_0.6/arg_all_0.01——0.6网络图.RData")

###########################################################################property_network

proall<-rbind(prold_forest,prold_cropland,prold_industry,prold_urban)
proall_network<-data.frame(proall)
proall_network$group<-c('Forest','Cropland','Industry','Urban')

##################################################################################Contribution to Degree

node_degree_df_forest$group<-rep("Forest",nrow(node_degree_df_forest))
node_degree_df_cropland$group<-rep("Cropland",nrow(node_degree_df_cropland))
node_degree_df_industry$group<-rep("Industry",nrow(node_degree_df_industry))
node_degree_df_urban$group<-rep("Urban",nrow(node_degree_df_urban))

degree_all<-rbind(node_degree_df_forest,node_degree_df_cropland,node_degree_df_industry,node_degree_df_urban)
degree_all$group<-factor(degree_all$group,levels = c('Forest','Cropland','Industry','Urban'))
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID_network.csv",header=T,row.names = 1)
degree_all$Type <- data2[degree_all$node, "Type"]


col11<- c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","darkgreen","#FFFF05","black","grey")

degree_all$Type <- factor(degree_all$Type, 
                          levels = c(setdiff(unique(degree_all$Type), "Others"), "Others"))
degree_all$Type <- gsub("Macrolide-lincosamide-streptogramin", "MLS", degree_all$Type)

f_degree_contri <- ggplot(degree_all, aes(x=group, y=degree, fill=Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x=" ", y = "Contribution to Degree") +
  labs(fill="  ") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 15, face="bold"),
    panel.spacing = unit(0, "lines"),
    legend.title = element_text(colour="black", size=15, face="bold"),
    legend.text = element_text(colour="black", size=15, face="bold"),
    axis.text = element_text(size=15, face="bold", colour="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 调整对齐
    axis.title = element_text(size=15, face="bold")
  ) +
  scale_fill_manual(values = col11)
f_degree_contri

#ggsave("/Users/yangy/Documents/D/test/arg0.01_0.6/contribution_to_degree_arg_all_0.1_0.06.pdf",f_degree_contri,heigh=6,width=8)

#pie
#正负连接饼图
#forest
library(ggplot2)
library(ggsci)

# 示例数据
pie_forest <- data.frame(
  links = c("positive", "negative"),
  value = c(6660, 3170)
)
library(dplyr)

pie_forest <- pie_forest %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_forest <- ggplot(pie_forest, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_forest

#cropland

# 示例数据
pie_cropland <- data.frame(
  links = c("positive", "negative"),
  value = c(7607, 762)
)

pie_cropland <- pie_cropland %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_cropland <- ggplot(pie_cropland, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_cropland



#industrty
pie_industry <- data.frame(
  links = c("positive", "negative"),
  value = c( 2614, 18)
)

pie_industry <- pie_industry %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_industry <- ggplot(pie_industry, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_industry


#urban
pie_urban <- data.frame(
  links = c("positive", "negative"),
  value = c(18757, 992)
)

pie_urban <- pie_urban %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_urban <- ggplot(pie_urban, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_urban
#ARG鲁棒性
rm(list=ls())


bacteria.network<-function(k){
  otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",sep=",",header=T,row.names=1)
  otu_all <- otu_all[, order(colnames(otu_all))]
  group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
  env0 <- group2[order(group2$sample.id), ]
  identical(group2$sample.id, colnames(otu_all))
  
  ID0 <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/all_arg_ID0.csv", check.names = F, header = T)
  # Set node color
  color <- c(
    "Aminoglycoside"                = "#1f77b4",  # 蓝
    "Beta_lactam"                   = "#ff7f0e",  # 橙
    "Macrolide-lincosamide-streptogramin" = "#2ca02c",  # 绿
    "Chloramphenicol"               = "#d62728",  # 红
    "Bacitracin"                    = "#9467bd",  # 紫
    "Multidrug"                     = "#8c564b",  # 棕
    "Fosfomycin"                    = "#e377c2",  # 粉
    "Tetracenomycin_C"              = "#7f7f7f",  # 灰
    "Tetracycline"                  = "#bcbd22",  # 黄绿
    "Other_peptide_antibiotics"     = "#17becf",  # 青
    "Trimethoprim"                  = "#aec7e8",  # 浅蓝
    "Polymyxin"                     = "#ffbb78",  # 浅橙
    "Vancomycin"                    = "#98df8a",  # 浅绿
    "Sulfonamide"                   = "#ff9896",  # 浅红
    "Quinolone"                     = "#c5b0d5",  # 浅紫
    "Florfenicol"                   = "#c49c94",  # 浅棕
    "Rifamycin"                     = "#f7b6d2",  # 淡粉
    "Novobiocin"                    = "#dbdb8d",  # 浅黄
    "Fusidic_acid"                  = "#9edae5",  # 浅青
    "Streptothricin"                = "#393b79",  # 深蓝
    "Antibacterial_fatty_acid"      = "#637939",  # 深绿
    "Defensin"                      = "#8c6d31",  # 深棕
    "Mupirocin"                     = "#843c39",  # 深红
    "Bleomycin"                     = "#7b4173",  # 深紫
    "Bicyclomycin"                  = "#5254a3",  # 靛蓝
    "Pleuromutilin_tiamulin"        = "#6b6ecf",  # 紫蓝
    "Factumycin"                    = "#9c9ede",  # 浅紫蓝
    "Puromycin"                     = "#17becf"   # 青绿
  )
  test <- as.data.frame(rowSums(otu_all))
  test$id <- rownames(test)
  test2 <- test[order(test$`rowSums(otu_all)`,decreasing = T),]
  otu_all <- otu_all[test2[test2$`rowSums(otu_all)` >= sum(test2$`rowSums(otu_all)`)*0.000,]$id,]
  # Subset data
  env.tmp <- env0[env0$group == k, ]
  otu_all2 <- otu_all[,colnames(otu_all) %in% env.tmp$sample.id]
  
  ID.tmp0<-ID0[ID0$gene %in% rownames(otu_all),]
  
  otu <- t(otu_all2)[, (colSums(t(otu_all2)) / sum(colSums(t(otu_all2)))) >= 0.0001]
  otu2 <- otu_all2[colnames(otu),]
  otu2$otusum <- rowSums(otu2)  
  otu2$abu <- otu2$otusum/sum(otu2$otusum)
  spman.r0 = corr.test(otu, use="pairwise",method="spearman",adjust="fdr", alpha = .05, ci=FALSE)
  
  r_value<-spman.r0$r #提取相关系数
  p_value<-spman.r0$p
  r_value[p_value>0.05|abs(r_value)<0.6] = 0
  r_value1<<-r_value
  Cor<-as.matrix(spman.r0$r)  
  Cor[lower.tri(Cor)] <- NA
  diag(Cor) <- NA
  Cor.df <- reshape2::melt(Cor, na.rm = TRUE) 
  
  P0<-as.matrix(spman.r0$p)
  P0[lower.tri(P0)] <- NA
  diag(P0) <- NA
  P.df <- reshape2::melt(P0, na.rm = TRUE)
  
  df <- merge(Cor.df, P.df, by=c("Var1","Var2"))
  names(df) <- c("row","col","Cor","p")
  
  da.tmp <- df.sig <- filter(df, abs(Cor) >= 0.6, p < 0.05)
  
  da.g1 <- data.frame(da.tmp)
  #da.g1 <- na.omit(da.g1)
  da.g2 <<- da.g1
  g <- graph_from_data_frame(da.g1, 
                             directed=FALSE,
                             vertices = NULL)
  g1 <<- g
  fc <- cluster_fast_greedy(g) 
  print(modularity(fc))  # >0.4 indicate modular structures
  modularity(fc,membership(fc))
  print(sizes(fc))
  
  # whole network properties
  net.metrics <- matrix(nrow=22,ncol=1)
  net.metrics[1,1] <- length(E(g))##length(curve_multiple(funteria. bacteria.archea.g))
  net.metrics[2,1] <- length(V(g))
  #net.fun.metrics[3,1] <- connectedness(as.matrix(as_adjacency_matrix(igraph, type = c("both"))))
  net.metrics[4,1] <- edge_density(g,loops=FALSE)
  net.metrics[5,1] <- mean(igraph::degree(g))
  net.metrics[6,1] <- average.path.length(g) 
  net.metrics[7,1] <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)
  net.metrics[8,1] <- edge_connectivity(g)
  net.metrics[9,1] <- vertex_connectivity(g)
  net.metrics[10,1] <- cohesion(g)
  net.metrics[11,1] <- transitivity(g)
  net.metrics[12,1] <- no.clusters(g)
  net.metrics[13,1] <- sum(graph.knn(g)$knn[!is.na(graph.knn(g)$knn)])/length(V(g))
  net.metrics[14,1] <- centralization.closeness(g)$centralization
  net.metrics[15,1] <- centralization.betweenness(g)$centralization
  net.metrics[16,1] <- centralization.degree(g)$centralization
  net.metrics[17,1] <- assortativity.degree(g)
  fc <- cluster_fast_greedy(g)
  net.metrics[18,1] <- modularity(g,membership(fc))
  net.metrics[19,1] <- mean(igraph::betweenness(g))
  net.metrics[20,1] <- diameter(g, directed = FALSE, unconnected = TRUE)
  net.metrics[21,1] <- sum(E(g)$Cor>0)  # number of positive correlation
  net.metrics[22,1] <- sum(E(g)$Cor<0)  # number of negative correlation
  
  rownames(net.metrics) <- c("num.edges","num.vertices","connectedness","connectance","average.degree","average.path.length","diameter",
                             "edge_connectivity","vertex.connectivity","cohesion","clustering.coefficient","no.clusters","avenei","centralization.closeness",
                             "centralization.betweenness","centralization.degree","assortativity.degree","modularity","average.betweenness","w.diameter","no.positive","no.negative")
  colnames(net.metrics) <- k
  #write.csv(net.bac.metrics,"16S.AU50.all.net.metricsHEL.csv")
  print(net.metrics)
  
  net.metrics <- as.data.frame(net.metrics)
  #net.metrics$name <- rownames(net.metrics)
  net.metrics1 <<- net.metrics
  
  # set vertices size
  v.size <- otu2[V(g)$name, ] 
  v.size1 = log(v.size$abu * 100000) 
  V(g)$size =  abs(v.size1)
  
  E.color <-  E(g)$Cor
  E.color = ifelse(E.color > 0, "red",ifelse(E.color < 0, "blue", "grey"))
  E(g)$color = as.character(E.color) 
  
  # color by Vertices/node with Type
  V(g)$Type <- ID0$Type[match(V(g)$name, ID0$gene)]
  #V(g)$color <- ID0$color[match(V(g)$name, ID0$OTU.ID)]
  
  E(g)$width = abs(E(g)$Cor)
  
  # title=paste("AMF",k,sep="-")
  name <- c("Forest" = "Forest", "Cropland" = "Cropland","Industry"="Industry","Urban"="Urban")
  
  set.seed(123)
  p1 <<-  ggraph(g, layout = "nicely") + 
    geom_edge_bend(aes(width = abs(E(g)$Cor), color = E(g)$color), strength = 0.2) + 
    scale_edge_colour_manual(values= c("blue"="blue", "red"="red")) +
    scale_edge_width_continuous(range = c(0.5, 1.2)) +
    geom_node_point(aes(color=factor(V(g)$Type), size = V(g)$size)) +
    scale_colour_manual(name="Type",values=color) +
    scale_size_continuous(range = c(1,5)) +
    labs(title=sprintf("ARGs-%s",name[k]),
         subtitle =  sprintf("Vertice = %d; Edge = %d; Average degree = %g", 
                             net.metrics["num.vertices", ], 
                             net.metrics["num.edges", ], 
                             round(net.metrics["average.degree", ], 2))) +
    guides(
      edge_colour = "none",
      edge_width = "none",
      size = "none",
      color = guide_legend(ncol = 1)   # ✅ 这里设置图例为一列
    ) +
    theme_void() +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
      plot.title = element_text(size = 20, face = "bold",hjust=0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 12, face = "bold",hjust=0.5, color="red", margin = margin(b = 10)), 
      legend.title = element_text(colour = "black", size = 15,face = "bold"),
      legend.text = element_text(colour = "black", size = 12,face = "bold")
    )
  p1
}



bacteria.network("Forest")
net.metrics <- net.metrics1
net.metrics_forest <- net.metrics1
da.g1_forest<-da.g2
r_forest<-r_value1
r_upper <- r_forest
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_forest<-g1
p_arg_forest <- p1

bacteria.network("Cropland")
net.metrics <- net.metrics1
net.metrics_cropland <- net.metrics1
net.metrics$Cropland <- net.metrics1
da.g1_cropland<-da.g2
r_cropland<-r_value1
r_upper <- r_cropland
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_cropland<-g1
p_arg_cropland <- p1

bacteria.network("Industry")
net.metrics <- net.metrics1
net.metrics_industry <- net.metrics1
net.metrics$industry <- net.metrics1
da.g1_industry<-da.g2
r_industry<-r_value1
r_upper <- r_industry
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_industry<-g1
p_arg_industry <- p1

bacteria.network("Urban")
net.metrics <- net.metrics1
net.metrics_urban <- net.metrics1
net.metrics$urban <- net.metrics1
da.g1_urban<-da.g2
r_urban<-r_value1
r_upper <- r_urban
r_upper[lower.tri(r_upper, diag = TRUE)] <- NA
library(reshape2)
r_df <- melt(r_upper, na.rm = TRUE)
colnames(r_df) <- c("Source", "Target", "Correlation")
r_df <- subset(r_df, Correlation != 0)
g_urban<-g1
p_arg_urban <- p1

arg_combined_plot <- (p_arg_forest | p_arg_cropland ) /
  (p_arg_industry | p_arg_urban) 
arg_combined_plot


save.image(file = "/Users/yangy/Documents/E/city/metagenome/arg/results/arg_co-occurence-network.RData")


rm(list=ls())
source("/Users/yangy/Documents/E/city/bacteria/bacteria1/code/info.centrality.R")#自定义函数
load("/Users/yangy/Documents/E/city/metagenome/arg/results/arg_co-occurence-network.RData")
cormatrix2<- r_forest
diag(cormatrix2)<- 0 # 去除自连接



#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv",sep=",",header=T,row.names=1)

otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_forest <- env0  %>% filter(group == "Forest")

otu1_forest <- otu1[group2_forest$sample.id,]
otu.rel <- sweep(otu1_forest, 1, rowSums(otu1_forest), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_forest) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_forest<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                          rbind(Weighted.simu, Unweighted.simu),
                          weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("Forest", 40)#根据自己的处理修改treat名称
)



###计算Cropland
cormatrix2<- r_cropland
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv", row.names = 1, check.names = F, header = T)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_cropland <- env0  %>% filter(group == "Cropland")

otu1_cropland <- otu1[group2_cropland$sample.id,]
otu.rel <- sweep(otu1_cropland, 1, rowSums(otu1_cropland), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_cropland) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_cropland<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("cropland", 40)#根据自己的处理修改treat名称
)


cormatrix2<- r_industry
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv", row.names = 1, check.names = F, header = T)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_industry <- env0  %>% filter(group == "Industry")

otu1_industry <- otu1[group2_industry$sample.id,]
otu.rel <- sweep(otu1_industry, 1, rowSums(otu1_industry), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_industry) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_industry<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("industry", 40)#根据自己的处理修改treat名称
)



cormatrix2<- r_urban
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/normalized_percell.csv", row.names = 1, check.names = F, header = T)
otu1<-t(otu_all)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_urban <- env0  %>% filter(group == "Urban")

otu1_urban <- otu1[group2_urban$sample.id,]
otu.rel <- sweep(otu1_urban, 1, rowSums(otu1_urban), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_urban) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_urban<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                         rbind(Weighted.simu, Unweighted.simu),
                         weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("urban", 40)#根据自己的处理修改treat名称
)
df_all <- rbind(dat1_forest, dat1_cropland, dat1_industry, dat1_urban)

df_all$treat <- dplyr::recode(df_all$treat,
                              cropland = "Cropland",
                              industry = "Industry",
                              urban    = "Urban")


#生成非加权网络的结果图
p_Robustness_unweighted <- ggplot(df_all[df_all$weighted == "unweighted",], 
                                  aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +  # 注意 theme_light 要放前面，避免覆盖自定义 theme
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),   # 加粗 x 轴刻度
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),  # 加粗 x 轴标题
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "ARGs-Robustness(Unweighted)")
p_Robustness_unweighted 

##生成加权网络的结果图
p_Robustness_weighted <- ggplot(df_all[df_all$weighted == "weighted",], 
                                aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "ARGs-Robustness(Weighted)")


p_Robustness_weighted
#108bac
library(GO.db)
library(vegan)
library(permute)        
library(lattice)
library(WGCNA)
library(igraph)
library(WGCNA)
library(psych)
library(reshape2)
library(igraph)
library(Hmisc)
library(GUniFrac)
library(sciplot)
library(ggpmisc)
library(indicspecies)
library(BiocManager)
library(patchwork)
library(agricolae)
library(ggplot2)
library(ggraph)
library(colorRamps)
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
rs.d<-otu[,colSums(otu)>0]
otuid <-read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F, header = T)
group<-group0[,-1]

rs_forest<-rs.d[group$group%in%c("Forest"),]
rs_forest1<-rs_forest[,colSums(rs_forest)>0]
rs_forest1<-as.data.frame(t(rs_forest1))
rs_forest1 <- t(rs_forest1)[, (colSums(t(rs_forest1)) / sum(colSums(t(rs_forest1)))) >= 0.0001]
trs_forest1<-t(rs_forest1)
trs_forest1<-data.frame(row.names(trs_forest1),trs_forest1)
names(trs_forest1)[1]<-c("OTU.ID")
forest<-merge(trs_forest1,otuid,by="OTU.ID")
forest_order<-forest[order(forest$Phylum),]
forest_order1<-forest_order[-29]
forest_order2<-data.frame(forest_order1)
#write.csv(forest_order2,file = "forest_order3.csv")
rownames(forest_order2)<-forest_order2[,1]
forest_order3<-forest_order2[,-1]


rs_cropland<-rs.d[group$group%in%c("Cropland"),]
rs_cropland1<-rs_cropland[,colSums(rs_cropland)>0]
rs_cropland1<-as.data.frame(t(rs_cropland1))
rs_cropland1 <- t(rs_cropland1)[, (colSums(t(rs_cropland1)) / sum(colSums(t(rs_cropland1)))) >= 0.0001]
trs_cropland1<-t(rs_cropland1)
trs_cropland1<-data.frame(row.names(trs_cropland1),trs_cropland1)
names(trs_cropland1)[1]<-c("OTU.ID")
cropland<-merge(trs_cropland1,otuid,by="OTU.ID")
cropland_order<-cropland[order(cropland$Phylum),]
cropland_order1<-cropland_order[-29]
cropland_order2<-data.frame(cropland_order1)
#write.csv(cropland_order2,file = "cropland_order3.csv")
rownames(cropland_order2)<-cropland_order2[,1]
cropland_order3<-cropland_order2[,-1]


rs_industry<-rs.d[group$group%in%c("Industry"),]
rs_industry1<-rs_industry[,colSums(rs_industry)>0]
rs_industry1<-as.data.frame(t(rs_industry1))
rs_industry1 <- t(rs_industry1)[, (colSums(t(rs_industry1)) / sum(colSums(t(rs_industry1)))) >= 0.0001]
trs_industry1<-t(rs_industry1)
trs_industry1<-data.frame(row.names(trs_industry1),trs_industry1)
names(trs_industry1)[1]<-c("OTU.ID")
industry<-merge(trs_industry1,otuid,by="OTU.ID")
industry_order<-industry[order(industry$Phylum),]
industry_order1<-industry_order[-29]
industry_order2<-data.frame(industry_order1)
#write.csv(industry_order2,file = "industry_order3.csv")
rownames(industry_order2)<-industry_order2[,1]
industry_order3<-industry_order2[,-1]

rs_urban<-rs.d[group$group%in%c("Urban"),]
rs_urban1<-rs_urban[,colSums(rs_urban)>0]
rs_urban1<-as.data.frame(t(rs_urban1))
rs_urban1 <- t(rs_urban1)[, (colSums(t(rs_urban1)) / sum(colSums(t(rs_urban1)))) >= 0.0001]
trs_urban1<-t(rs_urban1)
trs_urban1<-data.frame(row.names(trs_urban1),trs_urban1)
names(trs_urban1)[1]<-c("OTU.ID")
urban<-merge(trs_urban1,otuid,by="OTU.ID")
urban_order<-urban[order(urban$Phylum),]
urban_order1<-urban_order[-29]
urban_order2<-data.frame(urban_order1)
#write.csv(urban_order2,file = "urban_order3.csv")
rownames(urban_order2)<-urban_order2[,1]
urban_order3<-urban_order2[,-1]

############### forest
par(mfrow=c(1,4))

df1<-forest_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_forest<-igraph
#igraph<-igraph_forest
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F, header = T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Phylum = tax$Phylum

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Phylum) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))
E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)




data2=read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F,row.names = 1, header = T)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","grey","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)
num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_forest<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_forest<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                    diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)


forest <- plot(
  igraph,
  main = "Urban\n Nodes:  1698\n Edges: 35900\n Average degree:  42.3",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)

#cropland
df1<-cropland_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_cropland<-igraph
#igraph<-igraph_cropland
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F, header = T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Phylum = tax$Phylum

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Phylum) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))
E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)



####

data2=read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F,row.names = 1, header = T)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","grey","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)
num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_cropland<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_cropland<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)

cropland<-plot(igraph,main="Cropland\n Nodes:  30101\n Edges: 636\n Average degree: 31.3",vertex.frame.color=NA,
               edge.lty=1,edge.curved=TRUE,edge.width=0.01,margin=c(0,0,0,0),
               vertex.label=NA, layout=layout_in_circle,vertex.size=5)


cropland <- plot(
  igraph,
  main = "Cropland\n Nodes:  1963\n Edges: 30737\n Average degree:  31.3",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)


##industry
df1<-industry_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_industry<-igraph
#igraph<-igraph_industry
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F, header = T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Phylum = tax$Phylum

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Phylum) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight
E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))
#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####

data2=read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F,row.names = 1, header = T)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","grey","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_industry<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_industry<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)

industry <- plot(
  igraph,
  main = "Industry\n Nodes:  1708\n Edges: 6938\n Average degree:  8.1",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)




##urban
df1<-urban_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.01|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 
igraph_urban<-igraph
#igraph<-igraph_urban
tax = read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F, header = T)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Phylum = tax$Phylum

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Phylum) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E.color <- ifelse(E.color > 0,"#A9484B",
                  ifelse(E.color < 0, "grey","black"))
E(igraph)$color = as.character(E.color)
set.seed(123)




data2=read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F,row.names = 1, header = T)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","grey","darkgreen","#FFFF05","black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_urban<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_urban<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                   diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)


urban <- plot(
  igraph,
  main = "Urban\n Nodes:  1566\n Edges: 7444\n Average degree:  9.5",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 0.001,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.4)  # 设置透明度（0到1之间）
)


proall<-rbind(prold_forest,prold_cropland,prold_industry,prold_urban)
proall_network<-data.frame(proall)
proall_network$group<-c('Forest','Cropland','Industry','Urban')

write.csv(proall_network,"/Users/yangy/Documents/E/city/metagenome/arg/results/108bac/108.bac_0.01_network_property.csv")



node_degree_df_forest$group<-rep("Forest",nrow(node_degree_df_forest))
node_degree_df_cropland$group<-rep("Cropland",nrow(node_degree_df_cropland))
node_degree_df_industry$group<-rep("Industry",nrow(node_degree_df_industry))
node_degree_df_urban$group<-rep("Urban",nrow(node_degree_df_urban))

degree_all<-rbind(node_degree_df_forest,node_degree_df_cropland,node_degree_df_industry,node_degree_df_urban)
degree_all$group<-factor(degree_all$group,levels = c('Forest','Cropland','Industry','Urban'))
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/bacteria.OTU.ID1_network.csv", check.names = F,row.names = 1, header = T)
degree_all$Phylum <- data2[degree_all$node, "Phylum"]


col11<- c("red","orange" ,"green","cyan", "#0000FF","#F706DB","#7506F7","darkgreen","#FFFF05","black","grey")

degree_all$Phylum <- factor(degree_all$Phylum, 
                            levels = c(setdiff(unique(degree_all$Phylum), "Others"), "Others"))
f_degree_contri <- ggplot(degree_all, aes(x=group, y=degree, fill=Phylum)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x=" ", y = "Contribution to Degree") +
  labs(fill="  ") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 15, face="bold"),
    panel.spacing = unit(0, "lines"),
    legend.title = element_text(colour="black", size=15, face="bold"),
    legend.text = element_text(colour="black", size=15, face="bold"),
    axis.text = element_text(size=15, face="bold", colour="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 调整对齐
    axis.title = element_text(size=15, face="bold")
  ) +
  scale_fill_manual(values = col11)
f_degree_contri


ggsave("D:/test/10bbac0.6/contribution_to_degree_bac_all.pdf",f_degree_contri,heigh=6,width=8)
ggsave("D:/test/10bbac0.6/contribution_to_degree_bac_all.jpg",f_degree_contri,heigh=6,width=8)

save.image(file = "D:/test/10bbac0.6/bac_all_网络图展示版本_dong.RData")


#正负连接饼图
#forest
library(ggplot2)
library(ggsci)

# 示例数据
pie_forest <- data.frame(
  links = c("positive", "negative"),
  value = c(33119, 2781)
)
library(dplyr)

pie_forest <- pie_forest %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_forest <- ggplot(pie_forest, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_forest
ggsave("D:/test/10bbac0.6/positive_negative_link_forest_108bac.pdf",p1_forest,heigh=5,width=6)

#cropland

# 示例数据
pie_cropland <- data.frame(
  links = c("positive", "negative"),
  value = c(30101, 636)
)

pie_cropland <- pie_cropland %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_cropland <- ggplot(pie_cropland, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_cropland

ggsave("D:/test/10bbac0.6/positive_negative_link_cropland_108bac.pdf",p1_cropland,heigh=5,width=6)


#industrty
pie_industry <- data.frame(
  links = c("positive", "negative"),
  value = c(6938, 0)
)

pie_industry <- pie_industry %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_industry <- ggplot(pie_industry, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_industry

ggsave("D:/test/10bbac0.6/positive_negative_link_industry_108bac.pdf",p1_industry,heigh=5,width=6)

#urban
pie_urban <- data.frame(
  links = c("positive", "negative"),
  value = c(7375, 69)
)

pie_urban <- pie_urban %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_urban <- ggplot(pie_urban, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_urban

ggsave("D:/test/10bbac0.6/positive_negative_link_urban_108bac.pdf",p1_urban,heigh=5,width=6)


##bacteria 鲁棒性
rm(list=ls())
source("/Users/yangy/Documents/E/city/bacteria/bacteria1/code/info.centrality.R")#自定义函数
load("/Users/yangy/Documents/E/city/metagenome/arg/results/108bac_co-occurence-network.RData")
cormatrix2<- r_forest
diag(cormatrix2)<- 0 # 去除自连接



#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]


otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
otu_all0 <- otu_all0[order(rownames(otu_all0)),]

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]

otu1 <- otu[order(rownames(otu)), ]

library(dplyr)
group2 <-read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T,row.names = 1)
env0 <- group2[order(group2$sample.id), ]
group2_forest <- env0  %>% filter(group == "Forest")

otu1_forest <- otu1[group2_forest$sample.id,]

otu.rel <- sweep(otu1_forest, 1, rowSums(otu1_forest), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_forest) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_forest<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                          rbind(Weighted.simu, Unweighted.simu),
                          weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("Forest", 40)#根据自己的处理修改treat名称
)



###计算Cropland
cormatrix2<- r_cropland
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]


otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
otu_all0 <- otu_all0[order(rownames(otu_all0)),]

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]

otu1 <- otu[order(rownames(otu)), ]

library(dplyr)
group2 <-read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T,row.names = 1)
env0 <- group2[order(group2$sample.id), ]
group2_cropland <- env0  %>% filter(group == "Cropland")

otu1_cropland <- otu1[group2_cropland$sample.id,]

otu.rel <- sweep(otu1_cropland, 1, rowSums(otu1_cropland), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_cropland) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_cropland<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("Cropland", 40)#根据自己的处理修改treat名称
)






cormatrix2<- r_industry
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]


otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
otu_all0 <- otu_all0[order(rownames(otu_all0)),]

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]

otu1 <- otu[order(rownames(otu)), ]

library(dplyr)
group2 <-read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T,row.names = 1)
env0 <- group2[order(group2$sample.id), ]
group2_industry <- env0  %>% filter(group == "Industry")

otu1_industry <- otu1[group2_industry$sample.id,]

otu.rel <- sweep(otu1_industry, 1, rowSums(otu1_industry), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_industry) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_industry<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("Industry", 40)#根据自己的处理修改treat名称
)



cormatrix2<- r_urban
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]


otu_all0 <- read.csv("/Users/yangy/Documents/E/city/bacteria/bacteria_flattening.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))
otu_all0 <- otu_all0[order(rownames(otu_all0)),]

# Read metadata file
group0 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T)
group0 <- na.omit(group0)
group <- group0[order(group0$sample.id),]
otu_all0$sample.id <- rownames(otu_all0)
otu_sub <- otu_all0[otu_all0$sample.id %in% group$sample.id,]
otu<-otu_sub[,1:6927]

otu1 <- otu[order(rownames(otu)), ]

library(dplyr)
group2 <-read.csv("/Users/yangy/Documents/E/city/metagenome/group1.csv",header = T,row.names = 1)
env0 <- group2[order(group2$sample.id), ]
group2_urban <- env0  %>% filter(group == "Urban")

otu1_urban <- otu1[group2_urban$sample.id,]

otu.rel <- sweep(otu1_urban, 1, rowSums(otu1_urban), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_urban) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_urban<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                         rbind(Weighted.simu, Unweighted.simu),
                         weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("urban", 40)#根据自己的处理修改treat名称
)

df_all <- rbind(dat1_forest, dat1_cropland, dat1_industry, dat1_urban)

df_all$treat <- dplyr::recode(df_all$treat,
                              cropland = "Cropland",
                              industry = "Industry",
                              urban    = "Urban")


#生成非加权网络的结果图
p_Robustness_unweighted <- ggplot(df_all[df_all$weighted == "unweighted",], 
                                  aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +  # 注意 theme_light 要放前面，避免覆盖自定义 theme
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),   # 加粗 x 轴刻度
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),  # 加粗 x 轴标题
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "Bacteria-Robustness(Unweighted)")

p_Robustness_unweighted 

#ggsave("/Users/yangy/Documents/E/city/metagenome/arg/results/108bac/108bac_network_鲁棒性unweighted.pdf",p_Robustness_unweighted,width = 7, height = 6)


p_Robustness_weighted <- ggplot(df_all[df_all$weighted == "weighted",], 
                                aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "Bacteria-Robustness(Weighted)")


p_Robustness_weighted

##high risk arg
library(vegan)
library(psych)
library(igraph)
library(ggplot2)
library(colorRamps)
library(ape)
library(stringr)
library(reshape2)
library(splitstackshape)
library(emmeans)
library(graphics)
library(phyloseq)
library(dplyr)
library(ggraph)
library(pals)
library(RColorBrewer)
library(patchwork)
library(scales) 
rm(list=ls())


arg <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
arg<-arg[,-1]
arg<-t(as.data.frame(arg))
rs.d<-arg[,colSums(arg)>0]
group<-read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv", header=TRUE)
argid<-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE)


rs_forest<-rs.d[group$group%in%c("Forest"),]
rs_forest1<-rs_forest[,colSums(rs_forest)>0]
rs_forest1<-as.data.frame(t(rs_forest1))
rs_forest1 <- t(rs_forest1)[, (colSums(t(rs_forest1)) / sum(colSums(t(rs_forest1)))) >= 0.0001]
trs_forest1<-t(rs_forest1)
trs_forest1<-data.frame(row.names(trs_forest1),trs_forest1)
names(trs_forest1)[1]<-c("subtype")
forest<-merge(trs_forest1,argid,by="subtype")
forest_order<-forest[order(forest$Type),]
forest_order1<-forest_order[-29]
forest_order2<-data.frame(forest_order1)
#write.csv(forest_order2,file = "forest_order3.csv")
rownames(forest_order2)<-forest_order2[,1]
forest_order3<-forest_order2[,-1]


rs_cropland<-rs.d[group$group%in%c("Cropland"),]
rs_cropland1<-rs_cropland[,colSums(rs_cropland)>0]
rs_cropland1<-as.data.frame(t(rs_cropland1))
rs_cropland1 <- t(rs_cropland1)[, (colSums(t(rs_cropland1)) / sum(colSums(t(rs_cropland1)))) >= 0.0001]
trs_cropland1<-t(rs_cropland1)
trs_cropland1<-data.frame(row.names(trs_cropland1),trs_cropland1)
names(trs_cropland1)[1]<-c("subtype")
cropland<-merge(trs_cropland1,argid,by="subtype")
cropland_order<-cropland[order(cropland$Type),]
cropland_order1<-cropland_order[-29]
cropland_order2<-data.frame(cropland_order1)
#write.csv(cropland_order2,file = "cropland_order3.csv")
rownames(cropland_order2)<-cropland_order2[,1]
cropland_order3<-cropland_order2[,-1]


rs_industry<-rs.d[group$group%in%c("Industry"),]
rs_industry1<-rs_industry[,colSums(rs_industry)>0]
rs_industry1<-as.data.frame(t(rs_industry1))
rs_industry1 <- t(rs_industry1)[, (colSums(t(rs_industry1)) / sum(colSums(t(rs_industry1)))) >= 0.0001]
trs_industry1<-t(rs_industry1)
trs_industry1<-data.frame(row.names(trs_industry1),trs_industry1)
names(trs_industry1)[1]<-c("subtype")
industry<-merge(trs_industry1,argid,by="subtype")
industry_order<-industry[order(industry$Type),]
industry_order1<-industry_order[-29]
industry_order2<-data.frame(industry_order1)
#write.csv(industry_order2,file = "industry_order3.csv")
rownames(industry_order2)<-industry_order2[,1]
industry_order3<-industry_order2[,-1]


rs_urban<-rs.d[group$group%in%c("Urban"),]
rs_urban1<-rs_urban[,colSums(rs_urban)>0]
rs_urban1<-as.data.frame(t(rs_urban1))
rs_urban1 <- t(rs_urban1)[, (colSums(t(rs_urban1)) / sum(colSums(t(rs_urban1)))) >= 0.0001]
trs_urban1<-t(rs_urban1)
trs_urban1<-data.frame(row.names(trs_urban1),trs_urban1)
names(trs_urban1)[1]<-c("subtype")
urban<-merge(trs_urban1,argid,by="subtype")
urban_order<-urban[order(urban$Type),]
urban_order1<-urban_order[-29]
urban_order2<-data.frame(urban_order1)
#write.csv(urban_order2,file = "urban_order3.csv")
rownames(urban_order2)<-urban_order2[,1]
urban_order3<-urban_order2[,-1]


############### forest
df1<-forest_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.05|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 

tax = read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))
E.color <- rep("grey", length(E.color))
E(igraph)$color = as.character(E.color)
set.seed(123)



####
par(mfrow=c(1,4))

data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","green","skyblue" ,"#0000FF","pink","black")
("Aminoglycoside" ="red"
  "Beta_lactam"  ="green"
  "Chloramphenicol" ="skyblue"
  "Multidrug"   ="#0000FF"   
  "Trimethoprim"  ="pink"
  "Vancomycin"  ="black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)


num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_forest<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_forest<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                    diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)


forest <- plot(
  igraph,
  main = "Forest\n Nodes:  17\n Edges: 16\n Average degree:  1.9",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 2,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.6)  # 设置透明度（0到1之间）
)



############### cropland
df1<-cropland_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.05|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 

tax = read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight
E.color <- rep("grey", length(E.color))

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)




data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
#Aminoglycoside="red"
"Chloramphenicol"="skyblue"
"Forfenicol"="brown"
"Multidrug"="#0000FF"
"Quinolone"="#00a8f3"
"Trimethoprim"   ="pink"
color = c("red","brown","skyblue","#0000FF","#00a8f3","pink")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)

num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_cropland<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_cropland<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)

cropland <- plot(
  igraph,
  main = "Cropland\n Nodes:  13\n Edges: 13\n Average degree:  2",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 2,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.6)  # 设置透明度（0到1之间）
)

############### industry
df1<-industry_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.05|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 

tax = read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight
E.color <- rep("grey", length(E.color))

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)


data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","skyblue","#0000FF")



levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)





num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_industry<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_industry<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                      diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)

industry <- plot(
  igraph,
  main = "Industry\n Nodes:  6\n Edges: 3\n Average degree:  1",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 2,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.6)  # 设置透明度（0到1之间）
)

#urban
df1<-urban_order3
df_corr = corr.test(t(df1), use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
df_corr_r = df_corr$r
df_corr_p = df_corr$p
df_corr_r[df_corr_p>0.05|abs(df_corr_r)<0.6] = 0


igraph <- graph_from_adjacency_matrix(df_corr_r,mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs <- V(igraph )[igraph::degree(igraph ) == 0]
igraph  = delete.vertices(igraph , bad.vs)
igraph 

tax = read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE)
rownames(tax)<-tax[,1]
tax = tax[as.character(V(igraph)$name), ]

V(igraph)$Type = tax$Type

node_list = data.frame(
  label = names(V(igraph)),
  genus=V(igraph)$Type) 

igraph.weight = E(igraph)$weight
E(igraph)$weight = NA
sum(igraph.weight>0)
sum(igraph.weight<0)
average.degree = mean(igraph::degree(igraph))
average.degree
E.color = igraph.weight
E.color <- rep("grey", length(E.color))

#E.color <- ifelse(E.color > 0,"grey",
#ifelse(E.color < 0, "#A9484B","black"))

E(igraph)$color = as.character(E.color)
set.seed(123)



####

data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE,row.names = 1)

igraph_col = data2[V(igraph)$name,]
igraph_col2=as.factor(igraph_col)
levels(igraph_col2)

set.seed(123)
color = c("red","orange","green","skyblue" ,"brown","#0000FF","#00a8f3","pink","black")
("Aminoglycoside" ="red"
  "Bacitracin"="orange",
  "Beta_lactam"  ="green"
  "Chloramphenicol" ="skyblue"
  "Forfenicol"="brown"
  "Multidrug"   ="#0000FF"   
  "Quinolone"  ="#00a8f3"
  "Trimethoprim"  ="pink"
  "Vancomycin"  ="black")

levels(igraph_col2) = color
V(igraph)$color = as.character(igraph_col2)


num.edges = length(E(igraph)) 
num.edges
connectance = edge_density(igraph,loops=FALSE)
connectance
average.degree = mean(igraph::degree(igraph))
average.degree
average.path.length = average.path.length(igraph)
average.path.length
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter
edge.connectivity = edge_connectivity(igraph)
edge.connectivity
clustering.coefficient = transitivity(igraph)
clustering.coefficient
no.clusters = no.clusters(igraph)
no.clusters
centralization.betweenness = centralization.betweenness(igraph)$centralization
centralization.betweenness
centralization.degree = centralization.degree(igraph)$centralization
centralization.degree
fc = cluster_fast_greedy(igraph,weights =NULL)
modularity = modularity(igraph,membership(fc))
modularity
positive<-sum(igraph.weight>0)
negative<-sum(igraph.weight<0)
num.vertices = length(V(igraph))
num.vertices
node_degree <- igraph::degree(igraph)
node_degree_df_urban<- data.frame(
  node = names(node_degree),
  degree = node_degree,
  row.names = names(node_degree)  # 设置行名
)
prold_urban<-cbind(modularity,num.edges,num.vertices,connectance,average.degree,average.path.length,
                   diameter,edge.connectivity,clustering.coefficient,centralization.betweenness,centralization.degree,positive,negative)


urban <- plot(
  igraph,
  main = "Urban\n Nodes:  37\n Edges: 49\n Average degree:  2.6",
  vertex.frame.color = NA,
  edge.lty = 1,
  edge.curved = TRUE,
  edge.width = 2,
  margin = c(0,0,0,0),
  vertex.label = NA,
  layout = layout_in_circle,
  vertex.size = 5,
  edge.color = adjustcolor("grey", alpha.f = 0.6)  # 设置透明度（0到1之间）
)


###########################################################################property_network

proall<-rbind(prold_forest,prold_cropland,prold_industry,prold_urban)
proall_network<-data.frame(proall)
proall_network$group<-c('Forest','Cropland','Industry','Urban')

write.csv(proall_network,"/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg/high_risk_arg_property.csv")


##################################################################################Contribution to Degree

node_degree_df_forest$group<-rep("Forest",nrow(node_degree_df_forest))
node_degree_df_cropland$group<-rep("Cropland",nrow(node_degree_df_cropland))
node_degree_df_industry$group<-rep("Industry",nrow(node_degree_df_industry))
node_degree_df_urban$group<-rep("Urban",nrow(node_degree_df_urban))

degree_all<-rbind(node_degree_df_forest,node_degree_df_cropland,node_degree_df_industry,node_degree_df_urban)
degree_all$group<-factor(degree_all$group,levels = c('Forest','Cropland','Industry','Urban'))
data2=read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg_ID.csv", header=TRUE,row.names = 1)
degree_all$Type <- data2[degree_all$node, "Type"]
unique(degree_all$Type)
("Aminoglycoside" ="red"
  "Bacitracin"="orange",
  "Beta_lactam"  ="green"
  "Chloramphenicol" ="skyblue"
  "Forfenicol"="brown"
  "Multidrug"   ="#0000FF"   
  "Quinolone"  ="#00a8f3"
  "Trimethoprim"  ="pink"
  "Vancomycin"  ="black")


col11<- c("red","green","skyblue","#0000FF","pink","black","brown","#00a8f3","orange")

degree_all$Type <- factor(degree_all$Type, 
                          levels = c(setdiff(unique(degree_all$Type), "Others"), "Others"))
f_degree_contri <- ggplot(degree_all, aes(x=group, y=degree, fill=Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x=" ", y = "Contribution to Degree") +
  labs(fill="  ") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 15, face="bold"),
    panel.spacing = unit(0, "lines"),
    legend.title = element_text(colour="black", size=15, face="bold"),
    legend.text = element_text(colour="black", size=15, face="bold"),
    axis.text = element_text(size=15, face="bold", colour="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 调整对齐
    axis.title = element_text(size=15, face="bold")
  ) +
  scale_fill_manual(values = col11)
f_degree_contri





##正负相关饼图

#pie
#正负连接饼图
#forest
library(ggplot2)
library(ggsci)

# 示例数据
pie_forest <- data.frame(
  links = c("positive", "negative"),
  value = c(16, 0)
)
library(dplyr)

pie_forest <- pie_forest %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_forest <- ggplot(pie_forest, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_forest
ggsave("D:/test/high0.01_0.6/positive_negative_link_forest_high_arg_0.6_0.01.pdf",p1_forest,heigh=5,width=6)

#cropland

# 示例数据
pie_cropland <- data.frame(
  links = c("positive", "negative"),
  value = c(13, 0)
)

pie_cropland <- pie_cropland %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_cropland <- ggplot(pie_cropland, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_cropland



#industrty
pie_industry <- data.frame(
  links = c("positive", "negative"),
  value = c( 3, 0)
)

pie_industry <- pie_industry %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_industry <- ggplot(pie_industry, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_industry


#urban
pie_urban <- data.frame(
  links = c("positive", "negative"),
  value = c(49, 0)
)

pie_urban <- pie_urban %>%
  arrange(desc(links)) %>%
  mutate(
    fraction = value / sum(value),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

# 绘制饼图并加上数值标签
p1_urban <- ggplot(pie_urban, aes(x = "", y = value, fill = links)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(y = label_pos * sum(value), label = value), color = "black", size = 10) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  labs(fill = " ", x = NULL, y = NULL, title = " ") +
  scale_fill_manual(values = c("positive" = "#A9484B", "negative" = "grey"))

p1_urban


#high risk 鲁棒性
rm(list=ls())

source("/Users/yangy/Documents/E/city/bacteria/bacteria1/code/info.centrality.R")#自定义函数
load("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_co-occurence-network.RData")
cormatrix2<- r_forest
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
otu1<-otu_all[,-1]
otu1<-t(otu1)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_forest <- env0  %>% filter(group == "Forest")

otu1_forest <- otu1[group2_forest$sample.id,]
otu.rel <- sweep(otu1_forest, 1, rowSums(otu1_forest), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_forest) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_forest<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                          rbind(Weighted.simu, Unweighted.simu),
                          weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("Forest", 40)#根据自己的处理修改treat名称
)



###计算Cropland
cormatrix2<- r_cropland
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
otu1<-otu_all[,-1]
otu1<-t(otu1)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_cropland <- env0  %>% filter(group == "Cropland")

otu1_cropland <- otu1[group2_cropland$sample.id,]
otu.rel <- sweep(otu1_cropland, 1, rowSums(otu1_cropland), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_cropland) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_cropland<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("cropland", 40)#根据自己的处理修改treat名称
)


cormatrix2<- r_industry
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
otu1<-otu_all[,-1]
otu1<-t(otu1)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_industry <- env0  %>% filter(group == "Industry")

otu1_industry <- otu1[group2_industry$sample.id,]
otu.rel <- sweep(otu1_industry, 1, rowSums(otu1_industry), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_industry) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_industry<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                            rbind(Weighted.simu, Unweighted.simu),
                            weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("industry", 40)#根据自己的处理修改treat名称
)



cormatrix2<- r_urban
diag(cormatrix2)<- 0 # 去除自连接

#计算网络连接数量和节点数量
sum(abs(cormatrix2)> 0) / 2 # 连接数
sum(colSums(abs(cormatrix2))> 0) # 至少有一个连接的节点数
#提取连接网络的矩阵
network.raw<- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) >0]

otu_all <-read.csv("/Users/yangy/Documents/E/city/metagenome/arg/results/high_risk_arg.csv", row.names = 1, check.names = F, header = T)
otu1<-otu_all[,-1]
otu1<-t(otu1)
otu1 <- otu1[order(rownames(otu1)), ]
library(dplyr)
group2 <- read.csv("/Users/yangy/Documents/E/city/metagenome/group.csv",sep=",",header=T)
env0 <- group2[order(group2$sample.id), ]
group2_urban <- env0  %>% filter(group == "Urban")

otu1_urban <- otu1[group2_urban$sample.id,]
otu.rel <- sweep(otu1_urban, 1, rowSums(otu1_urban), "/")

sp.ra <- colMeans(otu.rel)

#sp.ra<- colMeans(otu1_urban) / 2892 # 每种物种的相对丰度
otus_in_network <- colnames(cormatrix2)[colSums(abs(cormatrix2)) > 0]
sp.ra2 <- sp.ra[otus_in_network]
sum(row.names(network.raw)== names(sp.ra2)) # 检查是否匹配

#鲁棒性模拟函数#随机移除部分物种后计算剩余物种比例
rand.remov.once<- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) { 
  id.rm <- sample(1:nrow(netRaw),
                  round(nrow(netRaw) * rm.percent)) 
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0;  net.Raw[, id.rm] <- 0  # 移除物种 
  if (abundance.weighted) { 
    net.stength <- net.Raw * sp.ra 
  } else { 
    net.stength <- net.Raw 
  } 
  sp.meanInteration <- colMeans(net.stength) 
  id.rm2 <- which(sp.meanInteration <=0)  # 移除无连接的物种 
  remain.percent <- (nrow(netRaw) -length(id.rm2)) / nrow(netRaw) 
  remain.percent
}


#鲁棒性模拟
rmsimu<- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm =100) 
{ t(sapply(rm.p.list, function(x) { 
  remains <- sapply(1:nperm, function(i) { 
    rand.remov.once(netRaw = netRaw,
                    rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted) }) 
  remain.mean <- mean(remains) 
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm) 
  result <- c(remain.mean, remain.sd,remain.se) 
  names(result) <-c("remain.mean", "remain.sd", "remain.se") 
  result }))
}
#加权和非加权模拟
Weighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = TRUE, nperm = 100)

Unweighted.simu<- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra =sp.ra2, abundance.weighted = FALSE, nperm = 100)

#整合结果
dat1_urban<- data.frame( Proportion.removed = rep(seq(0.05, 1, by =0.05), 2), 
                         rbind(Weighted.simu, Unweighted.simu),
                         weighted = rep(c("weighted","unweighted"), each = 20), treat = rep("urban", 40)#根据自己的处理修改treat名称
)
df_all <- rbind(dat1_forest, dat1_cropland, dat1_industry, dat1_urban)

df_all$treat <- dplyr::recode(df_all$treat,
                              cropland = "Cropland",
                              industry = "Industry",
                              urban    = "Urban")


#生成非加权网络的结果图
p_Robustness_unweighted <- ggplot(df_all[df_all$weighted == "unweighted",], 
                                  aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +  # 注意 theme_light 要放前面，避免覆盖自定义 theme
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),   # 加粗 x 轴刻度
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),  # 加粗 x 轴标题
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "High risk ARGs-Robustness(Unweighted)")
p_Robustness_unweighted 



##生成加权网络的结果图
p_Robustness_weighted <- ggplot(df_all[df_all$weighted == "weighted",], 
                                aes(x = Proportion.removed, y = remain.mean, group = treat, color = treat)) + 
  geom_line(size=1) + 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 1) + 
  scale_color_manual(
    name = "Ecosystem",
    limits = c("Forest", "Cropland", "Industry", "Urban"),
    values = c("Forest" = "#6cb509", 
               "Cropland" = "#fdb218", 
               "Industry" = "#5e9ceb", 
               "Urban" = "#ff7f27")
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(colour = "black", size = 12, face = "bold"),
    strip.text = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),  # 黑色panel边框
  ) +
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+
  labs(title = "High risk ARGs-Robustness(Weighted)")


p_Robustness_weighted


