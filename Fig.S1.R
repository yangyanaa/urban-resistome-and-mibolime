#Fig.S1
library(RColorBrewer)
library(cowplot)
library(colorRamps)
library(ggrepel)
library(ggspatial)
library(sf)
rm(list=ls())
group <- read.csv("/Users/yangy/Documents/E/city/bacteria/group.csv",sep=",",header=T,row.names = 1) 

unique_locations <- group %>%
  dplyr::group_by(Location) %>%
  dplyr::summarise(
    Latitude  = dplyr::first(Latitude),
    Longitude = dplyr::first(Longitude))


# 查看结果
print(unique_locations)
unique_locations<-as.data.frame(unique_locations)
group<-unique_locations[,c(2,3)]
rownames(group)<-unique_locations$Location

# 将 `group` 表转换为 sf 对象
site.CFB2 <- st_as_sf(group, coords = c("Longitude", "Latitude"), crs = 4326)
site.CFB2$sample_name<-row.names(group)

# 确保 sample_name 列按顺序排列
site.CFB2$Site <- factor(site.CFB2$sample_name, levels = c("Nanning","Guangzhou","Kunming","Fuzhou","Guiyang","Changsha","Hangzhou","Wuhan","Chengdu","Shanghai","Hefei","Nanjing","Xian","Zhengzhou","Lanzhou","Xining","Jinan","Taiyuan","Shijiazhuang","Yinchuan","Tianjin","Beijing","Hohhot","Shenyang","Urumqi","Changchun","Harbin"))

china <- st_read("/Users/yangy/Documents/E/city/map/China.json")
colors1 <- colorRampPalette(c("#0000FF", "#00FFFF"))(14)
colors2 <- colorRampPalette(c("#FFEB00", "#FF0202"))(13)

colors <- c(colors1, colors2)

mapCN<-ggplot() + 
  geom_sf(data = china,fill="NA",color="grey",size=0.5) + 
  #geom_sf(data = nine_line,color="black",size=0.5) + 
  scale_size(range = c(1,5))+
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering)+
  geom_sf(data=site.CFB2,aes(color=Site),size=5)+
  geom_text_repel(data = site.CFB2,aes(label=sample_name,geometry=geometry),stat = "sf_coordinates",color="blue",size=5)+
  #annotate("text",x=-2556175.2,y=-2387082,label="GS京(2022)1061号",size=3)+
  scale_colour_manual(values= colors,name="Site.latitude")+
  coord_sf(ylim = c(-2387082,1654989),xlim = c(-2556175.2,2816095),crs = "+proj=laea +lat_0=40 +lon_0=104")+
  theme_bw()+xlab("Longitude")+ylab("Latitude")+
  theme(plot.margin=unit(c(0,0,0,0),"mm"),
        axis.text = element_text(size = 18,colour = "black",face = "bold"),
        axis.title = element_text(size = 18,colour = "black",face = "bold"))

mapCN_nineline<-ggplot() + 
  geom_sf(data = china,fill="NA",color="grey",size=0.5) +
  scale_size(range = c(1,5))+
  annotation_scale(location = "br") +
  coord_sf(ylim = c(-4228017,-1877844),xlim = c(117131.4,2115095),crs="+proj=laea +lat_0=40 +lon_0=104")+
  theme_bw()+theme(aspect.ratio = 1.5,plot.margin=unit(c(0,0,0,0),"mm"))

ggdraw()+
  draw_plot(mapCN)+
  draw_plot(mapCN_nineline,x = 0.65, y = 0.19, width = 0.15, height = 0.230)
