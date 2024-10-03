#######################Fig1
library(ggplot2)
library(reshape2)
library(vegan)
library(ape)
library(picante) 
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson ָ??
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

otu <- read.csv("w1.csv",sep = ",",row.names = 1)
otu <-t(otu)
alpha_all <- alpha(otu, base = 2)

p1<-ggplot(data=da, aes(x = State, y = Richness,
                        color=State)) +
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA)+
  scale_color_manual(limits=c("Healthy","Diseased"), 
                     values=c("#203378","#e3bf2b"))+
  stat_compare_means(method = "wilcox.test",paired = F, 
                     comparisons=list(c("Healthy","Diseased")))+
  geom_jitter(size=0.5)+
scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  facet_wrap(~Compartment, ncol=4)+
  theme_bw(base_size = 16) +
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=20,vjust = 1.5),
        axis.title.y=element_text(colour='black', size=20,vjust = 1.5),#face = "bold"),
        axis.text.y=element_text(colour='black',size=18),#face = "bold"),
        axis.text.x=element_text(colour='black',size=18,angle = 90,hjust = 1),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = 'right')+
  labs(y='Bacteria richness')
p1


########Beta

otu <- read.csv("wr.csv",sep = ",",row.names = 1)
otu <- data.frame(t(otu))
write.csv(otu,"t_wr.csv")
group<-read.csv("Group.csv",sep = ",")
otu <- read.csv("nosoil_wr.csv",sep = ",",row.names = 1)
otu <- data.frame(t(otu))
group<-read.csv("Group.csv",sep = ",")
library(dplyr)
group <- group %>% 
  filter(State != "Soil")
distance <- vegdist(otu, method = 'bray')

mod <- betadisper(d = distance, group = group$State, type = 'centroid')
mod
set.seed(123)
p_disease<-permutest(mod, pairwise = TRUE, permutations = 999)
otuput <- data.frame(p_disease$tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'N.Perm', 'Pr (>F)')
disper<-data.frame(mod$distances)
write.csv(disper,"Whole_disper.csv")
write.csv(otuput,"Wholedisper_test.csv")

pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')

summary(pcoa)

point <- data.frame(pcoa$point)  

pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

sample_site <- data.frame({pcoa$point})[1:2]
sample_site$SampleID <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

sample_site <- merge(sample_site, group, by = 'SampleID', all.x = TRUE)

library(dplyr)
sample_site <- sample_site %>% 
  filter(State != "Soil")
sample_site$Compartment <- factor(sample_site$Compartment, levels = c('Rhizosphere',"Root" ,"Stem" ,"Seed" ))
sample_site$State<- factor(sample_site$State, levels = c('Healthy','Disease'))

cbbPalette <- c("#203378","#e3bf2b","#e31a1c","#C5199E","#FF5900", "#0D0DBB","#C5199E",
                "#65472F","#A4C9CC", "#DD7694", "#04376E")
cbbPalette <- c("#203378","#6989b9","#a9c1a3","#e3bf2b")



p1<- ggplot(sample_site, aes(PCoA1, PCoA2, group = Compartment)) +
  geom_point(aes(color = Compartment,shape =State), size = 3, alpha = 1)+
  scale_color_manual(values = cbbPalette) +
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  #???????????????ᡢͼ???ĸ?ʽ
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank())+
  labs(x = paste('PCoA 1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA 2: ', round(100 * pcoa_eig[2], 2), '%'))+
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18,vjust = 1.5),
        axis.title.y=element_text(colour='black', size=18,vjust = 1.5),
        axis.text.y=element_text(colour='black',size=15),
        axis.text.x=element_text(colour = "black",size = 15,
                                 hjust = 1,vjust = 0.5))+
  stat_ellipse(aes(color = Compartment), level = 0.95, linetype = 6, show.legend = FALSE)+
  annotate("text",-0.30,0.52,label = "Site:R2=0.1076** 
     Compartment:R2=0.078**  
     Site:R2=0.0013**",size=5)
p1

######################DDRs
library(geosphere)
site <- read.csv('Dis.csv', row.names = 1)
library(dplyr)
site<-site %>% 
  filter(State == "Soil")
site1<-site[-1:-9]
site1<-scale(site1)
site_dis<- as.matrix(vegan::vegdist(site1, method = 'euclidean'))

rownames(site_dis) <- rownames(site)
colnames(site_dis) <- rownames(site)

site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)
head(site_dis)

spe <- read.csv('SoilSoil.csv', row.names = 1)
spe <- data.frame(t(spe))

comm_sim <- 1 - as.matrix(vegan::vegdist(spe, method = 'bray'))

diag(comm_sim) <- 0  
comm_sim[upper.tri(comm_sim)] <- 0  
comm_sim <- reshape2::melt(comm_sim)
comm_sim <- subset(comm_sim, value != 0)
head(comm_sim)

comm_dis <- merge(comm_sim, site_dis, by = c('Var1', 'Var2'))
names(comm_dis) <- c('site1', 'site2', 'comm_sim', 'site_dis')
comm_dis$com <- 'Soil' 
comm_dis$sta <- 'Soil' 
so<-comm_dis
###################
###############Rhh
library(geosphere)
site <- read.csv('Dis.csv', row.names = 1)
library(dplyr)
site<-site %>% 
  filter(gro3 == "RhizosphereHealthy")

site1<-site[-1:-9]
site1<-scale(site1)
site_dis<- as.matrix(vegan::vegdist(site1, method = 'euclidean'))

rownames(site_dis) <- rownames(site)
colnames(site_dis) <- rownames(site)

site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)
head(site_dis)
#####################################plot
ge <- subset(net, Envrionment == "Geography")

ge$Compartment<-factor(ge$Compartment,levels = c("Soil","Rhizosphere","Root","Stem","Seed"))

p1<-ggplot(ge, aes(x = Compartment, y = Slop, fill = State)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#e3bf2b","#203378","#C75B5B")) +
  theme_bw(base_size = 16) +
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=20,vjust = 1.5),
        axis.title.y=element_text(colour='black', size=20,vjust = 1.5),#face = "bold"),
        axis.text.y=element_text(colour='black',size=18),#face = "bold"),
        axis.text.x=element_text(colour='black',size=18,angle = 90),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = 'right')+
  labs(y='Geography DDR Slop')

#################Fig2
############stem microbiome composition
library(reshape2)
library(vegan)
library(ggplot2)
library(RColorBrewer)
tax<- read.table("tax.csv",header=T,sep=",")
otu <-read.csv("wr.csv",row.names = 1)
tax$phylum<-factor(tax$Phylum)
tax$class<-factor(tax$Class)
tax$order<-factor(tax$Order)
tax$family<-factor(tax$Family)
tax$genus<-factor(tax$Genus)
c<- aggregate(otu, by=list(tax$class), sum) 
o<- aggregate(otu, by=list(tax$order), sum) 
f<- aggregate(otu, by=list(tax$family), sum) 
g<-aggregate(otu, by=list(tax$genus), sum) 

write.csv(p,"phy.csv",row.names = F)
write.csv(c,"class.csv",row.names = F)
write.csv(o,"order.csv",row.names = F)
write.csv(f,"family.csv",row.names = F)
write.csv(g,"genus.csv",row.names = F)

taxon=tax
c<-read.csv("phy.csv",row.names = 1)

top_10<- names(head(sort(rowSums(c), decreasing = T), 10)) 
top_10
top_10<-c("Proteobacteria"  , "Actinobacteriota" ,"Firmicutes"    ,   "Bacteroidota"    
          , "Chloroflexi"  ,    "Patescibacteria" , "Acidobacteriota" , "Gemmatimonadota" 
          , "Myxococcota"  ,    "Unclassified" ) 
taxon$Phylum[!(taxon$Phylum)%in%top_10]<-"Others" 
c_top<- aggregate(otu, by=list(taxon$Phylum), sum)
rownames(c_top)<-c_top[,1]
c_top=c_top[,-1]
c_top=c_top[order(-rowSums(c_top)),]
c_top1<-c_top
c_top2<-t(c_top1)
c_top2 <- cbind(rownames(c_top2), c_top2)
rownames(c_top2) <- NULL
colnames(c_top2)[1]<-c("SampleID")
group<-read.csv("Group.csv")
c_top4<-merge(group,c_top2,by="SampleID")
c_top5<-melt(c_top4,id=1:6)
colnames(c_top5)[7:8]<-c("Phylum","Abundance")
g<-melt(g)
write.csv(c_top4,"p_top10.csv",row.names = F)
data<-c_top5
library(dplyr)
data <- data %>% 
  filter(State != "Soil")
data$Abundance<-as.numeric(data$Abundance)

da1<-data%>% 
  select(Compartment,State,Phylum,Abundance) %>%
  group_by(Compartment,State,Phylum) %>% 
  summarise_all(mean)
dim(da1)
head(da1)

library(ggalluvial)
library(imputeTS)

cols <- c("#203378","#e3bf2b","#e31a1c","#1f78b4","#00B76D","#FF5900","#C5199E", "#C2855E",
          "#A4C9CC","#65472F", "#DD7694")

ggplot(data = da1,
            aes(x = State, y = Abundance, fill = reorder(Phylum,-Abundance), 
                stratum = reorder(Phylum,-Abundance),
                alluvium = reorder(Phylum,-Abundance))) +
  geom_alluvium()+
  geom_stratum(width=0.45, size=0.1) +
  geom_bar(aes(fill = reorder(Phylum,-Abundance)),stat='identity', width=0.7) +
  scale_y_continuous(expand=c(0, 0))+
  theme_bw() +
  facet_grid( ~ Compartment,scales = "fixed")+
  scale_fill_manual(values = cols,name="Phylum") +
  scale_color_manual(values = cols) +
  theme(legend.position = "right",
        panel.grid=element_blank(),
        panel.spacing.x = unit(0,units = "cm"),
        strip.background = element_rect(
          color="white", fill="white", 
          linetype="solid",size=0.8),
        strip.placement = "outside",
        axis.line.y.left = element_line(color="black",size=0.8),
        axis.line.x.bottom = element_line(color="black",size=0.8),
        strip.text.x = element_text(size = 14,face = "bold"),
        axis.text = element_text(face = "bold", 
                                 size = 12,color="black"),
        axis.title = element_text(face = "bold", 
                                  size = 14,colour = "black"),
        legend.title = element_text(face = "bold", 
                                    size =12,color="black"),
        legend.text = element_text(face = "bold", size =12,color="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.3),
  )+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18,vjust = 1.5),
        axis.title.y=element_text(colour='black', size=18,vjust = 1.5),
        axis.text.y=element_text(colour='black',size=15),
        axis.text.x=element_text(colour = "black",size = 15,
                                 hjust = 1,vjust = 1,angle = 45))+
  labs(x = "State",y= "Relative Abundance of Phylum (%)")

#############edge
###########Stem
counts <- read.csv("Stem_w.csv")
colnames(counts)[1] <- "ID"
metadata <- read.csv("Group_Stem.csv")
counts = data.frame(counts)
row_name = counts$ID
row.names(counts) <- row_name
counts = counts[,-1]
counts <- counts[which(rowSums(counts) > 20), ]
dgelist <- DGEList(counts = counts, group = metadata$State)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~metadata$State)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE) 
lrt <- glmLRT(fit)
aa<-lrt$table
aa <- cbind(ID = rownames(aa), aa)
rownames(aa) <- NULL
colnames(aa)[1] <- "ID"
aa$Compartment <- c("Stem")
aa3<-merge(aa,tax,by="ID")
write.csv(aa3, "Stem_result.csv",row.names = F)

#########################RF
fu<-read.csv("fusarium.csv",row.names = 1)
core<-as.data.frame(t(read.csv("Core14.csv",row.names = 1)))
core$SampleID<-rownames(core)
co<-merge(fu,core,by="SampleID")
co<-co[-1]
library(rfPermute)
set.seed(123)

otu_rfP <- rfPermute(Fusarium~., data = co, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 4)
otu_rfP$rf

importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = T)
importance_otu.scale
core<-importance_otu.scale[1]
names(core)[1]<-c("lncMSE")
core$tax<-c("Core")
nc<-read.csv("noCore3018.csv",row.names = 1)
nc <- nc[which(rowSums(nc) > 10), ] 
nco<-as.data.frame(t(nc))
nco$SampleID<-rownames(nco)
nco<-merge(fu,nco,by="SampleID")
nco<-nco[-1]
library(rfPermute)
set.seed(123)
otu_rfP <- rfPermute(Fusarium~., data = nco, importance = TRUE, ntree = 500, nrep = 10, num.cores = 4)
otu_rfP$rf
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = T)
importance_otu.scale
ncore<-importance_otu.scale[1]
names(ncore)[1]<-c("lncMSE")
ncore$tax<-c("Others")
ta<-rbind(core,ncore)

####################heatmap


mat <- as.data.frame(t(read.csv('heatmap.csv', row.names = 1)))
mat$SampleID<-rownames(mat)

ge<-read.csv('Group.csv')
library(dplyr)
ge <- ge %>% 
  filter(Compartment == "Stem")
ge<-ge[-2:-3]
ge<-ge[1:2]

aa<-merge(ge,mat,by="SampleID")

row.names(aa)<-aa$SampleID
aa<-aa[-1]
mat <-t(aa[-1])
mat <- as.matrix(mat)

mat <- mat - rowMeans(mat)  
mat <- mat / apply(mat, 1, max) * 2 - 1  

sample_group<-aa[1]
sample_group <- as.matrix(sample_group)
gene_anno <- as.matrix(read.csv('anno.csv', row.names = 1))


Heatmap(
  mat, name = 'ZRA', 
  col = colorRampPalette(c('grey30', 'grey50', 'grey80', '#FFA48C', '#FF1A09'))(100), 
  cluster_rows = T, cluster_columns = F,  
  row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10), 
  row_split = gene_anno, column_split = sample_group, 
  row_title_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 10), 

  right_annotation = HeatmapAnnotation(
    Class = gene_anno,  which = 'row', show_annotation_name = FALSE, 
    col = list(Class= c('Pantoea' = '#FFAD30', 'Klebsiella' = '#634FB8', 
                        'Enterobacter' = '#68AD30',"Lactococcus"="#f49128",
                        "Rahnella" ="#194a55","Stenotrophomonas"="#187c65",
                        "Escherichia"="#f26115","Pseudomonas"="#c29f62"))
  ),
  
  top_annotation = HeatmapAnnotation(
    Group = sample_group, which = 'column', show_annotation_name = FALSE, 
    col = list(Group = c('Healthy' = "#203378", 'Diseased' = "#e3bf2b"))
  )
)



#######################################SourceTracker
source('SourceTracker.r')
metadata <- read.table('gro6.txt',sep='\t',h=T,row.names=1,check=F,comment='')
otus <- read.table('w6.txt',sep='\t', header=T,row.names=1,check=F,skip=1,comment='')
otus<-round(otus);
otus <- t(as.matrix(otus));
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,];
metadata <- metadata[common.sample.ids,]
train.ix <- which(metadata$SourceSink=='source');
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env;
alpha1 <- alpha2 <- 0.001;
st <- sourcetracker(otus[train.ix,], envs[train.ix])
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)
a<-data.frame(results$proportions)
a

########################################Fig3
########################################automap
library(automap)
require(ggplot2)
require(RColorBrewer)
library(sp)
library(maptools)
library(rgdal)
library(mapproj)
library(sf)
library(raster)
library(reshape2)
library(tidyverse)
library(dplyr)



#####Geo: The geographic coordinates of each location 
#####value: The observed value of each location
#####df_grid: The prediction grid of the interpolation regions
china_map = readShapePoly("bou2_4p.shp")
china_map1<-fortify(china_map)
#Regions defined for each Polygons
china_map1$long<-round(china_map1$long,2)
china_map1$lat<-round(china_map1$lat,2)
NAME<-c("安徽省",XXXXXXXXX)
df_map<-china_map[as.vector(china_map@data$NAME)%in%NAME,]
df_grid1<-spsample(df_map, n =20461, "regular")
gridded(df_grid1)<-T
rich<-read.csv("b_a.fusarium",row.names = 1)
Geo<-read.csv("Dis.csv")
names(rich)[1]<-c("SampleID")
wh<-merge(Geo,rich,by="SampleID")
library(dplyr)
soil<- wh %>% 
 filter(State == "Soil")
        
soil<-soil[8:11]
soi<-melt(soil,id=4)
cbbPalette <- c("#203378","#6989b9","#a9c1a3","#e3bf2b")
library(ggpubr)
p1<-ggplot(soil, aes( x = log(valule),y  = Fusarium,
                        color = variable)) +
          geom_point( aes( color = variable), size=4)+  
          scale_color_manual(values = cbbPalette) +
          geom_smooth(method = 'lm', formula = y ~ x, se = T) + 
          stat_cor(data=soi, method = "spearman",size=6,label.y.npc = "top", label.x.npc = "center")+
          theme_bw()+
          theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
                axis.line = element_line(colour = "black"), 
                axis.title.x=element_text(colour='black', size=20,vjust = 1.5),
                axis.title.y=element_text(colour='black', size=20,vjust = 1.5),#face = "bold"),
                axis.text.y=element_text(colour='black',size=18),#face = "bold"),
                axis.text.x=element_text(colour = "black",size = 18,#face = "bold",
                                         angle = 0,hjust = 0.5,vjust = 0.5),
                strip.text = element_text(colour = "black",size = 18,face = "bold"),
                legend.background = element_rect(fill = "transparent"),
                legend.position = "top")

data<-soil
coordinates(data) =~ long+lat
kriging_result = autoKrige(Fusarium~1,input_data=data,new_data=df_grid1,data_variogram=data)
prediction_spdf = kriging_result$krige_output
a<-coordinates(prediction_spdf)

        
