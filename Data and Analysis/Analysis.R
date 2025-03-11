rm()

library(tidyverse)
library(readxl)

help()

#-----------Heavy metal------------
setwd("path")

data.raw <- read_excel('path/Heavy metal data.xls',sheet = "nosign")
data.raw <- column_to_rownames(data.raw, "CID2") 
data.env <- read_excel('path/Heavy metal data.xls',sheet = "nosign")

##PCA
data <- data.raw[,3:10]
indoor <- subset(data.raw, data.raw$Sample == "indoor")[,3:10]
outdoor <- subset(data.raw, data.raw$Sample == "outdoor")[,3:10]
Air <- subset(data.raw, data.raw$Sample == "Air-conditioning")[,3:10]
------
data <-as.matrix(scale(indoor))
rm1<-cor(data)
rs1<-eigen(rm1)
val <- rs1$values
St <- sqrt(abs(val))
Proportion <- round((val/sum(val))*100,digits = 2)
Cumulative_Proportion <- cumsum(Proportion)
Loading<-as.matrix(rs1$vectors)
PC <-data %*% Loading
-------
write.table(PC,"indoor.xls",sep = "\t",row.names = TRUE, fileEncoding = "UTF-8") 
write.table(Loading,"all_loading.xls",sep = "\t",row.names = TRUE, fileEncoding = "UTF-8")
print(Proportion)

##Permdisp
library(vegan) 
indoor <- subset(data.raw, data.raw$Sample == "indoor")
group <- read_excel('path/Risk.xls',sheet = "PLI")[1:67,c(1,12)]
df.betadisper <- betadisper(d = vegdist(data.raw[,-1:-2]), group = data.raw$Sample  , type = "centroid")
df.betadisper <- betadisper(d = vegdist(indoor[,-1:-2]), group = group$group  , type = "centroid")
anova(df.betadisper)

#Anosim
library(vegan) 
group <- read_excel('path/Risk.xls',sheet = "PLI")[1:67,c(1,12)]
df.dist=vegdist(data.raw[,-1:-2])
df.dist=vegdist(indoor[,-1:-2])
df.ano1=anosim(df.dist, data.raw$Sample, permutations = 999)
df.ano2=anosim(df.dist, group$group, permutations = 999)
summary(df.ano1)
summary(df.ano2)

##GI
data.GI <- data.env[,c(1:2)] 
data.GI$CdGI <- log2(data.env$Cd/(1.5*0.09))
data.GI$CoGI <- log2(data.env$Co/(1.5*10))
data.GI$CrGI <- log2(data.env$Cr/(1.5*58))
data.GI$CuGI <- log2(data.env$Cu/(1.5*20))
data.GI$MnGI <- log2(data.env$Mn/(1.5*521))
data.GI$NiGI <- log2(data.env$Ni/(1.5*25))
data.GI$PbGI <- log2(data.env$Pb/(1.5*19))
data.GI$VGI <- log2(data.env$V/(1.5*71))
write.table(data.GI,"GI.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")

##EF【ref:Sr】
data.EF <- data.env[,c(1:2)]
data.EF$Sr <- read_excel('path/Heavy metal data.xls',sheet = "Sr")
data.EF$CdEF <- (data.env$Cd/data.EF$Sr)/(0.09/271)
data.EF$CoEF <- (data.env$Co/data.EF$Sr)/(10/271)
data.EF$CrEF <- (data.env$Cr/data.EF$Sr)/(58/271)
data.EF$CuEF <- (data.env$Cu/data.EF$Sr)/(20/271)
data.EF$MnEF <- (data.env$Mn/data.EF$Sr)/(521/271)
data.EF$NiEF <- (data.env$Ni/data.EF$Sr)/(25/271)
data.EF$PbEF <- (data.env$Pb/data.EF$Sr)/(19/271)
data.EF$VEF <- (data.env$V/data.EF$Sr)/(71/271)
write.table(data.EF,"EF.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")

##PLI
data.PLI <- data.env[,c(1,3:11)]
data.PLI$CdPLI <- data.PLI$Cd/0.09
data.PLI$CoPLI <- data.PLI$Co/10
data.PLI$CrPLI <- data.PLI$Cr/58
data.PLI$CuPLI <- data.PLI$Cu/20
data.PLI$MnPLI <- data.PLI$Mn/521
data.PLI$NiPLI <- data.PLI$Ni/25
data.PLI$PbPLI <- data.PLI$Pb/19
data.PLI$VPLI <- data.PLI$V/71
data.PLI$PLI <- apply(data.PLI[, 11:18], 1, prod)^(1/8)
write.table(data.PLI,"PLI.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")

##RI  
data.RI <- data.env[,c(1,3:11)]
data.RI$CdEr <- (data.RI$Cd/0.09)*30
data.RI$CoEr <- (data.RI$Co/10)*5
data.RI$CrEr <- (data.RI$Cr/58)*2
data.RI$CuEr <- (data.RI$Cu/20)*5
data.RI$MnEr <- (data.RI$Mn/521)*1
data.RI$NiEr <- (data.RI$Ni/25)*5
data.RI$PbEr <- (data.RI$Pb/19)*5
data.RI$VEr <- (data.RI$V/71)*2
data.RI$RI <- apply(data.RI[, 11:18], 1, sum)
write.table(data.RI,"RI.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")

#-----------Metagenomics------------
library(vegan)
setwd("path")

data.raw <- read_excel('path/env.micbio.database.p.xls')
data.raw <- column_to_rownames(data.raw, "phylum") 
-------
data.raw <- read_excel('path/env.micbio.database.s.xls')
data.raw <- column_to_rownames(data.raw, "species")
data.raw <- t(data.raw)
-------
indoor <- subset(data.raw, data.raw$Sample == 1)
outdoor <- subset(data.raw, data.raw$Sample == 2)
Air <- subset(data.raw, data.raw$Sample == 3)

##Rank-abundance
options(scipen = 999)
all.otu <- read_excel('path/mic_database.xls')
all.otu <- column_to_rownames(all.otu, "ID") 
all.otu <- t(all.otu)
sp <- specaccum(all.otu, method = 'random')
tiff(filename="species_accumulation_curves.tiff",width = 11, height = 11, units = "cm", res = 600) 
plot(sp, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'lightblue')
boxplot(sp, col = 'yellow', add = TRUE, pch = '+')
dev.off()

install.packages("BiodiversityR")
library(BiodiversityR) 
otu_relative <- all.otu/rowSums(all.otu)
rank_dat <- data.frame()
for (i in rownames(all.otu)) {
  rank_dat_i <- data.frame(rankabundance(subset(all.otu, rownames(all.otu) == i), digits = 6))[1:2]
  rank_dat_i$sample <- i
  rank_dat <- rbind(rank_dat, rank_dat_i)
}
rank_dat <- subset(rank_dat, abundance != 0)

#Alpha analysis
data.raw <- data.raw[,c(3:204)]
indoor <- indoor[,c(3:204)]
indoor <- read_excel('path/Indoor&SHAP.xls',sheet = "indoor")
indoor <- column_to_rownames(indoor, "species")
indoor <- t(indoor)
shannon=diversity(data.raw,"shannon")
simpson=diversity(data.raw,"simpson")
shannon=diversity(indoor,"shannon")
simpson=diversity(indoor,"simpson")
write.table(shannon,"Shannon.xls",sep = "\t",row.names = TRUE, fileEncoding = "UTF-8")
write.table(simpson,"Simpson.xls",sep = "\t",row.names = TRUE, fileEncoding = "UTF-8")

##PCA
data <- data.raw[,c(3:204)]
data <-as.matrix(scale(data))
rm1<-cor(data)
missing_values <- colSums(is.na(rm1))
columns_to_delete <- which(missing_values == 201) 
rm1 <- subset(rm1, select = -columns_to_delete) 
data <- subset(data, select = -columns_to_delete) 
rm1 <- na.omit(rm1)  
rs1<-eigen(rm1)
val <- rs1$values
St <- sqrt(abs(val))
Proportion <- val/sum(val)
Cumulative_Proportion <- cumsum(Proportion)
Loading<-as.matrix(rs1$vectors)
PC <-data %*% Loading
write.table(PC,"PCA.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

##PCoA
data <- data.raw[,c(3:204)]
indoor <- indoor[,c(3:204)]
-------
indoor <- read_excel('path/Indoor&SHAP.xls',sheet = "indoor_s")
indoor <- column_to_rownames(indoor, "species")
indoor <- t(indoor)
-------
otu.distance <- vegdist(indoor,method='bray')
df.pcoa <- cmdscale(otu.distance,eig=TRUE,k=4)
pcoa <- df.pcoa$points[,1:4]
pc <- round(df.pcoa$eig/sum(df.pcoa$eig)*100,digits=2)
pcoa <- as.data.frame(pcoa)
pcoa$samples <- row.names(pcoa)
write.table(pcoa,"PCoA_indoor.xls",sep = "\t",row.names = FALSE, fileEncoding = "UTF-8")

#Anosim
group <- read_excel('path/group.xls',sheet = "time1")
data <- data.raw[,3:204] 
data <- data.raw 
df.dist=vegdist(data)
df.ano1=anosim(df.dist, group$Sample, permutations = 999)
df.ano2=anosim(df.dist, group$School, permutations = 999)
summary(df.ano2)
plot(df.ano) 
df.ano <- data.frame(Sample_x=df.ano1$class.vec, Sample_y=df.ano1$dis.rank, 
                     School_x=df.ano2$class.vec, School_y=df.ano2$dis.rank,
                     Sample_statistic=df.ano1$statistic, Sample_signif=df.ano1$signif,
                     Sample_statistic=df.ano2$statistic, Sample_signif=df.ano2$signif) #输出结果
write.table(df.ano,"Anosim1.xls",sep = "\t",row.names = FALSE, fileEncoding = "UTF-8")
---------
group <- read_excel('path/group.xls',sheet = "indoor")
indoor <- indoor[,3:204] 
indoor <- read_excel('path/Indoor&SHAP.xls',sheet = "indoor_s")
indoor <- column_to_rownames(indoor, "species")
indoor <- t(indoor)
df.dist=vegdist(indoor)
df.ano1=anosim(df.dist, group$Group, permutations = 999)
summary(df.ano1)
plot(df.ano1) 

##Permdisp
group <- read_excel('path/group.xls',sheet = "time1")
data <- data.raw[,-2] 
data <- data.raw[,-1] 
df.betadisper <- betadisper(d = vegdist(data[,-1]), group = data$Sample  , type = "centroid")
df.betadisper <- betadisper(d = vegdist(data[,-1]), group = data$School  , type = "centroid")
-------
df.betadisper <- betadisper(d = vegdist(data[,-1]), group = group$Sample  , type = "centroid")
df.betadisper <- betadisper(d = vegdist(data[,-1]), group = group$School  , type = "centroid")
anova(df.betadisper)
plot(df.betadisper)

##Differential analysis
###FC method
BiocManager::install("limma")
library(limma)
data.raw <- read_excel('path/Indoor&SHAP.xls',sheet = "indoor")
data.raw <- column_to_rownames(data.raw, "ID")
data = log2(data.raw)
data[data == -Inf] = 0
-------
group <- read_excel('path/group.xls', sheet = "indoor")
design <- model.matrix(~0+factor(group$Group))
colnames(design) <- levels(factor(group$Group))
rownames(design) <- colnames(data)
-------
contrast.matrix <- makeContrasts(Slightly_pollution-Moderately_pollution,levels = design)
contrast.matrix <- makeContrasts(Slightly_pollution-Heavy_pollution,levels = design)
-------
fit <- lmFit(data,design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$adj.P.Val > 0.05, "unchanged",
                       ifelse(DEG$logFC > 2, "up-regulated",
                              ifelse(DEG$logFC < -2, "down-regulated", "unchanged")))
write.table(DEG,"FC_SM.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")
write.table(DEG,"FC_SH.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

###diff pathway
output_dir=("path")
install.packages("ReporterScore")
BiocManager::install("S4Vectors")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")
library(ReporterScore)
library(dplyr)
library(pathview)
otu <- read_excel('path/Env.micbio.function.xls', sheet="MAT.kegg.ko")
otu <- column_to_rownames(otu, "KO_ID") 
group <- read_excel('path/group.xls')
group <- column_to_rownames(group, "Sample") 
ko_pvalue=ko.test(otu,"Group",group)
head(ko_pvalue)
ko_stat=pvalue2zs(ko_pvalue,mode="directed")
reporter_s=get_reporter_score(ko_stat)
reporter_res <- combine_rs_res(otu, "Group",group, ko_stat, reporter_s)
write.table(reporter_s,"KeggMap_data.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")
plot_KEGG_map(ko_stat,
              map_id = c("map02030",
                         "map02010",
                         "map02020"), type = "pathway",
              feature = "ko", color_var = "Z_score", 
              save_dir = output_dir,
              color = c("seagreen", "grey", "orange"))

#-----------Heavy metal and Microorganisms------------
library(vegan)
setwd("path")

#input data
data.otu <- read_excel('path/env.micbio.database.p.xls',sheet = "database_p")
data.otu <- column_to_rownames(data.otu, "phylum") 
all.otu <- data.otu[-c(8,32,33),c(3:204)]
all.otu<- decostand(all.otu, method = 'hellinger')
indoor.otu <- subset(data.otu, data.otu$Sample == 1)
indoor.otu <- indoor.otu[,c(3:204)]
indoor.otu<- decostand(indoor.otu, method = 'hellinger')
-------
data.otu <- read_excel('path/env.micbio.database.s.xls',sheet = "database_s")
data.otu <- column_to_rownames(data.otu, "species")
data.otu <- t(data.otu)
all.otu <- data.otu[-c(8,32,33),]
all.otu<- decostand(all.otu, method = 'hellinger')
indoor.otu <- subset(data.otu, grepl("D1I.", rownames(data.otu)))
indoor.otu<- decostand(indoor.otu, method = 'hellinger')
-------
data.env <- read_excel('path/Heavy metal data.xlsx',sheet = "Veen")
all.env <- data.env[,c(1,3:10)] 
all.env <- column_to_rownames(all.env, "ID") 
indoor.env <- data.env[1:41,c(1,3:10)]  
indoor.env <- column_to_rownames(indoor.env, "ID") 

#DCA|RDA
##Mic
otu <- data.otu[,c(3:204)]
ord <- decorana(otu)  
ord <- decorana(data.otu) 
ord
-------
RDA <- rda(all.otu,all.env,scale = T)  
vif.cca(RDA) 
permutest(RDA,permu=999)
r2<-RsquareAdj(RDA)
rda_noadj<-r2$r.squared
rda_adj <- r2$adj.r.squared 
-------
df_rda <- data.frame(RDA$CCA$u[,1:2],rownames(all.env))
colnames(df_rda)=c("RDA1","RDA2","samples")
RDA1 =round(RDA$CCA$eig[1]/sum(RDA$CCA$eig)*100,3)
RDA2 =round(RDA$CCA$eig[2]/sum(RDA$CCA$eig)*100,3)
write.table(df_rda,"indoor_RDA.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")
-------
df_rda_env <- envfit(RDA,all.env,permu=999) 
rda_env <- as.data.frame(RDA$CCA$biplot[,1:2])
rda_env$r2 <- df_rda_env$vectors$r
rda_env$P <- df_rda_env$vectors$pvals
write.table(rda_env,"indoor_env.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")
##Pathway
data.otu <- read_excel('path/Env.micbio.function.xls', sheet="MAT.kegg.ko") 
data.otu <- read_excel('path/Env.micbio.function.xls', sheet="MAT.cazy.level2")
data.otu <- column_to_rownames(data.otu, "KO_ID") 
data.otu <- column_to_rownames(data.otu, "CAZy_Family") 
data.otu <- t(data.otu)
data.otu<- decostand(data.otu, method = 'hellinger')
-------
data.otu <- read_excel('path/Env.micbio.function.xls', sheet="MAT.kegg.ko")
data.otu <- read_excel('path/Env.micbio.function.xls', sheet="MAT.cazy.level2")[,-2]
data.otu <- column_to_rownames(data.otu, "ID") 
data.otu <- column_to_rownames(data.otu, "CAZy_Family") 
data.otu <- t(data.otu)
data.otu<- decostand(data.otu, method = 'hellinger')
-------
ord <- decorana(data.otu)  
RDA <- rda(data.otu,indoor.env,scale = T) 
vif.cca(RDA) 
permutest(RDA,permu=999) 
r2<-RsquareAdj(RDA)
rda_noadj<-r2$r.squared 
rda_adj <- r2$adj.r.squared 
-------
df_rda <- data.frame(RDA$CCA$u[,1:2],rownames(indoor.env))
colnames(df_rda)=c("RDA1","RDA2","samples")
RDA1 =round(RDA$CCA$eig[1]/sum(RDA$CCA$eig)*100,2)
RDA2 =round(RDA$CCA$eig[2]/sum(RDA$CCA$eig)*100,2)
write.table(df_rda,"indoor_RDA_ko.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")
-------
df_rda_env <- envfit(RDA,indoor.env,permu=999) 
rda_env <- as.data.frame(RDA$CCA$biplot[,1:2])
rda_env$r2 <- df_rda_env$vectors$r
rda_env$P <- df_rda_env$vectors$pvals
write.table(rda_env,"indoor_env_ko.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

#Mantel test
devtools::install_github("Hy4m/linKEt",force = TRUE)
library(dplyr)
library(linkET)
data.otu <- read_excel('path/Env.micbio.function.xls',sheet="MAT.kegg.level1")
data.otu <- column_to_rownames(data.otu, "KO_Pathway_Level1") 
all.otu <- data.otu[,c(2:8)]
indoor.otu <- subset(data.otu, data.otu$Sample == "Indoor")
indoor.otu <- indoor.otu[,c(2:8)]
mantel01 <- mantel_test(all.otu, all.env,
                        spec_select = list(Cellular_Processes = 1,Environmental_Information_Processing = 2,Genetic_Information_Processing = 3,Human_Diseases = 4,
                                           Metabolism = 5,Organismal_Systems = 6,Others = 7)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, -0.2, 0, 0.2, Inf),labels = c("< -0.2", "-0.2 to 0","0 to 0.2", ">= 0.2")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),labels = c("< 0.01", "0.01 to 0.05", ">= 0.05")))
mantel02 <- mantel_test(indoor.otu, indoor.env,
                        spec_select = list(Cellular_Processes = 1,Environmental_Information_Processing = 2,Genetic_Information_Processing = 3,Human_Diseases = 4,
                                           Metabolism = 5,Organismal_Systems = 6,Others = 7)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, -0.2, 0, 0.2, Inf),labels = c("< -0.2", "-0.2 to 0","0 to 0.2", ">= 0.2")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),labels = c("< 0.01", "0.01 to 0.05", ">= 0.05")))
write.table(mantel01,"Mantel.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

#-----------Heavy metal and 16s------------
library(vegan)
setwd("path")

#HQ index
data <- read_excel('path/Risk.xlsx',sheet="HQ")
data <- column_to_rownames(data, "ID") 
data <- data[,c(2:11)] 
data_Ing <- data %>%  
  mutate(across(everything(), ~ (.x*200*230*10^-6) / (weight*365), .names = "{col}_Ing")) 
data_InH <- data %>%  
  mutate(across(everything(), ~ (.x*7.6*230) / (weight*365*1.36*10^9), .names = "{col}_Inh"))  
data_dermal <- data %>%  
  mutate(across(everything(), ~ (.x*SA_cm2*0.2*0.001*190*10^-6) / (weight*365), .names = "{col}_dermal"))
data_Hr <- cbind(data_Ing[13:20],data_InH[13:20],data_dermal[13:20])
write.table(data_Hr,"HR2.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

#Alpha analysis
data <- read_excel('path/16s_otu_table.xls',sheet = "p")[,1:9]
data <- read_excel('path/16s_otu_table.xls',sheet = "s")[,-126]
data <- column_to_rownames(data, "ID")
shannon=diversity(data,"shannon")
simpson=diversity(data,"simpson")
write.table(shannon,"Shannon.xls",sep = "\t",row.names = TRUE, fileEncoding = "UTF-8")
write.table(simpson,"Simpson.xls",sep = "\t",row.names = TRUE, fileEncoding = "UTF-8")

#PCoA
data <- read_excel('path/16s_otu_table.xls',sheet = "p")[,1:9]
data <- read_excel('path/16s_otu_table.xls',sheet = "s")[,-126]
data <- column_to_rownames(data, "ID")
otu.distance <- vegdist(data,method='bray')
df.pcoa <- cmdscale(otu.distance,eig=TRUE,k=4)
pcoa <- df.pcoa$points[,1:4]
pc <- round(df.pcoa$eig/sum(df.pcoa$eig)*100,digits=2)
pcoa <- as.data.frame(pcoa)
pcoa$samples <- row.names(pcoa)
write.table(pcoa,"PCoA_p.xls",sep = "\t",row.names = FALSE, fileEncoding = "UTF-8")

##Anosim
data <- read_excel('path/16s_otu_table.xls',sheet = "p")
data <- read_excel('path/16s_otu_table.xls',sheet = "s")
data <- column_to_rownames(data, "ID")
df.dist=vegdist(data[,-9])
df.dist=vegdist(data[,-125])
df.ano1=anosim(df.dist, data$group, permutations = 999)#p
df.ano2=anosim(df.dist, data$group, permutations = 999)#s
summary(df.ano2)

###FC methods
library(limma)
data.raw <- read_excel('path/16s_otu_table.xls',sheet = "otu")
data.raw <- column_to_rownames(data.raw, "#OTU ID")
data = log2(data.raw)
data[data == -Inf] = 0
-------
group <- read_excel('path/group.xls',sheet = "16s")
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(data)
-------
contrast.matrix <- makeContrasts(Slightly_pollution-Moderately_pollution,levels = design)
contrast.matrix <- makeContrasts(Slightly_pollution-Heavy_pollution,levels = design)
-------
fit <- lmFit(data,design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 2, "up-regulated",
                              ifelse(DEG$logFC < -2, "down-regulated", "unchanged")))
write.table(DEG,"FC_SM.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")
write.table(DEG,"FC_SH.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

###FC methods-pathway
library(limma)
####GO
data.raw <- read_excel('path/Meta.function.xls',sheet="meta_GO")
data.raw <- column_to_rownames(data.raw, "#OTU ID")
data = log2(data.raw)
data[data == -Inf] = 0
-------
group <- read_excel('path/group.xls', sheet="meta")
design <- model.matrix(~0+factor(group$Group))
colnames(design) <- levels(factor(group$Group))
rownames(design) <- colnames(data)
contrast.matrix <- makeContrasts(Moderate_pollution-Heavy_pollution,levels = design)
-------
fit <- lmFit(data,design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 2, "up-regulated",
                              ifelse(DEG$logFC < -2, "down-regulated", "unchanged")))
write.table(DEG,"GO_res.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")
####kegg
data.raw <- read_excel('path/Meta.function.xls',sheet = "meta_ko")
data.raw <- column_to_rownames(data.raw, "#ko")
data = log2(data.raw)
data[data == -Inf] = 0
-------
group <- read_excel('path/group.xls',sheet="meta")
design <- model.matrix(~0+factor(group$Group))
colnames(design) <- levels(factor(group$Group))
rownames(design) <- colnames(data)
contrast.matrix <- makeContrasts(Moderate_pollution-Heavy_pollution,levels = design)
-------
fit <- lmFit(data,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 2, "up-regulated",
                              ifelse(DEG$logFC < -2, "down-regulated", "unchanged")))
write.table(DEG,"KEGG_res.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

#DCA|RDA
##16s
otu <- read_excel('path/16s_otu_table.xls',sheet = "p")[,-10]
otu <- read_excel('path/16s_otu_table.xls',sheet = "s")[,-126]
otu <- column_to_rownames(otu, "ID") 
otu <- decostand(otu, method = 'hellinger')
env <- read_excel('path/HR.xls',sheet = "HI")
env <- column_to_rownames(env, "ID") 
env <- log10(env)
-------
ord <- decorana(otu)  
RDA <- rda(otu,env,scale = T) 
vif.cca(RDA) 
RDA <- rda(otu,env[,-5:-6],scale = T) 
permutest(RDA,permu=999) 
r2<-RsquareAdj(RDA)
rda_noadj<-r2$r.squared 
rda_adj <- r2$adj.r.squared 
-------
df_rda <- data.frame(RDA$CCA$u[,1:2],rownames(env))
colnames(df_rda)=c("RDA1","RDA2","samples")
RDA1 =round(RDA$CCA$eig[1]/sum(RDA$CCA$eig)*100,2)
RDA2 =round(RDA$CCA$eig[2]/sum(RDA$CCA$eig)*100,2)
write.table(df_rda,"RDAp.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")
-------
df_rda_env <- envfit(RDA,env[,-5:-6],permu=999)
rda_env <- as.data.frame(RDA$CCA$biplot[,1:2])
rda_env$r2 <- df_rda_env$vectors$r
rda_env$P <- df_rda_env$vectors$pvals
write.table(rda_env,"indoor_env.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")
##pathway
sample <- c("choosed ID")
data.go <- read_excel('path/function.xls',sheet = "meta_GO")
data.go <- column_to_rownames(data.go, "#OTU ID")
colnames(data.go) <- sample
data.go <- t(data.go)
data.go<- decostand(data.go, method = 'hellinger')
-------
data.kegg <- read_excel('path/Meta.function.xls',sheet="meta_ko")
data.kegg <- column_to_rownames(data.kegg, "#ko")
colnames(data.kegg) <- sample
data.kegg <- t(data.kegg)
data.kegg<- decostand(data.kegg, method = 'hellinger')
-------
env <- read_excel('path/HR.xls',sheet = "HI")
meta.env <- subset(env, ID %in% sample)
meta.env <- column_to_rownames(meta.env, "ID") 
meta.env <- log10(meta.env)
-------
ord <- decorana(data.kegg)  
RDA <- rda(data.go,meta.env,scale = T) 
vif.cca(RDA) 
RDA <- rda(data.go,meta.env[,c(-5,-6,-8)],scale = T) 
permutest(RDA,permu=999) 
r2<-RsquareAdj(RDA)
rda_noadj<-r2$r.squared 
rda_adj <- r2$adj.r.squared 
-------
df_rda <- data.frame(RDA$CCA$u[,1:2],rownames(meta.env))
colnames(df_rda)=c("RDA1","RDA2","samples")
RDA1 =round(RDA$CCA$eig[1]/sum(RDA$CCA$eig)*100,2)
RDA2 =round(RDA$CCA$eig[2]/sum(RDA$CCA$eig)*100,2)
write.table(df_rda,"RDAp.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")
-------
df_rda_env <- envfit(RDA,meta.env[,c(-5,-6,-8)],permu=999)
rda_env <- as.data.frame(RDA$CCA$biplot[,1:2])
rda_env$r2 <- df_rda_env$vectors$r
rda_env$P <- df_rda_env$vectors$pvals
write.table(rda_env,"indoor_env.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")

###diff pathway
output_dir=setwd("path")
library(ReporterScore)
library(dplyr)
library(pathview)
otu <- read_excel('path/Meta.function.xls',sheet="meta_ko")
otu <- column_to_rownames(otu, "#ko") 
group <- read_excel('path/group.xls',sheet="meta")
group <- column_to_rownames(group, "Sample") 
####diff kegg map
ko_pvalue=ko.test(otu,"Group",group)
head(ko_pvalue)
ko_stat=pvalue2zs(ko_pvalue,mode="directed")
reporter_s=get_reporter_score(ko_stat)
reporter_res <- combine_rs_res(otu, "Group",group, ko_stat, reporter_s)
write.table(reporter_s,"KeggMap_data.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")
plot_KEGG_map(ko_stat,
              map_id = c(
                         "map02060",
                         "map01230" ),type = "pathway",
              feature = "ko", color_var = "Z_score", 
              save_dir = output_dir,
              color = c("seagreen", "grey", "orange"))

#Correlation analysis 
library(dplyr)
library(linkET)
##function
sample <- c("choosed ID")
otu <- read_excel('path/Meta.function.xls',sheet="meta_pathway")
otu <- column_to_rownames(otu, "LEVEL1")
colnames(otu) <- sample
otu <- t(otu)
--------
env <- read_excel('path/Risk.xlsx',sheet = "HQ")[,c(1,5:12)]
env <- subset(env, ID %in% sample)
env <- column_to_rownames(env, "ID")
mantel01 <- mantel_test(otu, env,
                        spec_select = list(Cellular_Processes = 1,Environmental_Information_Processing = 2,Genetic_Information_Processing = 3,
                                           Human_Diseases = 4, Metabolism = 5,Organismal_Systems = 6)) %>%
            mutate(rd = cut(r, breaks = c(-Inf, -0.2, 0, 0.2, Inf),labels = c("< -0.2", "-0.2 to 0","0 to 0.2", ">= 0.2")),
                   pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),labels = c("< 0.01", "0.01 to 0.05", ">= 0.05")))
write.table(mantel01,"mantel.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")

##diff 16s
otu <- read_excel('path/16s_otu_table.xls',sheet = "DEdata")
otu <- column_to_rownames(otu, "#OTU ID") 
otu <- t(otu)
env <- read_excel('path/Risk.xlsx',sheet = "HQ")[,c(1,9:16)]
env <- column_to_rownames(env, "ID")
--------
mantel02 <- mantel_test(otu, env,
                          spec_select = list(g_Collinsella = 1,f_Bacteroidaceae = 2,g_Bacteroides = 3,
                                             f_Leuconostocaceae = 4, g_Weissella = 5,g_Intestinibacter = 6,
                                             g_Butyricicoccus = 7, g_Erysipelatoclostridium = 8,g_Holdemanella = 9)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, -0.2, 0, 0.2, Inf),labels = c("< -0.2", "-0.2 to 0","0 to 0.2", ">= 0.2")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),labels = c("< 0.01", "0.01 to 0.05", ">= 0.05")))  
write.table(mantel02,"Cor.xls",sep = "\t",row.names = F, fileEncoding = "UTF-8")

#-----------Multi omics analysis-----------
##WGCNA
setwd("path/WGCNA")
library(WGCNA)
###input data
allowWGCNAThreads()
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
gene <- read_excel('path/WGCNA.xls')
gene <- column_to_rownames(gene, "ID") 
sam <- read_excel('path/16s_otu_table.xlsx',sheet = "otu")[,-11]
sam <- column_to_rownames(sam, "ID") 
m.vars <- apply(gene,1,var)
gene.upper <- gene[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
dim(gene.upper)
datgene <- as.data.frame(t(gene.upper))
nGenes <- ncol(datgene)
nSamples <- nrow(datgene)
gsg <- goodSamplesGenes(datgene, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datgene), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
###powers
powers = c(seq(1,10,by = 1), seq(12, 20, by = 2))
sft = pickSoftThreshold(datgene, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0,col="red")
###self-correlation
sam_Tree = hclust(dist(sam), method = "average")
traitColors = numbers2colors(sam, signed = FALSE)
tiff(filename="Sam dendrogram.tiff",width = 14, height = 16, units = "cm", res = 600)
plotDendroAndColors(sam_Tree, traitColors, 
                    groupLabels = names(sam), 
                    dendroLabels = FALSE, main = "Gut microbiota dendrogram and trait heatmap")
dev.off() 
###network
net = blockwiseModules(datgene, power = 10, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)
tiff(filename="dendrogram.tiff",width = 14, height = 10, units = "cm", res = 600)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off() 
###Module-trait
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs0 = moduleEigengenes(datgene, moduleColors)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, sam, use = "p")
colnames(MEsWW)
modlues=MEsWW
----------
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
tiff(filename="Module-trait.tiff",width = 16, height = 14, units = "cm", res = 600)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(sam), yLabels = names(MEsWW), cex.lab = 0.9,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off() 
###Output Module
library(stringr)
tiff(filename="Module-Module.tiff",width = 14, height = 16, units = "cm", res = 600)
plotEigengeneNetworks(MEsWW, "Eigengene adjacency heatmap", 
                      marDendro = c(0.1,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,  # marDendro/marHeatmap 设置下、左、上、右的边距
                      xLabelsAngle = 90)
dev.off() 
TOM = TOMsimilarityFromExpr(datgene, power = 10 )
save.image(file="project_image.RData")
load(file="project_image.RData")
module = "red"
module = "black"
module = "purple"
module = "yellow"
module = "brown"
module = "green"
module = "blue"
module = "pink"
module = "magenta"
module = "turquoise"
probes  = colnames(datgene) 
inModule = (moduleColors==module)
modProbes = probes[inModule]
--------
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt <- exportNetworkToCytoscape(modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule])
write.table(cyt$edgeData,"edge_pink.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8") #其实不用这个步骤已经自动输出txt

##RLQ and Fourth-corner
setwd("path")
library(vegan)
library(ade4)
library(WGCNA)
###input data
spe <- read_excel('path/WGCNA.xls',sheet="WGCNA")
spe <- column_to_rownames(spe, "ID")
spe <- as.data.frame(t(spe))
env <- read_excel('path/HR.xls',sheet = "HI")
env <- column_to_rownames(env, "ID") 
data <- read_excel('path/16s_otu_table.xls',sheet = "p")[,-11] 
data <- column_to_rownames(data, "ID") 
traits <- as.data.frame (cor(spe,data,method="spearman"))
###CA
afcL.aravo <- dudi.coa(spe, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(env, row.w = afcL.aravo$lw, scannf = FALSE)
acpQ.aravo <- dudi.pca(traits, row.w = afcL.aravo$cw,scannf = FALSE)
###RLQ
rlq.aravo <- rlq(
  dudiR = acpR.aravo, 
  dudiL = afcL.aravo, 
  dudiQ = acpQ.aravo,
  scannf = FALSE)
site<-rlq.aravo$lR
otu<-rlq.aravo$lQ
chem<-rlq.aravo$l1*6
trait<-rlq.aravo$c1*6
write.table(trait,"trait.xls",sep = "\t",row.names = T, fileEncoding = "UTF-8")
AxcR1 <- round((rlq.aravo$eig[1] / sum(rlq.aravo$eig)),2)
AxcR2 <- round((rlq.aravo$eig[2] / sum(rlq.aravo$eig)),2)

##fourth-corner
fourth.aravo <- fourthcorner(
  tabR = env, 
  tabL = spe, 
  tabQ = traits,
  modeltype = 6,
  p.adjust.method.G = "none", 
  p.adjust.method.D = "none", 
  nrepet = 999)
###FDR
fourth.aravo.adj <- p.adjust.4thcorner(fourth.aravo,
                    p.adjust.method.G = "fdr", p.adjust.method.D = "fdr", p.adjust.D = "global") 
fourth.aravo.adj
plot(fourth.aravo.adj,
     alpha = 0.05,
     stat = "D2")
tiff(filename="fourth_net.tiff",width = 13, height = 13, units = "cm", res = 600)
plot(fourth.aravo.adj,x.rlq = rlq.aravo,alpha = 0.1,stat = "D2",type = "biplot")
dev.off() 

##procrustes
library(vegan)
env <- read_excel('path/HR.xls',sheet = "HI")
env <- column_to_rownames(env, "ID")
spe <- read_excel('path/16s_otu_table.xls',sheet = "p")[,-11] 
spe <- column_to_rownames(data, "ID") 
spe.dist <- vegdist(spe)  #物种一般用Bray-Curtis
env.dist <- vegdist(env)  #环境一般使用欧式距离
mds.s <- monoMDS(spe.dist)
mds.e <- monoMDS(env.dist)
pro.s.e <- procrustes(mds.s,mds.e, symmetric = TRUE)
set.seed(1) 
pro.s.e_t <- protest(mds.s,mds.e, permutations = 999)
Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
Pro_X <- data.frame(pro.s.e$rotation)

#-----------Machine learning-----------
##XGboost
setwd("path/XGboost")
install.packages("shapviz")
install.packages("xgboost")
install.packages("caret")
library(shapviz)
library(xgboost)
library(caret)
###input
all.otu1 <- read_excel('path/Indoor&SHAP.xlsx',sheet = "SHAP_p")[c(-8,-32,-33),c(1,4)]
all.otu2 <- read_excel('path/Indoor&SHAP.xlsx',sheet = "SHAP_s")[c(-8,-32,-33),c(1,4)]
-------
all.otu1 <- read_excel('path/Indoor&SHAP.xlsx',sheet = "SHAP_indoorp")[,c(1,4)]
all.otu2 <- read_excel('path/Indoor&SHAP.xls',sheet = "SHAP_indoors")[,c(1,4)]
-------
all.otu1 <- read_excel('path/Indoor&SHAP.xlsx',sheet = "SHAP_p")[,c(1,4)]
all.otu2 <- read_excel('path/Indoor&SHAP.xlsx',sheet = "SHAP_s")[,c(1,4)]
all.env <- read_excel('path/HR.xls',sheet = "HI")
-------
data.env <- read_excel('path/Heavy metal data.xlsx',sheet = "Veen")
all.env <- data.env[,c(1,3:10)] 
indoor.env <- data.env[1:41,c(1,3:10)]  
###merge
all.data <- merge(all.otu1,all.env,by="ID")
all.data <- merge(all.otu2,all.env,by="ID")
-------
all.data <- merge(all.otu1,indoor.env,by="ID")
all.data <- merge(all.otu2,indoor.env,by="ID")
###Divide the training set and testing set
all.data$Shannon <- all.data$Shannon*100
all.data <- all.data[,-1]
set.seed(123)
inTrain<-createDataPartition(y=all.data[,"Shannon"], p=0.8, list=F) 
traindata<-all.data[inTrain,]
testdata<-all.data[-inTrain,]
###Parameter Range
grid<-expand.grid(nrounds=c(75,100,150),  
                  colsample_bytree=1,  
                  min_child_weight=1, 
                  eta=c(0.01,0.1,0.3),  
                  gamma=c(0.5,0.25,0.1), 
                  subsample=0.5,
                  max_depth=c(2,3,4))
cntrl<-trainControl(method="cv",
                    number=5,
                    verboseIter=F,
                    returnData=F,
                    returnResamp="final")
set.seed(123)
train.xgb<-train(x=traindata[,2:9],
                 y=traindata[,1],  
                 trControl=cntrl,  
                 tuneGrid=grid,
                 method="xgbTree")
train.xgb
###building model
model_xgboost = xgboost(  #p
  data = as.matrix(traindata[,c(2:(ncol(traindata)))]),
  label = traindata$Shannon,
  nrounds = 75,
  max_depth = 4,
  eta  = 0.1,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective= "reg:squarederror" )
model_xgboost = xgboost(  #s
  data = as.matrix(traindata[,c(2:(ncol(traindata)))]),
  label = traindata$Shannon,
  nrounds = 100,
  max_depth = 2,
  eta  = 0.3,
  gamma = 0.25,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective= "reg:squarederror" )
model_xgboost = xgboost(  #indoor-p
  data = as.matrix(traindata[,c(2:(ncol(traindata)))]),
  label = traindata$Shannon,
  nrounds = 75,
  max_depth = 4,
  eta  = 0.1,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective= "reg:squarederror" )
model_xgboost = xgboost(  #indoor-s
  data = as.matrix(traindata[,c(2:(ncol(traindata)))]),
  label = traindata$Shannon,
  nrounds = 75,
  max_depth = 4,
  eta  = 0.3,
  gamma = 0.5,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective= "reg:squarederror" )
model_xgboost = xgboost(  #16s-p
  data = as.matrix(traindata[,c(2:(ncol(traindata)))]),
  label = traindata$Shannon,
  nrounds = 75,
  max_depth = 2,
  eta  = 0.1,
  gamma = 0.5,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective= "reg:squarederror" )
model_xgboost = xgboost(  #16s-s
  data = as.matrix(traindata[,c(2:(ncol(traindata)))]),
  label = traindata$Shannon,
  nrounds = 75,
  max_depth = 4,
  eta  = 0.1,
  gamma = 0.1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective= "reg:squarederror" )
###SHAP
shap_xgboost=shapviz(model_xgboost,X_pred=as.matrix(traindata[,c(2:(ncol(traindata)))]))
