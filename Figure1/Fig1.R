rm(list=ls())

library(tidyverse)
library(reshape2)
library(readxl)
library(ggplot2)

help(" ")

x <- colorRampPalette(c("#F9035A","#049AF5"))(16)
scales::show_col(x)

#F1.A
df <- read_excel('path/Heavy metal data.xls',sheet = "nosign")
df <- df[,c(1,3:11)]
df <- melt(df, id.vars = c("CID","Sample"))
p1 <- ggplot(df, aes(x=variable, y=value,fill=Sample))+
  stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.9))+
  geom_boxplot(position=position_dodge(0.9),size=0.3,width=0.8,color="black",
               outlier.color = "black",outlier.stroke = 0.3,
               notch = F, notchwidth = 0.5)+
  scale_fill_manual(values = c("#ff6100","#72C3A3","#5861AC"))+
  theme_bw()+
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10))+
  labs(x="", y = "Metal content (μg/g)")+
  scale_y_log10()
ggsave("Sample_bar.tiff", p1, width =20, height =10, units = "cm", dpi = 600)

#F1.B
data <- read_excel('path/Heavy metal data.xls',sheet = "nosign")
data <- data[,c(2,4:11)]
data <- column_to_rownames(data, "CID2") 
Cor.r <- as.matrix(cor(data, method = "spearman"))
Cor.p <- cor.mtest(data, conf.level = .95)
tiff(filename="Cor.tiff",width = 10, height = 10, units = "cm", res = 600)
corrplot(Cor.r, method = c('pie'),type = c('upper'), 
         col =c("#009ACD","#06a7cd","#37b8d6","#69c9e0","#9adbea",
                "#f5b5ac","#f09283","#eb6d5a","#e74a32","#CD3333"),
         outline = 'grey', 
         diag = T,
         tl.cex = 0.5, tl.col = 'black',tl.pos = 'tp',
         p.mat = Cor.p$p,
         sig.level = c(.001, .01, .05),
         insig = "label_sig", pch.cex = 0.5, pch.col = 'black')
corrplot(Cor.r, add = TRUE, method = c('number'), 
         type = c('lower'),
         col =c("#009ACD","#06a7cd","#37b8d6","#69c9e0","#9adbea",
                "#f5b5ac","#f09283","#eb6d5a","#e74a32","#CD3333"),
         diag = FALSE, 
         number.cex = 0.5,
         tl.pos = 'n', cl.pos = 'n',
         p.mat = Cor.p$p,
         insig = "pch",pch.cex = 0.5,pch.col = 'black')
dev.off() 

#F1.C
df <- read_excel("path/PCA.xls")
p1<-ggplot(data=df)+
  theme_bw()+
  geom_point(aes(x=V1,y=V2, color=Sample),size=3, shape=16)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x="1 43.93% ",y="2 18.95%",color="Group",fill="Group")+
  scale_color_manual(values = c("#ff6100","#72C3A3","#5861AC")) +
  scale_fill_manual(values = c("#ff6100","#72C3A3","#5861AC"))+
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10),
        legend.position = "bottom",
        panel.grid=element_blank())+
  stat_ellipse(data=df, geom="polygon", level=0.95,linetype=2, size=1,
               aes(x=V1,y=V2, fill=Sample),alpha=0.1,show.legend = T)
ggsave("Sample_PCA.tiff", p1, width =13, height = 15, units = "cm", dpi = 600)

#F1.D
df <- read_excel("path/Risk.xls",sheet = "PLI")
df$Sample <- factor(df$Sample,levels = c('Air-conditioning','Outdoor','Indoor')) 
p1<- ggplot(df, aes(x=Sample, y=PLI,fill=Sample))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 2, alpha = 0.15,fill="#d7ebce") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = 3, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3, ymax = Inf, alpha = 0.2,fill="#ffdc80")+
  stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.9))+
  geom_boxplot(position=position_dodge(0.9),size=0.3,width=0.8,color="black",
               outlier.color = "black",outlier.stroke = 0.3,
               notch = F, notchwidth = 0.5)+
  scale_fill_manual(values = c("#ff6100","#5861AC","#72C3A3"))+
  theme_bw()+
  theme(axis.title = element_text(size = 15), axis.text = element_text(size=11),
        legend.position = "bottom")+
  labs(x="", y = "The pollution load index")
ggsave("PLI_Metal.tiff", p1, width =7.7, height = 15, units = "cm", dpi = 600)
---------
df <- read_excel("path/PCA.xls",sheet = "PCA_indoor")
p1<-ggplot(data=df)+
  theme_bw()+
  geom_point(aes(x=V1,y=V2, color=group),size=3, shape=16)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x="1 52.93% ",y="2 14.99%",color="Group",fill="Group")+
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080")) +
  scale_fill_manual(values =c('#A5C2E2',"#f28080","#FFc080"))+
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10),
        legend.position = "bottom",
        panel.grid=element_blank())+
  stat_ellipse(data=df, geom="polygon", level=0.95,linetype=2, size=1,
               aes(x=V1,y=V2, fill=group),alpha=0.1,show.legend = T)
ggsave("PLI_PCA.tiff", p1, width =13, height = 15, units = "cm", dpi = 600)

#F1.E
df <- read_excel('path/filldata.xls',sheet = "data_p")
df <- melt(df, id.vars = c("phylum","Sample"))
---------
df <- read_excel('path/filldata.xls',sheet = "data_s")
df <- melt(df, id.vars = c("species","Sample"))
group <- read_excel('path/group.xls',sheet="indoor")
df <- subset(df, df$Sample == "Indoor")
names(df)=c("ID","Sample","variable","value")
df <- merge(df,group,by="ID")
p1 <- ggplot(df) +
  geom_col(aes(Group, value, fill = variable),position = 'fill',width = 0.9, size=0.2) +
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values =  rev(c("#FF6100",'#FA8072',"#FF7F00",'#FEA040','#F5DEB3',"#d6e2e2","#B2DFEE",'#A5C2E2',"#72C3A3","#6AB4CA",'#5861AC'))) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.02))) +
  labs(x = '', y = 'Relative Abundance') +
  theme_bw()+
  theme(axis.title.y = element_text(size = 15), axis.text = element_text(size=12),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = 11),legend.title = element_blank(), 
        legend.key.size = unit(0.25, "cm"), legend.key.height = unit(0.25, "cm"),legend.key.width = unit(0.25,"cm"))
p2 <- ggplot(df) +
  geom_col(aes(Group, value, fill = variable),position = 'fill',width = 0.9, size=0.2) +
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values =  rev(c("#FF6100",'#FA8072',"#FF7F00",'#FEA040','#F5DEB3',"#d6e2e2","#B2DFEE",'#A5C2E2',"#6AB4CA",'#5861AC',"#72C3A3"))) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.02))) +
  labs(x = '', y = 'Relative Abundance') +
  theme_bw()+
  theme(axis.title.y = element_text(size = 15), axis.text = element_text(size=12),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent'), 
        legend.text = element_text(size = 11),legend.title = element_blank(), 
        legend.key.size = unit(0.25, "cm"), legend.key.height = unit(0.25, "cm"),legend.key.width = unit(0.25,"cm"))
ggsave("Indoor_pfill.tiff", p1, width =17.9, height =13, units = "cm", dpi = 600)
ggsave("Indoor_sfill.tiff", p2, width =20.1, height =13, units = "cm", dpi = 600)

#F1.F
df <- read_excel("path/PCoA.xls", sheet = "indoor_p")
df <- read_excel("path/PCoA.xls", sheet = "indoor_s")
p1<-ggplot(data=df)+
  theme_bw()+
  geom_point(aes(x=V1,y=V2, color=Group),size=3, shape=16)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x="PCo1 86.18% ",y="PCo2 7.96%",color="Group",fill="Group")+ 
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080")) +
  scale_fill_manual(values = c('#A5C2E2',"#f28080","#FFc080"))+
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10),
        plot.title = element_text(size=13,hjust = 0.5),
        legend.position = "bottom",
        panel.grid=element_blank())+
  stat_ellipse(data=df, geom="polygon", level=0.9,linetype=2, size=1,
               aes(x=V1,y=V2, fill=Group),alpha=0.1,show.legend = T)
ggsave("Indoor_pPCoA.tiff", p1, width =13, height = 14.5, units = "cm", dpi = 600)

#F1.G
library(ggbreak)
data <- read_excel("path/FC.res.SvsM.xls")
data <- read_excel("path/FC.res.SvsH.xls")
p1<- ggplot(data, aes(x=site, y=Log10P))+
  geom_point(aes(fill =group, color=group, shape=regulate),alpha=0.9,size=3)+
  scale_shape_manual(values = c("down-regulated" = 25, "up-regulated" = 24, "unchanged" = 21)) +
  scale_fill_manual(values = c('#FEA040',"#d6e2e2","#B2DFEE",'#F5DEB3',"#FF7F00",'#A5C2E2'))+
  scale_color_manual(values = c('#FEA040',"#d6e2e2","#B2DFEE",'#F5DEB3',"#FF7F00",'#A5C2E2'))+
  labs(x = "OTU", y = "-log10(P-value)", fill = "level",color = "level")+
  geom_hline(yintercept= 1.3, linetype= "dashed",linewidth=0.5)+
  theme_bw()+
  theme(legend.title =element_text(size =14),legend.text = element_text(size =11),
        axis.text =element_text(size =12), axis.title =element_text(size =15),
        panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(0,4))+
  scale_y_break(c(0.4,0.4),space = 0.01, scales = 1, ticklabels = NULL)+ 
  scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), breaks = c(100,300,523,734,947,1180),labels = c("Phylum","Class","Order","Family","Genus","Spesis"))
ggsave("SM_manha.tiff", p1, width =20, height = 10, units = "cm", dpi = 600)
---------
  p2<- ggplot(data, aes(x=site, y=Log10P))+
  geom_point(aes(fill =group, color=group, shape=regulate),alpha=0.9,size=3)+
  scale_shape_manual(values = c("down-regulated" = 25, "up-regulated" = 24, "unchanged" = 21)) +
  scale_fill_manual(values = c('#FEA040',"#d6e2e2","#B2DFEE",'#F5DEB3',"#FF7F00",'#A5C2E2'))+
  scale_color_manual(values = c('#FEA040',"#d6e2e2","#B2DFEE",'#F5DEB3',"#FF7F00",'#A5C2E2'))+
  labs(x = "OTU", y = "-log10(P-value)", fill = "level",color = "level")+
  geom_hline(yintercept= 1.3, linetype= "dashed",linewidth=0.5)+
  theme_bw()+
  theme(legend.title =element_text(size =14),legend.text = element_text(size =11),
        axis.text =element_text(size =12), axis.title =element_text(size =15),
        panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(0,4))+
  scale_y_break(c(0.59,0.59),space = 0.01, scales = 1, ticklabels = NULL)+ 
  scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), breaks = c(100,300,523,734,947,1180),labels = c("Phylum","Class","Order","Family","Genus","Spesis"))
ggsave("SH_manha.tiff", p2, width =20, height = 10, units = "cm", dpi = 600)

#F1.H
data <- read_excel('path/LEfSe.res.xls',sheet = "DEotu")
p1 <- ggplot(data,aes(x=reorder(Biomarker, order)))+
  scale_fill_manual(values =  c('#A5C2E2',"#f28080","#FFc080"))+
  geom_bar(aes(y=LDAvalue,fill=Groups),stat = "identity",color="black", size=0.2)+
  coord_flip()+ 
  ylab("LDA SCORE (log 10)")+
  xlab("")+
  scale_y_continuous(expand = expansion(mult = c(0,0)), limits = c(0,4.5),breaks = seq(0,4.5,1))+
  geom_hline(yintercept= 2, linetype= "dashed",linewidth=0.5)+
  geom_text(aes(y = LDAvalue +0.2, label = sign), size =4, hjust =0.5, vjust =0.7) +
  theme_classic()+
  theme(axis.line = element_line(size=0.25),
        legend.title= element_blank(), legend.position = "top",
        legend.key.size = unit(0.25,"cm"), legend.key.height = unit(0.25,"cm"), legend.key.spacing.y = unit(0.25,"cm"))
ggsave("indoor_Lefse.tiff", p1, width =10, height = 7, units = "cm", dpi = 600)
--------
data <- read_excel('path/LEfSe.res.xls',sheet = "DEkegg")
p1 <- ggplot(data)+
  ggforce::geom_link(aes(x=0,y=reorder(KO_NAME, level),
                         xend =LDA_value,yend = KO_NAME,
                         color = Groups, size = after_stat(index)),
                     n= 1000, show.legend =F)+
  geom_point(aes(x= LDA_value,y= KO_NAME),color="black",fill="white",
             size=5.5,shape =21)+
  geom_line(aes(x=Gene_Num2,y=KO_NAME,group=1),
            orientation = "y",linewidth=1,color="#FFC24B")+
  geom_point(aes(x=Gene_Num2,y=KO_NAME,group=1),color="#FFC24B", size=2.5)+
  geom_vline(xintercept= 2, linetype= "dashed",linewidth=0.5)+
  geom_text(aes(x = LDA_value +0.4,y=reorder(KO_NAME,level), label = sign), size =3.5, hjust =0.5, vjust =0.7) +
  scale_x_continuous(sec.axis = sec_axis(~.))+
  labs(x="LDA SCORE(log 10)",y=NULL,title = "Gene number per 1000",vjust =0.5)+
  theme_bw()+
  theme(axis.text =element_text(size =12,color = "black"),
        axis.title.x=element_text(size = 14,color = "black"),
        axis.line = element_line(size=0.25),
        panel.grid = element_blank(), panel.background = element_blank())+
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080"))
ggsave("Kegg_LDA_indoor.tiff", p1, width =10, height = 7, units = "cm", dpi = 600)
--------
data <- read_excel('path/LEfSe.res.xls',sheet = "DEcazy")
p1 <- ggplot(data)+
  ggforce::geom_link(aes(x=0,y=reorder(Biomarker, level),
                         xend =LDA_value,yend = Biomarker,
                         color = Groups, size = after_stat(index)),
                     n= 1000, show.legend =F)+
  geom_point(aes(x= LDA_value,y= Biomarker),color="black",fill="white",
             size=5.5,shape =21)+
  geom_line(aes(x=Gene_Num2,y=Biomarker,group=1),
            orientation = "y",linewidth=1,color="#FFC24B")+
  geom_point(aes(x=Gene_Num2,y=Biomarker,group=1),color="#FFC24B", size=2.5)+
  geom_vline(xintercept= 2, linetype= "dashed",linewidth=0.5)+
  geom_text(aes(x = LDA_value +0.4,y=reorder(Biomarker,level), label = sign), size =3.5, hjust =0.5, vjust =0.7) +
  scale_x_continuous(sec.axis = sec_axis(~.))+
  labs(x="LDA SCORE(log 10)",y=NULL,title = "Gene number per 1000000",vjust =0.5)+
  theme_bw()+
  theme(axis.text =element_text(size =12,color = "black"),
        axis.title.x=element_text(size = 14,color = "black"),
        axis.line = element_line(size=0.25),
        panel.grid = element_blank(), panel.background = element_blank())+
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080"))
ggsave("CAZy_LDA_indoor.tiff", p1, width =8, height = 7.5, units = "cm", dpi = 600)
--------
data <- read_excel('path/LEfSe.res.xls',sheet = "DEnog")
p1 <- ggplot(data)+
  ggforce::geom_link(aes(x=0,y=reorder(Description, level),
                         xend =LDA_value,yend = Description,
                         color = Groups, size = after_stat(index)),
                     n= 1000, show.legend =F)+
  geom_point(aes(x= LDA_value,y= Description),color="black",fill="white",
             size=5.5,shape =21)+
  geom_line(aes(x=Gene_Num2,y=Description,group=1),
            orientation = "y",linewidth=1,color="#FFC24B")+
  geom_point(aes(x=Gene_Num2,y=Description,group=1),color="#FFC24B", size=2.5)+
  geom_vline(xintercept= 2, linetype= "dashed",linewidth=0.5)+
  geom_text(aes(x = LDA_value +0.4,y=reorder(Description,level), label = sign), size =3.5, hjust =0.5, vjust =0.7) +
  scale_x_continuous(sec.axis = sec_axis(~.))+
  labs(x="LDA SCORE(log 10)",y=NULL,title = "Gene number per 1000000",vjust =0.5)+
  theme_bw()+
  theme(axis.text =element_text(size =12,color = "black"),
        axis.title.x=element_text(size = 14,color = "black"),
        axis.line = element_line(size=0.25),
        panel.grid = element_blank(), panel.background = element_blank())+
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080"))
ggsave("NOG_LDA_indoor.tiff", p1, width =16, height = 8.5, units = "cm", dpi = 600)
