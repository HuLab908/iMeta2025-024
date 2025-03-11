rm(list=ls())

library(tidyverse)
library(reshape2)
library(readxl)
library(ggplot2)

help(" ")

x <- colorRampPalette(c("#F9035A","#049AF5"))(16)
scales::show_col(x)

#F2.A
df <- read_excel("path/RDA.xls",sheet = "RDA_pall")
df <- read_excel("path/RDA.xls",sheet = "RDA_sall")
p1<-ggplot(data=df)+
  theme_bw()+
  geom_point(aes(x=RDA1,y=RDA2, color=Sample),size=3, shape=16)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x="RD1 46.4% ",y="RD2 15.2%",color="Group",fill="Group")+  
  scale_color_manual(values = c("#FF6100","#257d8b","#5861AC")) +
  scale_fill_manual(values = c("#FF6100","#257d8b","#5861AC"))+
  ggtitle("DC1 Axis Lengths = 0.98")+  
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10),
        plot.title = element_text(size=13,hjust = 0.5),
        legend.position = "bottom",
        panel.grid=element_blank())+
  stat_ellipse(data=df, level=0.95,linetype=2, size=0.5,
               aes(x=RDA1,y=RDA2, color=Sample),show.legend = T)
p2 <- p1+ geom_segment(data=df,aes(x=0,y=0,xend=df$env_RDA1,yend=df$env_RDA2),
                       color="#4F4F4F",size=0.55,
                       arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$env_RDA1,y=df$env_RDA2, label=Metal),size=3.5, color="#00008B", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))+
  geom_segment(data=df,aes(x=0,y=0,xend=df$sig_RDA1,yend=df$sig_RDA2),
               color="#FF0000",size=0.55,
               arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$sig_RDA1,y=df$sig_RDA2, label=sig_Metal),size=3.5, color="#8B0000", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))
ggsave("All_pRDA.tiff", p2, width =13, height = 14.5, units = "cm", dpi = 600)
---------
p1<-ggplot(data=df)+
  theme_bw()+
  geom_point(aes(x=RDA1,y=RDA2, color=Sample),size=3, shape=16)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x="RD1 31.4% ",y="RD2 15.0%",color="Group",fill="Group")+  
  scale_color_manual(values = c("#FF6100","#257d8b","#5861AC")) +
  scale_fill_manual(values = c("#FF6100","#257d8b","#5861AC"))+
  ggtitle("DC1 Axis Lengths = 1.24")+  
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10),
        plot.title = element_text(size=13,hjust = 0.5),
        legend.position = "bottom",
        panel.grid=element_blank())+
  stat_ellipse(data=df, level=0.95,linetype=2, size=0.5,
               aes(x=RDA1,y=RDA2, color=Sample),show.legend = T)
p2 <- p1+ geom_segment(data=df,aes(x=0,y=0,xend=df$env_RDA1,yend=df$env_RDA2),
                       color="#4F4F4F",size=0.55,
                       arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$env_RDA1,y=df$env_RDA2, label=Metal),size=3.5, color="#00008B", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))+
  geom_segment(data=df,aes(x=0,y=0,xend=df$sig_RDA1,yend=df$sig_RDA2),
               color="#FF0000",size=0.55,
               arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$sig_RDA1,y=df$sig_RDA2, label=sig_Metal),size=3.5, color="#8B0000", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))
ggsave("All_sRDA.tiff", p2, width =13, height = 14.5, units = "cm", dpi = 600)

#F2.B
df <- read_excel("path/RDA.xls",sheet = "RDA_pindoor")
df <- read_excel("path/RDA.xls",sheet = "RDA_sindoor")
p1<-ggplot(data=df)+
  theme_bw()+
  geom_point(aes(x=RDA1,y=RDA2, color=Group),size=3, shape=16)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x="RD1 39.6% ",y="RD2 11.6%",color="Group",fill="Group")+  
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080")) +
  scale_fill_manual(values = c('#A5C2E2',"#f28080","#FFc080"))+
  ggtitle("DC1 Axis Lengths = 0.98")+  
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10),
        plot.title = element_text(size=13,hjust = 0.5),
        legend.position = "bottom",
        panel.grid=element_blank())+
  stat_ellipse(data=df, level=0.95,linetype=2, size=0.5,
               aes(x=RDA1,y=RDA2, color=Group),show.legend = T)
p2 <- p1+ geom_segment(data=df,aes(x=0,y=0,xend=df$env_RDA1,yend=df$env_RDA2),
                       color="#4F4F4F",size=0.55,
                       arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$env_RDA1,y=df$env_RDA2, label=Metal),size=3.5, color="#00008B", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))+
  geom_segment(data=df,aes(x=0,y=0,xend=df$sig_RDA1,yend=df$sig_RDA2),
               color="#FF0000",size=0.55,
               arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$sig_RDA1,y=df$sig_RDA2, label=sig_Metal),size=3.5, color="#8B0000", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))
ggsave("Indoor_pRDA.tiff", p2, width =13, height = 15, units = "cm", dpi = 600)
---------
p1<-ggplot(data=df)+
  theme_bw()+
  geom_point(aes(x=RDA1,y=RDA2, color=Group),size=3, shape=16)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x="RD1 30.9% ",y="RD2 13.5%",color="Group",fill="Group")+  
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080")) +
  scale_fill_manual(values = c('#A5C2E2',"#f28080","#FFc080"))+
  ggtitle("DC1 Axis Lengths = 1.24")+  
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=10),
        plot.title = element_text(size=13,hjust = 0.5),
        legend.position = "bottom",
        panel.grid=element_blank())+
  stat_ellipse(data=df, level=0.95,linetype=2, size=0.5,
               aes(x=RDA1,y=RDA2, color=Group),show.legend = T)
p2 <- p1+ geom_segment(data=df,aes(x=0,y=0,xend=df$env_RDA1,yend=df$env_RDA2),
                       color="#4F4F4F",size=0.55,
                       arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$env_RDA1,y=df$env_RDA2, label=Metal),size=3.5, color="#00008B", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))+
  geom_segment(data=df,aes(x=0,y=0,xend=df$sig_RDA1,yend=df$sig_RDA2),
               color="#FF0000",size=0.55,
               arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df,aes(x=df$sig_RDA1,y=df$sig_RDA2, label=sig_Metal),size=3.5, color="#8B0000", 
            hjust="inward",vjust=0.5*(1-sign(df$env_RDA2)))
ggsave("Indoor_sRDA.tiff", p2, width =13, height = 15, units = "cm", dpi = 600)

#F2.C and D
##shap_xgboost is the result file from SHAP analysis, the SHAP code showed in Analysis.R
p1 <- sv_importance(shap_xgboost,kind = "beeswarm",show_numbers= TRUE, size=3,
                    format_fun = format_max,
                    number_size = 3.2,
                    sort_features = TRUE)+
  xlab("SHAP value(impact on model output)")+
  theme_bw()+
  theme(legend.title =element_text(size =12),legend.text = element_text(size =10),
        axis.text =element_text(size =12), axis.title =element_text(size =15))+
  scale_color_gradient(low="#049AF5",high="#F9035A")
ggsave("SHAP.tiff", p1, width = 11.2, height = 14, units = "cm", dpi = 600)

#F2.E and F
library(corrplot)
library(WGCNA)
options(scipen=999) 
all.otu <- read_excel('path/Cor.res.xls',sheet = "all_sample")
all.otu <- column_to_rownames(all.otu, "ID") 
Cor.r <- cor(all.otu,all.env,method="spearman")
Cor.p <- corPvalueStudent(Cor.r, 61)
tiff(filename="DEotu_Cor.tiff",width = 25, height = 14, units = "cm", res = 600)
corrplot(Cor.r, method = c('circle'),diag = TRUE, 
         col =c("#049AF5","#168EE9","#2982DD","#3C77D1","#665DB7",
                "#B72B83","#C82179","#D31A71","#E60E65","#F9035A"),
         tl.cex = 0.9, 
         tl.col = 'black', 
         tl.pos = 'TP',
         tl.srt=0, 
         p.mat = Cor.p,
         sig.level = c(.001, .01, .05), insig = "label_sig",
         pch.cex = 1.5, pch.col = 'white')
dev.off() 
--------
indoor.otu <- read_excel('path/Cor.res.xls',sheet = "indoor_sample")
indoor.otu <- column_to_rownames(indoor.otu, "ID") 
indoor.otu <- as.data.frame(t(indoor.otu))
Cor.r <- cor(indoor.otu,indoor.env,method="spearman")
Cor.p <- corPvalueStudent(Cor.r, 41)
tiff(filename="indoor_DEotu_Cor.tiff",width = 25, height = 14, units = "cm", res = 600)
corrplot(Cor.r, method = c('circle'),diag = TRUE, 
         col =c("#049AF5","#168EE9","#2982DD","#3C77D1","#665DB7",
                "#B72B83","#C82179","#D31A71","#E60E65","#F9035A"),
         tl.cex = 0.9, 
         tl.col = 'black', 
         tl.pos = 'TP',
         tl.srt=0, 
         p.mat = Cor.p,
         sig.level = c(.001, .01, .05), insig = "label_sig",
         pch.cex = 1.5, pch.col = 'white')
dev.off() 

#F2.G
##heatmap showed in F2.G.xlsx
site <-read_excel("path/RLQ.res.xls",sheet = "site")
otu <-read_excel("path/RLQ.res.xls",sheet = "otu")
p1 <- ggplot() +
  geom_point(data = otu, aes(AxcQ1, AxcQ2), color="#E6E6FA",size=2, shape=5) +  
  geom_point(data = otu, aes(Q1, Q2), color="#CDB5CD",size=2, shape=5) +
  geom_point(data = site, aes(AxcR1, AxcR2, fill=group),color="white",size=4, shape=21) +
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x = "RLQ1 70%", y = "RLQ2 26%",color="Group") +
  scale_fill_manual(values = c('#A5C2E2',"#f28080","#FFc080")) +
  theme_bw() + 
  theme(panel.grid = element_blank())+
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=12),
        legend.title =element_text(size =15),legend.text = element_text(size =12),
        legend.position = "bottom",
        panel.grid=element_blank())
ggsave("site.tiff", p1, width =13, height = 14.5, units = "cm", dpi = 600) 

#F2.H
site <-read_excel("path/RLQ.res.xls",sheet = "site")
chem <-read_excel("path/RLQ.res.xls",sheet = "chem")
trait <-read_excel("path/RLQ.res.xls",sheet = "trait")
p1<-ggplot(data=site)+
  theme_bw()+
  geom_point(aes(x=AxcR1,y=AxcR2, color=group),size=4, shape=16,alpha=0.5)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x = "RLQ1 70%", y = "RLQ2 26%",color="Group")+
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080")) +
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=12),
        legend.title =element_text(size =15),legend.text = element_text(size =12),
        legend.position = "bottom",
        panel.grid=element_blank())
p2 <- p1+ geom_segment(data=chem,aes(x=0,y=0,xend = chem$RS1, yend = chem$RS2),
                       color="#00008B",size=0.55,
                       arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=chem,aes(x = chem$RS1, y = chem$RS2, label=name),size=3.5, color="#00008B", 
            hjust="inward",vjust=0.5*(1-sign(chem$RS2)))
ggsave("chem.tiff", p2, width =13, height = 14.5, units = "cm", dpi = 600) 
--------
p1<-ggplot(data=site)+
  theme_bw()+
  geom_point(aes(x=AxcR1,y=AxcR2, color=group),size=4, shape=16,alpha=0.5)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x = "RLQ1 70%", y = "RLQ2 26%",color="Group")+
  scale_color_manual(values = c('#A5C2E2',"#f28080","#FFc080")) +
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12,angle=90),
        axis.text=element_text(size=12),
        legend.title =element_text(size =15),legend.text = element_text(size =12),
        legend.position = "bottom",
        panel.grid=element_blank())
p2 <- p1+ geom_segment(data=trait,aes(x=0,y=0,xend = trait$CS1, yend = trait$CS2),
                       color="#CD0000",size=0.55,
                       arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=trait,aes(x = trait$CS1, y = trait$CS2, label=name),size=3.5, color="#CD0000", 
            hjust="inward",vjust=0.5*(1-sign(trait$CS2)))
ggsave("trait.tiff", p2, width =13, height = 14.5, units = "cm", dpi = 600)

#F2.I
##output from Fourth-corner analysis, and the code is showed in Analysis.R

#F2.J
group <- read_excel('path/RLQ.res.xls')[,4:5]
Pro_Y <- cbind(Pro_Y,group)
p1 <- ggplot(data=Pro_Y) +
  geom_segment(aes(x = X1, y = X2,
                   xend = (MDS1), yend = (MDS2),),
               arrow = arrow(length = unit(0, 'cm')),
               color = 'gray', size = 0.8) +
  geom_point(aes(X1, X2, color = Group2), size = 3, shape = 17) +
  geom_point(aes(MDS1, MDS2, color = group), size = 3, shape = 16) +
  scale_color_manual(values = c('#A5C2E2','#A5C2E2',"#f28080","#f28080","#FFc080","#FFc080"))+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.title =element_text(size =15),legend.text = element_text(size =12),
        axis.text =element_text(size =12), axis.title =element_text(size =12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0,  linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[1,2]/Pro_X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[2,2]/Pro_X[2,1], size = 0.3) 
ggsave("pro.tiff", p1, width =18, height = 12, units = "cm", dpi = 600) 
