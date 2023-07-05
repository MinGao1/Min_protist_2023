library(RColorBrewer)
library("ggpubr")
library(ggplot2)
library(ggthemes)
library(grid)
library(magrittr)
library(edgeR)
library("gplots")
library("limma")
library(dplyr) 
library(tidyr) 
library(splines) 
library(scales)

setwd("")
ma=read.csv(file=".csv",header=T,check.names=FALSE,sep=",")
ma$nich=factor(ma$nich,levels=c())
p2=ggplot(ma, aes(x=, y=, fill=nich,col=))+ 
  geom_boxplot(width=0.5,position=position_dodge(0.6),alpha=0.8)+geom_jitter(width = 0.2,alpha=0.8,size=1)

p2=p2+theme(legend.title = element_text(colour="black", size=9, face="bold"))
p2=p2+theme(legend.text = element_text(colour="black", size = 9))
p2=p2+scale_color_manual(name="Origin",
                         values =c())
p2=p2+scale_fill_manual(name="origin",
                        values =c())
p2=p2+theme(axis.title.x = element_text(face="bold", colour="black", size=9),
            axis.title.y = element_text(face="bold", colour="black", size=9),
            axis.text.x  = element_text(angle=0,size = 9,face="bold",colour="black"),
            axis.text.y  = element_text(angle=0,size =9))
p2=p2+ylab("")+xlab("")
p2






setwd("")
sp=read.csv(file=".csv",header=T,check.names=FALSE,sep=",")
head(sp)
df <- sp %>% gather("item",value,-1) %>% 
bind_cols(data.frame(item_id=rep(:,each=))) 
head(df)
df$item = factor(df$item, levels=c())   
p=ggplot(df,aes(sample,value))+
  geom_bar(aes(fill=item),stat = "identity",position="fill",width=0.8)
p=p+xlab(NULL)+ylab("")
p=p+theme(legend.title=element_blank())
p=p+scale_fill_manual(name=NULL,values = c('gray64',"red",'green4','green2','aquamarine2','magenta1',"blue3","orange",'green4',"cyan1",
                                           "gold","brown",'purple','yellow2'))
p=p+theme(legend.text = element_text(colour="black", size = 9,face="bold"))
p=p+theme(axis.title.x = element_text(face="bold", colour="black", size=9),
          axis.title.y = element_text(face="bold", colour="black", size=9),
          axis.text.x  = element_text(angle=0,size = 9,face="bold",colour="black"),
          axis.text.y  = element_text(angle=0,size = 9))
p=p + scale_y_continuous(labels=percent)
p



options(warn=-1)

rm(list=ls())
# parameters
threshold <- 0.5
scale <- 1000               
alpha <- 0.05              
p.adj.method <- "fdr"       




setwd("")
design = read.table(".txt",row.names= 1, header=T,check.names=FALSE, sep="\t")
head(design)
otu_table <- read.table(file=".txt",header=T,check.names=FALSE,row.names= 1,sep="\t")
head(otu_table)
dim(otu_table)

taxonomy <- read.table(file=".txt",header=F, fill=T,sep="")
head(taxonomy)
dim(taxonomy)

otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))
head(otu_table_norm)
idx <- rowSums(otu_table_norm * 100 >0.1) >= 1
otu_table <- otu_table[idx, ]
head(otu_table)
dim(otu_table)

idx=taxonomy$V1 %in% row.names(otu_table)
taxonomy2=taxonomy[idx,]
head(taxonomy2)
dim(taxonomy2)

idx <- design$HD %in% c("",'')
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]
head(otu_table_subset)
head(design_subset,20)
dim(otu_table_subset)

groups <- design_subset$HD

d <- DGEList(counts=otu_table_subset, group=groups)
d <- calcNormFactors(d,method='TMM')

# fit the GLM
design.mat <- model.matrix(~ 0 + d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
fit <- glmFit(d2, design.mat)

#
lrt_H <- glmLRT(fit, contrast=c(1,-1))
de_H <- decideTestsDGE(lrt_H, adjust.method=p.adj.method, p.value=alpha)
otu_table_subset <- read.table(file="Potutable.txt",header=T,check.names=FALSE,sep="\t")

H_otus <- otu_table_subset$'OTU ID'[de_H==1]
head(H_otus)

H_pvals <- lrt_H$table
head(H_pvals)
H_pvals$sig=de_H
head(H_pvals)
enriched1 = row.names(subset(H_pvals,sig==1))
depleted1 = row.names(subset(H_pvals,sig==-1))
He1=summary(enriched1) 
Hd1=summary(depleted1)
He1
Hd1

head(taxonomy)
idx=taxonomy$V1 %in% row.names(H_pvals)
taxonomy2=taxonomy[idx,]
head(taxonomy2)
dim(taxonomy2)

H_pvals$level = as.factor(ifelse(H_pvals$sig==1, "enriched",ifelse(H_pvals$sig==-1, "depleted","nosig")))
head(H_pvals,10)

