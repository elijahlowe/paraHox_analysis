source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
a
library("DESeq2")

setwd("/Users/elijahklowe/Desktop/collaborations/claudia/express_counts/sp/overlap/gene_level/")
samples <- c("CONTR_PM_LOX_N2_5_AGTTCC_L002", "CONTR_PM_LOX_N3_4_CCGTCC_L002", #"PM_CONT_LOX_N8_CGATGT_L001",
             "PM_LOX_N2_5_ATGTCA_L002", "PM_LOX_N3_4_GTCCGC_L002") # "PM_LOX_N8_TGACCA_L001")
read.sample <- function(sample.name) {
  file.name <- paste(sample.name, "_htseq_counts.txt", sep="")
  result <- read.delim(file.name, col.names=c("gene","count"),sep="\t", colClasses=c("character", "numeric"))
}

PmCdxCtrl.1 <- read.sample(samples[1])
head(PmCdxCtrl.1)
nrow(PmCdxCtrl.1)

#Read the second sample
PmCdxCtrl.2 <- read.sample(samples[2])

#Let's make sure the first and second samples have the same number of rows and the same genes in each row
nrow(PmCdxCtrl.1) == nrow(PmCdxCtrl.2)
all(PmCdxCtrl.1$gene == PmCdxCtrl.2$gene)

#Now let's combine them all into one dataset
pm.pm.all.data <- PmCdxCtrl.1
pm.all.data <- cbind(PmCdxCtrl.1, PmCdxCtrl.2$count)
for (c in 3:length(samples)) {
  temp.data <- read.sample(samples[c])
  pm.all.data <- cbind(pm.all.data, temp.data$count)
}

#We now have a data frame with all the data in it:
head(pm.all.data)

colnames(pm.all.data)[2:ncol(pm.all.data)] <- samples

#Now look:
head(pm.all.data)

#Let's look at the bottom of the data table
tail(pm.all.data)

pm.all.data <- pm.all.data[1:(nrow(pm.all.data)-5),]
tail(pm.all.data)

rownames(pm.all.data) <-pm.all.data$gene
pm.all.data <- pm.all.data[,2:ncol(pm.all.data)]

pm_cdx_lox.design <- data.frame(
  row.names=samples,
  batch=rep( c("N2", "N3", "N2", "N3")),
  condition=c(rep("ctrl", 2),rep("lox", 2)),
  libType=rep("paired-end", 4)
)
#Double check it...
pm_cdx_lox.design

#cdxMF.design$condition <- relevel(cdxMF.design$condition, ref="ctrl")

#DESeq data object
(Pm.deseq.data <- DESeqDataSetFromMatrix(countData = pm.all.data, 
                                        colData = pm_cdx_lox.design,
                                        design = ~ batch+condition))

#design(Pm.deseq.data) <- formula(~ condition + batch)
Pm.deseq.data$condition <- relevel(Pm.deseq.data$condition, "ctrl")

dds.pm <- DESeq(Pm.deseq.data,test="LRT", reduced=~condition)
plotDispEsts(dds.pm)
resultsNames(dds.pm)
res.pm<-results(dds.pm)
rld <- rlog(dds.pm)
plotPCA( rld, intgroup = c("condition","batch"), title(main="PCA CdxMO sample 3"))

lox_ctrl.pm <- results( dds.pm, name = "condition_lox_vs_ctrl", test="Wald")

solure<-c("PMI_013730","PMI_001719")
digestive<-c("PMI_018190","PMI_004527","PMI_004764","PMI_010194","PMI_000507")
dna_binding<-c("PMI_004511","PMI_004263")
RA_path<-c("PMI_023031","PMI_000148")
fos<-("PMI_015901")

plotMA(lox_ctrl.pm, ylim=c(-5,5),lwd=1,cex=0.75,alpha=0.05,main="PmLox MO 66hpf")

with(lox_ctrl.pm[solure, ], {
points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2, pch=18)
})
with(lox_ctrl.pm[digestive, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2, pch=15)
})
with(lox_ctrl.pm[dna_binding, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2, pch=17)
})
with(lox_ctrl.pm[RA_path, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2, pch=20)
})
with(lox_ctrl.pm["PMI_022837", ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2, pch="+")
})


summary(lox_ctrl.pm, alpha=0.05)

sum(cdx_ctrl.pm$padj<=0.01 & (cdx_ctrl.pm$log2FoldChange>=1), na.rm = TRUE)
sum(cdx_ctrl.pm$padj<=0.01 & (cdx_ctrl.pm$log2FoldChange<=-1), na.rm = TRUE)
sum(cdx_ctrl.pm$padj<0.01 & abs(cdx_ctrl.pm$log2FoldChange)>=1, na.rm = TRUE)

Pm.lox_ctrl.resOrdered <- lox_ctrl.pm[order(lox_ctrl.pm$padj),]
Pm.lox_ctrl.resSig <- subset(Pm.lox_ctrl.resOrdered, padj <= 0.05 & abs(log2FoldChange)>=0.5)
write.csv(as.data.frame(Pm.lox_ctrl.resSig),file="Pm-lox_ctrl.adjp05l2fc05.csv", quote=FALSE)


library( "genefilter" )
library("gplots")
library("RColorBrewer")
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 500 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column",
           cexCol = 1.8,
           labCol = c("CN2", "CN3","CN8","CN8","N2", "N3", "N8"),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( ctrl="gray", lox="darkgreen", cdx = "orange" )[
             colData(rld)$condition ] )

library(ggplot2)

plotMA(cdx_ctrl.pm, ylim=c(-5,5))
with(cdx_ctrl.pm["PMI_029444-tr", ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Hox11/13a", pos=4, col="purple",lwd=6)
})
with(cdx_ctrl.pm["PMI_028602-tr", ], {
  points(baseMean, log2FoldChange, col="green", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Hox11/13b", pos=4, col="purple",lwd=6)
})

(data <- plotPCA(rld_reduced, intgroup = c( "batch","condition"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
qplot(PC1, PC2, color=batch, shape=condition, data=data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

#####################################################
############## Examine known genes ##################
#####################################################

library("ggplot2")

data <- plotCounts(dds.pm, "PMI_020722-tr", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=4) +
  ggtitle('Hox11/13a: PMI_020722') +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

data <- plotCounts(dds.pm, "PMI_029444-tr", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=4) +
  ggtitle('Hox11/13a: PMI_029444') +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

data <- plotCounts(dds.pm, "PMI_028602-tr", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=4) +
  ggtitle('Hox11/13b: PMI_028602') + 
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

data <- plotCounts(dds.pm, "PMI_028777-tr", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=4) +
  ggtitle('Hox11/13b: PMI_028777') + 
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

#####################################################
############## Create HTML reports ##################
##################################################### 
library("ReportingTools")
des2Report <- HTMLReport(shortName = "cdxMOvscontrolR",
                         title = "DESeq2 RNA-seq analysis of CdxMO vs Control reduced",
                         reportDirectory = "./reports")
publish(ddsMF_reduced,des2Report, pvalueCutoff=0.1,contrast = c("condition", "cdx", "ctrl"),
        factor = colData(Pm.deseq.data)$condition,
        reportDir="./reports")
finish(des2Report)

des2Report <- HTMLReport(shortName = "cdxMOvsfluoR",
                         title = "DESeq2 RNA-seq analysis of CdxMO vs Fluo reduced",
                         reportDirectory = "./reports")
publish(ddsMF_reduced,des2Report, pvalueCutoff=0.1,contrast = c("condition", "cdx", "fluo"),
        factor = colData(Pm.deseq.data)$condition,
        reportDir="./reports")
finish(des2Report)

des2Report <- HTMLReport(shortName = "fluovscontrolR",
                         title = "DESeq2 RNA-seq analysis of Fluo vs Control reduced",
                         reportDirectory = "./reports")
publish(ddsMF_reduced,des2Report, pvalueCutoff=0.1,contrast = c("condition", "fluo", "ctrl"),
        factor = colData(Pm.deseq.data)$condition,
        reportDir="./reports")
finish(des2Report)

