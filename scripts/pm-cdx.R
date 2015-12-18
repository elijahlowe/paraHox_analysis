source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
a
library("DESeq2")

setwd("/Users/elijahklowe/Desktop/collaborations/claudia/express_counts/sp/overlap/gene_level/")
samples <- c("PM_CONTR_CDX_N8_ACAGTG_L001", "PM_CONTR_CDX_N9_GCCAAT_L001",  "PM_CONTR_CDX_N10_CTTGTA_L001",
             "PM_CDX_N8_ATCACG_L001", "PM_CDX_N9_ACTTGA_L001", "PM_CDX_N10_AGTCAA_L002"
             )
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
  batch=c("N8", "N9", "N10","N8", "N9", "N10"),
  condition=c(rep("ctrl", 3), rep("cdx", 3)),
  libType=rep("paired-end", 6)
)
#Double check it...
pm_cdx_lox.design

#DESeq data object
(Pm.deseq.data <- DESeqDataSetFromMatrix(countData = pm.all.data, 
                                         colData = pm_cdx_lox.design,
                                         design = ~ batch+condition))

Pm.deseq.data$condition <- relevel(Pm.deseq.data$condition, "ctrl")

dds.pm <- DESeq(Pm.deseq.data,test="LRT", reduced=~batch)
plotDispEsts(dds.pm)
resultsNames(dds.pm)
res.pm<-results(dds.pm)
cdx_ctrl.pm <- results( dds.pm, name = "condition_cdx_vs_ctrl", test="Wald")
summary(cdx_ctrl.pm, alpha=0.05)
rld <- rlog(dds.pm)
plotPCA( rld, intgroup = c("condition","batch"), title(main="PCA CdxMO sample 3"))

unknown<-c("PMI_011495","PMI_019393","PMI_026189")
solute_cdx<-c("PMI_014343","PMI_001510","PMI_007866")
digestive_cdx<-c("PMI_024140","PMI_016340","PMI_015084")
dna_binding_cdx<-c("PMI_004263","PMI_000693","PMI_010356")
RA_path_cdx<-c("PMI_008803")
stomach_tdg<-c("PMI_009038")


plotMA(cdx_ctrl.pm, ylim=c(-5,5),lwd=1,cex=0.75,alpha=0.1,main="PmCdx MO 90hpf")
with(cdx_ctrl.pm["PMI_028602", ], {
  points(baseMean, log2FoldChange, col="green", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Hox11/13b", pos=4, col="purple",lwd=6)
})
with(cdx_ctrl.pm[unknown, ], {
  points(baseMean, log2FoldChange, col="pink", cex=1.5, lwd=2, pch=18)
})
with(cdx_ctrl.pm[solute_cdx, ], {
  points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=2, pch=18)
})
with(cdx_ctrl.pm[digestive_cdx, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2, pch=15)
})
with(cdx_ctrl.pm[dna_binding_cdx, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2, pch=17)
})
with(cdx_ctrl.pm[RA_path_cdx, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2, pch=20)
})
with(cdx_ctrl.pm[stomach_tdg, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2, pch="+")
})

cdx_ctrl_reduced.resOrdered <- cdx_ctrl.pm[order(cdx_ctrl.pm$padj),]
cdx_ctrl_reduced.resSig <- subset(cdx_ctrl_reduced.resOrdered, padj <= 0.05 & abs(log2FoldChange)>=0.05)
write.csv(as.data.frame(cdx_ctrl_reduced.resSig),file="Pm-cdx_ctrl.adjp05l2fc05.csv", quote=FALSE)

library("ggplot2")

data <- plotCounts(dds.pm, "PMI_028602-tr", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=4) +
  ggtitle('Hox11/13b: PMI_028777') +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

library( "genefilter" )
library("gplots")
library("RColorBrewer")
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 500 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column",
           cexCol = 1.8,
           labCol = c("CN8", "CN9", "CN10","N8", "N9", "N10"),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( ctrl="gray", lox="darkgreen", cdx = "orange" )[
             colData(rld)$condition ] )
