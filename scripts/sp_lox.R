#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library("DESeq2")

setwd("/Users/elijahklowe/Desktop/collaborations/claudia/express_counts/sp/overlap/gene_level/")
samples <- c("48h_C_ACAGTG_L004", "48h_C_ACAGTG_L005","48h_C_ACAGTG_L006",
             "48h_L_GCCAAT_L004", "48h_L_GCCAAT_L005","48h_L_GCCAAT_L006",
             "72h_C_CTTGTA_L004","72h_C_CTTGTA_L005","72h_C_CTTGTA_L006",
             "72h_L_GTGAAA_L004","72h_L_GTGAAA_L005","72h_L_GTGAAA_L006"
             )
read.sample <- function(sample.name) {
  file.name <- paste(sample.name, "_htseq_counts.txt", sep="")
  result <- read.delim(file.name, col.names=c("gene","count"),sep="\t", colClasses=c("character", "numeric"))
}

SpLoxCtrl.1 <- read.sample(samples[1])
head(SpLoxCtrl.1)
nrow(SpLoxCtrl.1)

#Read the second sample
SpLoxCtrl.2 <- read.sample(samples[2])

#Let's make sure the first and second samples have the same number of rows and the same genes in each row
nrow(SpLoxCtrl.1) == nrow(SpLoxCtrl.2)
all(SpLoxCtrl.1$gene == SpLoxCtrl.2$gene)

#Now let's combine them all into one dataset
SpLox.all.data <- SpLoxCtrl.1
SpLox.all.data <- cbind(SpLoxCtrl.1, SpLoxCtrl.2$count)
for (c in 3:length(samples)) {
  temp.data <- read.sample(samples[c])
  SpLox.all.data <- cbind(SpLox.all.data, temp.data$count)
}

#We now have a data frame with all the data in it:
head(SpLox.all.data)

colnames(SpLox.all.data)[2:ncol(SpLox.all.data)] <- samples

#Now look:
head(SpLox.all.data)

#Let's look at the bottom of the data table
tail(SpLox.all.data)

SpLox.all.data <- SpLox.all.data[1:(nrow(SpLox.all.data)-5),]
tail(SpLox.all.data)

rownames(SpLox.all.data) <- SpLox.all.data$gene
SpLox.all.data <- SpLox.all.data[,2:ncol(SpLox.all.data)]

loxMF_reduce.design <- data.frame(
  row.names=samples,
  batch=c(rep(c("sample1","sample2","sample3"),2),rep(c("sample4","sample5","sample6"),2)), #rep( c("sample1", "sample2", "sample3"),4),
  condition=rep(c("ctrl","ctrl","ctrl","lox","lox","lox"),2),
  time=c(rep("48h",6), rep("72h",6)),
  libType=rep("paired-end", 12)
)

#Double check it...
loxMF_reduce.design

#DESeq data object
(deseqMF.data <- DESeqDataSetFromMatrix(countData = SpLox.all.data,
                                         colData = loxMF_reduce.design,
                                         design = ~ condition+time+condition:time))

full_model <- formula(~ condition+time+condition:time )
deseqMF.data$condition <- relevel(deseqMF.data$condition, "ctrl")

ddsMF_reduced <- DESeq(deseqMF.data, test="LRT",full=full_model, reduced=~time)
plotDispEsts(ddsMF_reduced)

rld <- rlog(ddsMF_reduced)
plotPCA( rld, intgroup = c("condition","batch"), title(main="PCA loxMO sample 3"))

resTC <- results(ddsMF_reduced)
resTC$symbol <- mcols(ddsMF_reduced)$symbol
head(resTC[order(resTC$pvalue),],4)

lox_48h<-results(ddsMF_reduced, name="condition_lox_vs_ctrl", test="Wald")
solute48<-c("WHL22.52290")
digestive48<-c("WHL22.699183","WHL22.531433","WHL22.716369","WHL22.70457","WHL22.131957")
dna_binding48<-c("WHL22.261957","WHL22.699183","WHL22.349501","WHL22.21056","WHL22.699375","WHL22.257532","WHL22.483798","WHL22.194357","WHL22.255599","WHL22.719670","WHL22.450320","WHL22.399764","WHL22.699183","WHL22.699183","WHL22.280477")
RA_path48<-c()
fos<-("WHL22.538480")
plotMA(lox_48h, ylim=c(-5,5),lwd=1,cex=0.75,alpha=0.05, main="SpLox MO 48hpf")
with(lox_48h["WHL22.23530", ], {
  points(baseMean, log2FoldChange, col="blue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Wnt10", pos=2, col="blue",cex=1.5)
})
with(lox_72h[solute48, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1.5, lwd=2,pch=16)
})
with(lox_48h[digestive48, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2, pch=15)
})
with(lox_48h[dna_binding48, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2, pch=17)
})
with(lox_48h[RA_path48, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2, pch=20)
})
with(lox_48h[fos, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2, pch="+")
})

############################################################
######################### Lox 72 ###########################
############################################################

lox_72h<-results(ddsMF_reduced, name="conditionlox.time72h", test="Wald")
solure<-c("WHL22.419541","WHL22.282777")
digestive<-c("WHL22.117270","WHL22.237453","WHL22.100800","WHL22.327805","WHL22.117771")
dna_binding<-c("WHL22.541378","WHL22.502178")
RA_path<-c("WHL22.311466","WHL22.45927")
fos<-("WHL22.538480")

plotMA(lox_72h, ylim=c(-5,5),lwd=1,cex=0.75,alpha=0.05, main="SpLox MO 72hpf")

with(lox_72h[solure, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1.5, lwd=2,pch=16)
})
with(lox_72h[digestive, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2, pch=15)
})
with(lox_72h[dna_binding, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2, pch=17)
})
with(lox_72h[RA_path, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2, pch=20)
})
with(lox_72h[fos, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2, pch="+")
})


with(lox_72h["WHL22.23530", ], {
  points(baseMean, log2FoldChange, col="blue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Wnt10", pos=2, col="blue",cex=1.5)
})
with(lox_72h["WHL22.311466.1", ], {
  points(baseMean, log2FoldChange, col="blue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Rdh8", pos=1, col="blue",cex=1.5)
})


with(lox_48h["WHL22.124787.0", ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Rio3", pos=2, col="dodgerblue",cex=1.5)
})
with(lox_48h["WHL22.349571.0", ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Rdh8_8", pos=2, col="dodgerblue")
})

ctrl_72hv48h<-results(ddsMF_reduced, contrast=c("time","72h","48h"), test="Wald")
ctrl_72hv48h<-results(ddsMF_reduced, name="time_72h_vs_48h", test="Wald")
plotMA(ctrl_72hv48h, ylim=c(-5,5))
with(lox_72h["WHL22.169409", ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "xlox", pos=2, col="dodgerblue")
})


lox_48h.resOrdered <- lox_48h[order(lox_48h$padj),]
lox_48h.resSig <- subset(lox_48h.resOrdered, padj <= 0.05 & abs(log2FoldChange)>=0.5)
write.csv(as.data.frame(lox_48h.resSig),file="Sp-lox_48.adjp05l2fc05.csv", quote=FALSE)

lox_72h.resOrdered <- lox_72h[order(lox_72h$padj),]
lox_72h.resSig <- subset(lox_72h.resOrdered, padj <= 0.05 & abs(log2FoldChange)>=0.5)
write.csv(as.data.frame(lox_72h.resSig),file="Sp-lox_72.adjp05l2fc05.csv", quote=FALSE)

#################################################
##################################################################################################
#################################################

library( "genefilter" )
library("gplots")
library("RColorBrewer")
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( ctrl="gray", lox="darkgreen", cdx = "orange" )[
             colData(rld)$condition ] )
