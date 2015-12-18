source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library("DESeq2")

setwd("/Users/elijahklowe/Desktop/collaborations/claudia/express_counts/sp/overlap/gene_level/")
samples <- c("SpCdx1_ATCACG_L002","SpCdx2_GATCAG_L002","SpCdx3_CCGTCC_L002",
             "SpCdxContr1_TTAGGC_L002", "SpCdxContr2_TAGCTT_L002",  "SpCdxContr3_GTCCGC_L002")
read.sample <- function(sample.name) {
  file.name <- paste(sample.name, "_htseq_counts.txt", sep="")
  result <- read.delim(file.name, col.names=c("gene","count"),sep="\t", colClasses=c("character", "numeric"))
}

SpCdx.1 <- read.sample(samples[1])
head(SpCdx.1)
nrow(SpCdx.1)

#Read the second sample
SpCdx.2 <- read.sample(samples[2])

#Let's make sure the first and second samples have the same number of rows and the same genes in each row
nrow(SpCdx.1) == nrow(SpCdx.2)
all(SpCdx.1$gene == SpCdx.2$gene)

#Now let's combine them all into one dataset
all.data <- SpCdx.1
all.data <- cbind(SpCdx.1, SpCdx.2$count)
for (c in 3:length(samples)) {
  temp.data <- read.sample(samples[c])
  all.data <- cbind(all.data, temp.data$count)
}

#We now have a data frame with all the data in it:
head(all.data)

colnames(all.data)[2:ncol(all.data)] <- samples

#Now look:
head(all.data)

#Let's look at the bottom of the data table
tail(all.data)

all.data <- all.data[1:(nrow(all.data)-5),]

tail(all.data)

#Remove the first column
raw.deseq.data <- all.data[,2:ncol(all.data)]
#Set row names to the gene names
rownames(raw.deseq.data) <- all.data$gene

cdx.design <- data.frame(
  row.names=samples,
  batch=c("sample1","sample2","sample3","sample1","sample2","sample3"),
  condition=rep( c("cdx","cdx","cdx", "ctrl","ctrl","ctrl"), 1 ),
  libType=rep("paired-end", 6)
)
#Double check it...
cdx.design

cdx.design$condition <- relevel(cdx.design$condition, ref="ctrl")
#cdx.design$background <- relevel(cdx.design$background, ref="ORE")

#DESeq data object
(deseq.data <- DESeqDataSetFromMatrix(countData = raw.deseq.data,
                                  colData = cdx.design,
                                  design = ~batch + condition))

ddsMF <- DESeq(deseq.data, reduced=~condition, test="Wald")
ddsMF <- DESeq(deseq.data,test="LRT", reduced=~batch)
plotDispEsts(ddsMF)

rld_reduced <- rlog(ddsMF)
plotPCA( rld_reduced, intgroup = c("condition","batch"), title(main="PCA CdxMO sample 3"))

res <- results( ddsMF )
summary(res, alpha=0.05)
cdx_ctrl.sp <- results( ddsMF, name = "condition_cdx_vs_ctrl", test="Wald")
unknown<-c("WHL22.242989.0","WHL22.699183.1","WHL22.356292.0")
solute_cdx<-c("WHL22.536473.1","WHL22.485270.1","WHL22.402096.0")
digestive_cdx<-c("WHL22.65858.1","WHL22.47793.0","WHL22.494353.0")
dna_binding_cdx<-c("WHL22.502178.0","WHL22.408813.0","WHL22.699375.0")
RA_path_cdx<-c("WHL22.99793.1")
stomach_tdg<-c("WHL22.621623.1")


plotMA(cdx_ctrl.sp, ylim=c(-5,5),lwd=1,cex=0.75,alpha=0.05,main="SpCdx MO 66hpf")
with(cdx_ctrl.sp["WHL22.23530", ], {
  points(baseMean, log2FoldChange, col="blue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, "Wnt10", pos=2, col="blue",cex=1.5)
})
with(cdx_ctrl.sp[unknown, ], {
  points(baseMean, log2FoldChange, col="pink", cex=1.5, lwd=2,pch=18)
})
with(cdx_ctrl.sp[solute_cdx, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1.5, lwd=2,pch=16)
})
with(cdx_ctrl.sp[digestive_cdx, ], {
  points(baseMean, log2FoldChange, col="green", cex=1.5, lwd=2, pch=15)
})
with(cdx_ctrl.sp[dna_binding_cdx, ], {
  points(baseMean, log2FoldChange, col="purple", cex=1.5, lwd=2, pch=17)
})
with(cdx_ctrl.sp[RA_path_cdx, ], {
  points(baseMean, log2FoldChange, col="orange", cex=1.5, lwd=2, pch=20)
})
with(cdx_ctrl.sp[stomach_tdg, ], {
  points(baseMean, log2FoldChange, col="black", cex=1.5, lwd=2, pch="+")
})


cdx_ctrl_reduced.resOrdered <- cdx_ctrl.sp[order(cdx_ctrl.sp$padj),]
cdx_ctrl_reduced.resSig <- subset(cdx_ctrl_reduced.resOrdered, padj <= 0.05 & abs(log2FoldChange)>=0.05)
write.csv(as.data.frame(cdx_ctrl_reduced.resSig),file="Sp-cdx_ctrl.adjp05l2fc05.csv", quote=FALSE)

data <- plotCounts(ddsMF, "WHL22.23530.1", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=4) +
  ggtitle('Wnt-10: WHL22.23530.1') +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

data <- plotCounts(ddsMF, "WHL22.642075.0", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=4) +
  ggtitle('Cdx: WHL22.642075.0') +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

data <- plotCounts(ddsMF, "WHL22.169409.0", intgroup=c("condition","batch"), returnData=TRUE)
ggplot(data, aes(x=condition, y=count, color=batch)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size= 4) +
  ggtitle('lox: WHL22.169409.0') +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
#currently testing
dds.lrt <- DESeq(deseq.data,test = nbinomLRT)
dds.lrt <- estimateSizeFactors(deseq.data)
dds.lrt <- estimateDispersions(dds.lrt)
dds.lrt <- nbinomLRT(dds.lrt,reduced=~ 1)

#Above here works
fit.condition <- fitNbinomGLMs(deseq.data, count ~ condition)
fit.null <- fitNbinomGLMs(deseq.data, count ~ 1)

#Generate raw p-values for the first comparison: full model vs. reduced model without an interaction term; significant p-values will tell you that the full/more complex model does a "significantly" better job at explaining variation in gene expression than the reduced/less complex model (in this case, the one without the interaction term)
pvals.interaction <- nbinomGLMTest(fit.full, fit.nointeraction)
#Generate p-values adjusted for multiple comparisons using the Benjamini-Hochberg approach
padj.interaction <- p.adjust(pvals.interaction, method="BH")
#Look at the genes that have a significant adjusted p-value
fit.full[(padj.interaction <= 0.05) & !is.na(padj.interaction),]

pvals.mutant <- nbinomGLMTest(fit.genotype, fit.null)
padj.interaction <- p.adjust(pvals.interaction, method="BH")
fit.full[(padj.interaction <= 0.05) & !is.na(padj.interaction),]

volcano_mutant <- cbind(pvals.mutant, fit.genotype[,2])
plot(x=volcano_mutant[,2], y=-log10(volcano_mutant[,1]), xlab = "FOLD CHANGE", ylab="-log10(p)")

identify(x=volcano_mutant[,2], y=-log10(volcano_mutant[,1]))
