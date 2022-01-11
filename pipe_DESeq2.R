
###############################################################################################################################
########################                     EDGE DESeq2           ######################################################
dir <- "/path/to/RNAseqs/DESeq-EDGE/"
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
rownames(samples) <- samples$sample
samples[,c("code","sample","condition")]
edge_data <- read.table(file.path(dir, "deseqFile"),row.names = "gene" ,header = TRUE)

library("DESeq2")
dds_edge <- DESeqDataSetFromMatrix(countData = edge_data,
                              colData = samples,
                              design = ~ condition)

dds_edge$condition <- relevel(dds_edge$condition, ref = "SL1344")
dds_edge <- DESeq(dds_edge)

res <- results(dds_edge)
res
dds_edge
resultsNames(dds_edge)

resLFC <- lfcShrink(dds_edge, coef=2)
resLFC

resOrdered <- res[order(res$pvalue),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

resOrdered2 <- res[order(res$log2FoldChange),]


library("IHW")
resIHW <- results(dds_edge, filterFun=ihw)
summary(resIHW)

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

idx <- identify(resLFC$baseMean, resLFC$log2FoldChange)

resApe <- lfcShrink(dds_edge, coef=2, type="apeglm")
resAsh <- lfcShrink(dds_edge, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="normal")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
rownames(resLFC)[idx]

write.csv(as.data.frame(resOrdered), file="DESeq2_results_EDGE.csv")


################################################################################################################################
#########################################################           htseq_count DESeq2             #################################################

directory <- "/path/to/htseq/output/"
setwd(directory)
sampleFiles <- grep("zika",list.files(directory),value=TRUE)
sampleCondition <- sub("zika_(.*)_.*.txt","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition)
ddsHTSeq

ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "control")

dds <- DESeq(ddsHTSeq) #, fitType='local')  
res <- results(dds)
res
dds
resultsNames(dds)

#resLFC <- lfcShrink(dds, coef=2, type="normal")
#resLFC

resOrdered <- res[order(res$pvalue),]
summary(res)
resOrdered

sum(res$padj < 0.1, na.rm=TRUE)

resOrdered2 <- res[order(res$log2FoldChange),]


library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

idx <- identify(resLFC$baseMean, resLFC$log2FoldChange)

resApe <- lfcShrink(dds, coef=2, type="apeglm")
#resAsh <- lfcShrink(dds, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="normal")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
rownames(resLFC)[idx]

write.csv(as.data.frame(resOrdered), file="control_x_zikv.csv")

###################################################################################################################
###########################                          Salmon DESeq2           #############################################

dir <- "/path/to/salmon/output/"
list.files(dir)
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
rownames(samples) <- samples$sample
samples[,c("code","sample","condition")]
files <- file.path(dir, "quant", samples$code, "quant.sf")
names(files) <- samples$sample
tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
head(tx2gene)

library(tximport)
library(readr)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon)
library(DESeq2)

ddsTxi <- DESeqDataSetFromTximport(txi.salmon, colData = samples, design = ~ condition)
ddsTxi
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

dds2 <- DESeq(ddsTxi)
res2 <- results(dds2)
res2

resultsNames(dds2)

resLFC2 <- lfcShrink(dds2, coef=2)
resLFC2

resOrdered3 <- res2[order(res2$pvalue),]
summary(res2)
resOrdered3

sum(res2$padj < 0.1, na.rm=TRUE)

resOrdered4 <- res2[order(res2$log2FoldChange),]

write.csv(as.data.frame(resOrdered3), file="DESeq2_results_salmon.csv")

##################################################################################################################
#################################                       Ballgown             ###############################################

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
directory <- "/path/to/stringtie-e/output/"
sampleFiles <- grep("AST",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*)_AST.*","\\1",sampleFiles)
repnames <- sub(".*AST_(.*)","\\1",sampleFiles)
sampleTable <- data.frame(fileName = sampleFiles,
                          condition = sampleCondition,
                          replicates = repnames)

bg <- ballgown(dataDir = "/path/to/stringtie-e/output/", samplePattern="AST", pData=sampleTable, meas="all")
#structure(bg)$exon
#structure(bg)$intron
#structure(bg)$trans
#transcript_fpkm = texpr(bg, 'FPKM')
#transcript_cov = texpr(bg, 'cov')
#whole_tx_table = texpr(bg, 'all')
#exon_mcov = eexpr(bg, 'mcov')
#junction_rcount = iexpr(bg)
#whole_intron_table = iexpr(bg, 'all')
#gene_expression = gexpr(bg)
#bg@meas

## Filter low abundance genes
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset=TRUE)

## DE by transcript
results_transcripts <-  stattest(bg_filt, feature='transcript', covariate='condition', 
                                 getFC=TRUE, meas='FPKM')

## DE by gene
#results_genes <-  stattest(bg_filt, feature='gene', covariate='condition', 
#                           getFC=TRUE, meas='FPKM')

## Add gene name
results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_filt),
                                  geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)

## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
#results_genes <-  arrange(results_genes, pval)

## Write results to CSV
write.csv(results_transcripts, "transcripts_results.csv", row.names=FALSE)

## Filter for genes with q-val <0.05
subset(results_transcripts, results_transcripts$qval <=0.1)
#subset(results_genes, results_genes$qval <=0.05)

