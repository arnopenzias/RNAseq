if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("S4Vectors")
install.packages("png")
BiocManager::install("DESeq2")
################################################          DESeq2 for Riboseq       ####################################################
directory <- "/home/patrick/Riboseq/R/"
library(DESeq2)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

# Create a matrix with the Norm
expression3 <- read.csv(file.path(directory, "expression3.csv"))
head(expression3)
expression3$WT1_norm <- expression3$WT1_sum/expression3$WT1_rpc
expression3$WT2_norm <- expression3$WT2_sum/expression3$WT2_rpc
expression3$D1_norm <- expression3$D1_sum/expression3$D1_rpc
expression3$D2_norm <- expression3$D2_sum/expression3$D2_rpc
expression3$position_name <- paste(expression3$orf, expression3$position, sep="_")
all_data <- select(expression3, D1_norm, D2_norm, WT1_norm, WT2_norm)
row.names(all_data) <- expression3$position_name
head(all_data)
library(tidyr)
all_data$new_col <- row.names(all_data)
all_data <- all_data %>% separate(new_col, c("orf", "position"), sep = "_")
head(all_data)
all_data[] <- sapply(all_data, as.integer)
all_data = na.exclude(all_data)
head(all_data)

# Define your variables
condition <- factor(c('eIF5a','eIF5a','WT','WT'))

# Create the colData object
colData <- data.frame(condition = condition)
row.names(colData) <- c('D1_norm', 'D2_norm', 'WT1_norm', 'WT2_norm')

# Create DEseqset
dds <- DESeqDataSetFromMatrix(countData = all_data,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "eIF5a", "WT"), addMLE = T)


# This block will return a data.table with DESeq output
eIF5_deseq_gene <- as.data.frame(res)
setDT(eIF5_deseq_gene, keep.rownames = T)
colnames(eIF5_deseq_gene)[1] <- "orf"
setkeyv(eIF5_deseq_gene, c("orf"))
View(eIF5_deseq_gene[padj < 0.05])


write.csv(eIF5_deseq_gene, "/Users/danuzarossi/Desktop/Schuler_Green_Kevin/DESeq_gene.csv")


###############################################################################################################################
########################                     EDGE DESeq2           ######################################################
dir <- "/media/patrick/PASIQUIBAC/Bioinfo/RNAseqs/DESeq-EDGE/0/"
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


##########################################################################################################################################################
##################################                    Pearson correlation            #######################################################

#install.packages("corrplot")
#install.packages("PerformanceAnalytics")
#multmerge = function(mypath){
#  filenames=list.files(path=mypath, full.names=TRUE)
#  datalist = lapply(filenames, function(x)
#  {read.delim(file=x)})
#  Reduce(function(x,y) {merge(x,y)}, datalist)}
#mymergeddata = multmerge("/media/patrick/PASIQUIBAC/Bioinfo/Valeria_col/pearson_cor_folder/")

library("robustbase")
library("corrplot")
library("PerformanceAnalytics")
#BiocManager::install("genefilter")
library("genefilter")
directory2 <- "/media/patrick/GERVAZIO/Bioinfo/neuroimuno_lab/NPC_carol/datasets_ago/rnaseq/heatmap/"
#sampleFiles <- grep("*.txt", list.files(directory), value=TRUE)
setwd(directory2)
mydf <- read.delim("dataframe4.txt", row.names = 1)
mydf

mydfwocpt0 <- mydf[-grep("37", names(mydf))]
mydfwocpt <- mydfwocpt0[-grep("HDPC", names(mydfwocpt0))]
mydfwocpt
row_sub= apply(mydf, 1, function(row) all(row !=0 ))
mydf1 <- mydf[row_sub,] #### name modified, old mydf_nonzero
mydf1
write.csv(as.data.frame(mydf1), file="AD_all_nozero.csv")
write_tsv(mydf1, file="AD_all_nozero.txt")

mydf1 <- mydfwocpt+0.00001
mydf1 <- mydf+0.000000001
mydf1 <- mydf

#mydfmean <- mydf_nozero / colMeans(mydf_nozero)[col(mydf_nozero)]
mydf1rowmean0 <- data.frame(#ACBRI371=rowMeans(mydf1[grep("ACBRI371", names(mydf1))]),# HDPC=rowMeans(mydf1[grep("HDPC", names(mydf1))]), 
                           #R186=rowMeans(mydf1[grep("R186", names(mydf1))]), R259=rowMeans(mydf1[grep("R259", names(mydf1))]),
                           #R286=rowMeans(mydf1[grep("R286", names(mydf1))]), 
                           #U87=rowMeans(mydf1[grep("U87", names(mydf1))]),
                           #U138=rowMeans(mydf1[grep("U138", names(mydf1))]),
                           #T98=rowMeans(mydf1[grep("T98", names(mydf1))]),
                           #U251=rowMeans(mydf1[grep("U251", names(mydf1))]),
                           #U343=rowMeans(mydf1[grep("U343", names(mydf1))])#, 
                           #UW467=rowMeans(mydf1[grep("UW467", names(mydf1))])
                           #SL1344=rowMeans(mydf1[grep("SL1344", names(mydf1))]),
                           #STm11=rowMeans(mydf1[grep("STm11", names(mydf1))]),
                           #STm30=rowMeans(mydf1[grep("STm30", names(mydf1))]),
                           #Control=rowMeans(mydf1[grep("Control", names(mydf1))]),
                           #Zika_infected=rowMeans(mydf1[grep("PE243V", names(mydf1))]),
                           Control=rowMeans(mydf1[grep("control", names(mydf1))]),
                           Zika_infected=rowMeans(mydf1[grep("paraiba", names(mydf1))])
                           )
write.csv(as.data.frame(mydf1rowmean0), file="mydf1rowmean0.csv")
mydf1rowmean <- mydf1rowmean0 #+ 0.000000001
row_sub2= apply(mydf1rowmean0, 1, function(row) all(row !=0 ))
mydf1rowmean <- mydf1rowmean0[row_sub2,]


write.csv(as.data.frame(mydf1rowmean), file="row_mean_zika.csv")
mydf1sd <- data.frame(#ACBRI371=rowMeans(mydf1[grep("ACBRI371", names(mydf1))]),# HDPC=rowMeans(mydf1[grep("HDPC", names(mydf1))]), 
  #R186=rowMeans(mydf1[grep("R186", names(mydf1))]), R259=rowMeans(mydf1[grep("R259", names(mydf1))]),
  #R286=rowMeans(mydf1[grep("R286", names(mydf1))]), 
  #U87=rowSds(mydf1[grep("U87", names(mydf1))]),
  #U138=rowSds(mydf1[grep("U138", names(mydf1))]),
  #T98=rowSds(mydf1[grep("T98", names(mydf1))]),
  #U251=rowSds(mydf1[grep("U251", names(mydf1))])#,
  #U343=rowMeans(mydf1[grep("U343", names(mydf1))]), 
  #UW467=rowMeans(mydf1[grep("UW467", names(mydf1))])
  SL1344=rowSds(mydf1[grep("SL1344", names(mydf1))]),
  STm11=rowSds(mydf1[grep("STm11", names(mydf1))]),
  STm30=rowSds(mydf1[grep("STm30", names(mydf1))])
  )


mydfmean <- mydf1rowmean / colMeans(mydf1rowmean)[col(mydf1rowmean)]
mydfmedian <- mydf1rowmean / colMedians(data.matrix(mydf1rowmean))[col(mydf1rowmean)]
mydfmean
lg2mydfmean <- log2(mydfmedian)
lg2mydfmean <- log2(mydfmean)
lg2mydfmean <- log2(mydf1rowmean)
lg2mydfmean <- mydf1rowmean
lg2mydfmean <- mydfmean
lg2mydfmean <- mydfmedian


res <- cor(mydf1)
chart.Correlation(mydf1) #+ xlab("identified transcripts per gene")

library("pheatmap")
library("RColorBrewer")
library("gplots")
#dds_markers <- dds[rownames(dds) %in% tmarkers[1,], ]
#vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
#ntd <- normTransform(dds)
#rld <- rlog(dds, blind=FALSE)

#select <- order(rowMeans(counts(dds, normalized=TRUE)),
#                decreasing=TRUE)[1:4533]
#select <- order(res$pvalue)[1:284]

#df <- as.data.frame(colData(dds)[,"condition"])
#rownames(df) <- colnames(dds)
paletteLength <- 500
markers_rel <- read.delim("relative_values.txt", row.names = 1)
target_genes <- c("AKT1", "AKT2", "AKT3")
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
myBreaks1 <- c(seq(min(lg2mydfmean), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(lg2mydfmean)/paletteLength, max(lg2mydfmean), length.out=floor(paletteLength/2)))
myBreaks2 <- c(seq(min(markers), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(markers)/paletteLength, max(markers), length.out=floor(paletteLength/2)))
markers <- lg2mydfmean[rownames(lg2mydfmean) %in% target_genes, ]
markers2 <- mydf1rowmean[rownames(mydf1rowmean) %in% target_genes, ]
tmarkers2 <- t(t(markers2))
sdmarkers2 <- t(t(mydf1sd[rownames(mydf1sd) %in% target_genes, ]))
pheatmap(markers, cluster_rows = TRUE, show_rownames = TRUE, show_colnames = TRUE,
         cluster_cols = TRUE, annotation_names_col = TRUE, annotation_legend = TRUE, color = colors, breaks = myBreaks2)

barplot2(tmarkers2[3,], plot.ci = TRUE, ci.l = tmarkers2[3,]-sdmarkers2[3,], ci.u = tmarkers2[3,]+sdmarkers2[3,], ylab = "number of transcripts", main = "AKT3")







################################################################################################################################
#########################################################           htseq_count DESeq2             #################################################

directory <- "/media/patrick/GERVAZIO/Bioinfo/neuroimuno_lab/NPC_carol/datasets_ago/rnaseq/htseq/"
setwd(directory)
list.files()
mydf1 <- read.delim("GSE53697_RNAseq_AD.txt")
mydf1 <- mydf1 %>% mutate_if(is.numeric, ~round(., 0))
mydf1
row_sub= apply(mydf1, 1, function(row) all(row !=0 ))
mydf2 <- mydf1[row_sub,]
write_tsv(as.data.frame(mydf2), file="GSE53697_RNAseq_AD_round_nozero.txt")
#setwd("/media/patrick/PASIQUIBAC/Bioinfo/Amanda/RNAseq/ht-seq_bt/37/STm11xSTm30/")
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
dds1 <- estimateSizeFactors(dds)
dds2 <- counts(dds1, normalized = T)


library(ggplot2)
keep <- rowSums(counts(dds)) > 1
dds2 <- dds[keep,]
dds2 <- assay(dds2)
res2 <- results(dds2)
dds3 <- dds2 / colMedians(data.matrix(dds2))[col(dds2)]
dds3 <- dds2 / colMeans(dds2)[col(dds2)]
rld <- log2(dds3) #, blind = TRUE)
rld <- rlog(dds2)
rld_counts <- assay(rld)
head(rld_counts)
plotPCA(rld, intgroup="condition") + geom_text(aes(label=name))

#directory2 <- "/media/patrick/JERONIMO/projects/Doutorado/Pesquisa e laboratório/Projetos paralelos/valeria_col/DESeq2/cpt/"
#setwd(directory2)
#BiocManager::install("EnhancedVolcano")
library("EnhancedVolcano")
#res4 <- read.csv("DESeq2_htseq_ACBRI371_cpt.csv")
markers <- read.delim("imune_genes.txt")
tmarkers <- t(markers)
markers2 <- read.delim("paths.txt")
tmarkers2 <- t(markers2)
EnhancedVolcano(res,
                lab = rownames(res), #NA,
                selectLab = tmarkers[1,],
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.0071,
                #FCcutoff = 0.1
                #xlim = c(-8, 8)
                drawConnectors = TRUE,
                widthConnectors = 0.2
                )



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

#vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
#ntd <- normTransform(dds)
library("vsn")
#meanSdPlot(assay(ntd))
#meanSdPlot(assay(vsd))
#meanSdPlot(assay(rld))

#directory2 <- "/media/patrick/JERONIMO/projects/Doutorado/Artigos e manuscritos/wilson_col/heatmaps_clusters/"
#setwd(directory2)
markers <- read.delim("neut.activ.txt", header = F)
tmarkers <- t(markers)

library("pheatmap")
#dds_markers <- dds[rownames(dds) %in% tmarkers[1,], ]
#vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
#ntd <- normTransform(dds)
#rld <- rlog(dds, blind=FALSE)

#select <- order(rowMeans(counts(dds, normalized=TRUE)),
#                decreasing=TRUE)[1:4533]
#select <- order(res$pvalue)[1:284]

#install.packages("robustbase")
library("robustbase")
library("tidyverse")

directory <- "/media/patrick/GENESIO1/PASIQUIBAC/Bioinfo/Wilson_col/heatmaps/SW480/"
setwd(directory)
markers <- read.delim("target_genes.txt", header = T)
tmarkers <- t(markers)
mydf2 <- read.delim(("dataframe12.txt"), row.names = NULL) ##change to mydf3 if the line below will be needed
#mydf2 <- mydf3 %>% distinct(Gene, .keep_all = TRUE)
mydf2
mydf <- data.frame(mydf2, row.names = 1)
mydf


row_sub= apply(mydf, 1, function(row) all(row !=0 ))
mydf_nozero <- mydf[row_sub,]
mydf_nozero
#mydf1 <- mydf+0.00001
#mydf1 <- data.frame(ID=mydf[,1], mydf[,-1]+0.00001)
mydf1 <- mydf_nozero
#mydfmean <- mydf_nozero / colMeans(mydf_nozero)[col(mydf_nozero)]

mydf1rowmean <- data.frame(POS=rowMeans(mydf1[grep("POS", names(mydf1))]), 
                           SNAIL=rowMeans(mydf1[grep("SNAIL", names(mydf1))]), 
                           PLVX=rowMeans(mydf1[grep("PLVX", names(mydf1))])) #,
                           #NEG=rowMeans(mydf1[grep("NEG", names(mydf1))]))

mydf1rowmean
#mydf1SNAIL <- mydf1[grep("SNAIL", names(mydf1))]
#mydf1SNAILrmean <- data.frame(SNAIL=rowMeans(mydf1SNAIL))
#mydf1POS <- mydf1[grep("POS", names(mydf1))]
#mydf1POSrmean <- data.frame(POS=rowMeans(mydf1POS))
#mydf1NEG <- mydf1[grep("NEG", names(mydf1))]
#mydf1NEGrmean <- data.frame(NEG=rowMeans(mydf1NEG))
#mydf1PLVX <- mydf1[grep("PLVX", names(mydf1))]
#mydf1PLVXrmean <- data.frame(PLVX=rowMeans(mydf1PLVX))

#mydfmediancol <- mydf1rowmean / colMedians(mydf1rowmean)[col(mydf1rowmean)]
#colmeanWT <- colMeans(mydf1rowmean, na.rm = TRUE)
#colmeanWT
mydfmeancol <- mydf1rowmean / colMeans(mydf1rowmean, na.rm = TRUE)[col(mydf1rowmean)]
mydfmeancol
lg2mydfmean <- log2(mydfmeancol)
lg2mydfmean
#lg2mydfmean <- log2(mydf1rowmean)


library("RColorBrewer")
paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength) #brewer.pal(11,"RdBu")
colors <- scico(10, palette = "acton", direction = -1)

#df <- as.data.frame(colData(dds)[,"condition"])
#rownames(df) <- colnames(dds)
dld1_markers <- rld[rownames(rld) %in% tmarkers[1,], ]
dld1_markers <- lg2mydfmean[rownames(lg2mydfmean) %in% tmarkers[9,], ]
myBreaks <- c(seq(min(dld1_markers), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(dld1_markers)/paletteLength, max(dld1_markers), length.out=floor(paletteLength/2)))
colnames(dld1_markers) <- c("Healthy", "Zika infected")
tiff(filename = "heatmap.neut.activ.tiff", units = "px", compression = "lzw", res = 300, height = 3400, width = 1170)
pheatmap(dld1_markers, cluster_rows = TRUE, show_rownames = TRUE, show_colnames = TRUE,
         cluster_cols = F, annotation_names_col = FALSE, annotation_legend = TRUE, main = "Main targets", 
         color = colors)#, breaks = myBreaks)  #colorRampPalette(c("red","white","blue"))(256))
dev.off()


dld1_markers2 <- mydfmean[rownames(mydfmean) %in% tmarkers[1,], ]
a <- dld1_markers2[,1]
b <- dld1_markers2[,2]
c <- a/a
d <- b/a
dld1_markers_rel <- data.frame(c,d)
names(dld1_markers_rel) <- c('Control','Zika infected')
row.names(dld1_markers_rel) <- rownames(dld1_markers2)
dld1_markers_rel
myBreaks <- c(seq(min(dld1_markers_rel), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(dld1_markers_rel)/paletteLength, max(dld1_markers_rel), length.out=floor(paletteLength/2)))
pheatmap(dld1_markers_rel, cluster_rows = TRUE, show_rownames = TRUE, show_colnames = TRUE,
         cluster_cols = FALSE, annotation_names_col = FALSE, annotation_legend = TRUE, main = "Main Targets", 
         color = colors, breaks = myBreaks)  #colorRampPalette(c("red","white","blue"))(256))




sampleDists <- dist(t(assay(vsd)))
#library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(dds_stem_markers)
colnames(sampleDistMatrix) <- colnames(dds_stem_markers)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)




###################################################################################################################
###########################                          Salmon DESeq2           #############################################

dir <- "/media/patrick/GENESIO1/LINUX_inspiron/RNAseqs_old_backup/WTxDygiV/SL1344_as_ref/salmon-fmd_SL1344/"
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
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "SL1344")

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

vsd <- vst(dds2, blind=FALSE)
rld <- rlog(dds2, blind=FALSE)
ntd <- normTransform(dds2)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
library("pheatmap")
#select <- order(rowMeans(counts(dds, normalized=TRUE)),
#                decreasing=TRUE)[1:4533]
#select <- order(res$pvalue)[1:4533]
df <- as.data.frame(colData(dds2)[,"condition"])
rownames(df) <- colnames(dds2)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(dds2)
colnames(sampleDistMatrix) <- colnames(dds2)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)


##################################################################################################################
#################################                       Ballgown             ###############################################

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
directory <- "/media/patrick/PASIQUIBAC/Bioinfo/Valeria_col/stringtie-e/"
sampleFiles <- grep("AST",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*)_AST.*","\\1",sampleFiles)
repnames <- sub(".*AST_(.*)","\\1",sampleFiles)
sampleTable <- data.frame(fileName = sampleFiles,
                          condition = sampleCondition,
                          replicates = repnames)

bg <- ballgown(dataDir = "/media/patrick/PASIQUIBAC/Bioinfo/Valeria_col/stringtie-e/", samplePattern="AST", pData=sampleTable, meas="all")
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
write.csv(results_transcripts, "U87_x_AST_transcripts_results.csv", row.names=FALSE)
#write.csv(results_genes, "chrX_genes_results.csv", row.names=FALSE)

## Filter for genes with q-val <0.05
subset(results_transcripts, results_transcripts$qval <=0.1)
#subset(results_genes, results_genes$qval <=0.05)





##########################################   original ballgown script from paper ###################################################################
#!/usr/bin/env Rscript
# run this in the output directory for rnaseq_pipeline.sh
# passing the pheno data csv file as the only argument 
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # assume no output directory argument was given to rnaseq_pipeline.sh
  pheno_data_file <- paste0(getwd(), "/chrX_data/geuvadis_phenodata.csv")
} else {
  pheno_data_file <- args[1]
}

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

## Read phenotype sample data
pheno_data <- read.csv(pheno_data_file)

## Read in expression data
bg_chrX <- ballgown(dataDir = "ballgown", samplePattern="ERR", pData=pheno_data)

## Filter low abundance genes
bg_chrX_filt <- subset(bg_chrX, "rowVars(texpr(bg_chrX)) > 1", genomesubset=TRUE)

## DE by transcript
results_transcripts <-  stattest(bg_chrX_filt, feature='transcript', covariate='sex', 
                                 adjustvars=c('population'), getFC=TRUE, meas='FPKM')

## DE by gene
results_genes <-  stattest(bg_chrX_filt, feature='gene', covariate='sex', 
                           adjustvars=c('population'), getFC=TRUE, meas='FPKM')

## Add gene name
results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),
                                  geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)

## Write results to CSV
write.csv(results_transcripts, "chrX_transcripts_results.csv", row.names=FALSE)
write.csv(results_genes, "chrX_genes_results.csv", row.names=FALSE)

## Filter for genes with q-val <0.05
subset(results_transcripts, results_transcripts$qval <=0.05)
subset(results_genes, results_genes$qval <=0.05)

## Plotting setup
#tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
#palette(tropical)

## Plotting gene abundance distribution
#fpkm <- texpr(bg_chrX, meas='FPKM')
#fpkm <- log2(fpkm +1)
#boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')

## Plot individual transcripts
#ballgown::transcriptNames(bg_chrX)[12]
#plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
#     main=paste(ballgown::geneNames(bg_chrX)[12], ' : ',ballgown::transcriptNames(bg_chrX)[12]),
#     pch=19, xlab="Sex", ylab='log2(FPKM+1)')
#points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))

## Plot gene of transcript 1729
#plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX,
#                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

## Plot average expression
#plotMeans(ballgown::geneIDs(bg_chrX)[203], bg_chrX_filt, groupvar="sex", legend=FALSE)





##################################################################################################################
###############################                            GOplot           ##################################################

library(GOplot)
setwd("/home/patrick/Riboseq/bubble.plot")
circ <- circle_dat(read.csv('terms.csv'), read.csv('genes.csv'))
write.csv(as.data.frame(circ), file="terms_zscore.csv")
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
GOBubble(reduced_circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), ID = FALSE, labels = 2, table.legend = TRUE, table.col = TRUE, display = 'multiple')
