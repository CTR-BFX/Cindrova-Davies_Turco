## !/usr/local/bin/Rscript
## R 3.6.2
#------------------------------------------------------------------------------
# Menstrual data analysis
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/Cindrova-Davies
#
#
# Analysis Performed by Xiaohui Zhao
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library('DESeq2')
  library('ggplot2')
  library('RColorBrewer')
  library("cowplot")
  library("pheatmap")
  library("ggrepel")
  library("reshape2")
  library("biomaRt")
  library("matrixStats")
  library("plyr")
  library("BiocParallel")
  library("dplyr")
  library("ggalt")
  library("limma")
  library("apeglm")
  library("ComplexHeatmap")
  library("readxl")
})


register(MulticoreParam(2))

Project         <- "CTR_gjb2_0010"
significance    <- 0.05
l2fc            <- 1 
elementTextSize <- 10


message("+-------------------------------------------------------------------------------+")
message("+                  Set up some constants e.g. base directories                  +")
message("+-------------------------------------------------------------------------------+")

Project  <- "CTR_gjb2_0010"
Base.dir <- "/storage/CTR-Projects/CTR_gjb2/CTR_gjb2_0010"
setwd(Base.dir)
HTSeq.dir       <- paste0(Base.dir,"/Menstrual")
out.dir <- paste0(Base.dir,"/Menstrual")
elementTextSize <- 10

message("+-------------------------------------------------------------------------------+")
message("+                       Retrieve ensEMBL annotations                            +")
message("+-------------------------------------------------------------------------------+")

ensembl    =  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'chromosome_name'), mart = ensembl)          
head(ensEMBL2id)

message("+------------------------------------------------------------------------------+")
message("+            Set up the sample table                                           +")
message("+------------------------------------------------------------------------------+")

sampleFiles      <- grep('*htseq_counts.txt',list.files(HTSeq.dir),value=TRUE)

sampleNames      <- gsub("_merged_trimmed_GRCh38.star.bam_htseq_counts.txt", "", sampleFiles)
sampleNames      <- gsub("_R1_001", "", sampleNames)

sampleLabels     <- gsub("_S.*", "", sampleNames)
sampleLabels     <- gsub("_", "-", sampleLabels)

sampleIndividual <- substr(sampleLabels,1,3)

sampleCondition  <- sampleNames
sampleCondition  <- gsub("B55_S.*",  "Scratch", sampleCondition)
sampleCondition  <- gsub("B60_S.*",  "Scratch", sampleCondition)
sampleCondition  <- gsub("B61_S.*",  "Scratch", sampleCondition)
sampleCondition  <- gsub("B70_S.*",  "Scratch", sampleCondition)
sampleCondition  <- gsub("B72_S.*",  "Scratch", sampleCondition)
sampleCondition  <- gsub("B74_S.*",  "Scratch", sampleCondition)
sampleCondition  <- gsub("B75_S.*",  "Scratch", sampleCondition)

sampleCondition  <- gsub("B55M_S.*",  "Menstrual", sampleCondition)
sampleCondition  <- gsub("B60M_S.*",  "Menstrual", sampleCondition)
sampleCondition  <- gsub("B61M_S.*",  "Menstrual", sampleCondition)
sampleCondition  <- gsub("B70M_S.*",  "Menstrual", sampleCondition)
sampleCondition  <- gsub("B72M_S.*",  "Menstrual", sampleCondition)
sampleCondition  <- gsub("B74M_S.*",  "Menstrual", sampleCondition)
sampleCondition  <- gsub("B75M_S.*",  "Menstrual", sampleCondition)


sampleTable <- data.frame(sampleName=sampleNames, fileName=sampleFiles, Labels=sampleLabels, 
                          individual=sampleIndividual, condition=sampleCondition)

print(sampleTable)
nrow(sampleTable)
str(sampleTable)

write.table(sampleTable, file = paste0(out.dir, "/", Project, "-Control_Menstrual_SampleTable.txt"), row.names=F)

message("+-------------------------------------------------------------------------------+")
message("+                  Produce ddHTSeq object                                       +")
message("+-------------------------------------------------------------------------------+")

#
#   Control vs Menstrual
#

ddsHTSeq.menstrual <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, 
                                                 design=~individual+condition)
dds.menstrual <- collapseReplicates(ddsHTSeq.menstrual, ddsHTSeq.menstrual$Labels, renameCols = TRUE)

dds.menstrual <- DESeq(dds.menstrual, parallel=TRUE)
resultsNames(dds.menstrual)
design(dds.menstrual)
colData(dds.menstrual)

rld.menstrual <- rlogTransformation(dds.menstrual, blind=F) 


message("+-------------------------------------------------------------------------------+")
message("+                       Produce  PCA plot                                       +")
message("+-------------------------------------------------------------------------------+")

# PCA plot

customPCA <- function(sampleTBL, RLD, TOPNUM, model, ensEMBL2id) {
  
  rv     <- rowVars(RLD)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=rownames(sampleTBL), pca$x, individual=sampleTBL$individual, 
                          condition=sampleTBL$condition)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    theme(text = element_text(size=elementTextSize)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3 ) + 
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  loadings                 <- as.data.frame(pca$rotation)
  loadings$ensembl_gene_id <- rownames(loadings)
  loadings                 <- merge(loadings, ensEMBL2id, by="ensembl_gene_id")
  
  pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(external_gene_name,levels=unique(external_gene_name)), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(external_gene_name,levels=unique(external_gene_name)), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  return(list(plt.pca, plt.pca.nl, pca.1.25.plot, pca.2.25.plot) )
  
}

sampleTBL <- as.data.frame(colData(dds.menstrual))
RLD <- assay(rld.menstrual)
TOPNUM <- 2000
model <- "rld.menstrual" 


pca.plt.rld.menstrual.100 <- customPCA(sampleTBL, RLD, TOPNUM, model, ensEMBL2id)

pdf(paste0(out.dir, "/", Project, "-PCAplot_rld.menstrual_TopN", TOPNUM, ".pdf"), height = 8, width = 10)
plot_grid(pca.plt.rld.menstrual.100[[1]],
          pca.plt.rld.menstrual.100[[2]],
          pca.plt.rld.menstrual.100[[3]],
          pca.plt.rld.menstrual.100[[4]])
dev.off()



message("+-------------------------------------------------------------------------------+")
message("+                       DESeq Analysis                                          +")
message("+-------------------------------------------------------------------------------+")

res.menstrual <- lfcShrink(dds.menstrual,coef=2, type="apeglm", parallel=TRUE)

topN <- 10;
l2fc <- 1;
significance <- 0.05;

nrow(subset(res.menstrual, padj <= significance & abs(log2FoldChange) >= l2fc))
head(subset(res.menstrual, padj <= significance & abs(log2FoldChange) >= l2fc))

## 103 significants

res.menstrual.dat <- as.data.frame(res.menstrual)
res.menstrual.dat$ensembl_gene_id <- rownames(res.menstrual.dat)

res.menstrual.mer <- merge(res.menstrual.dat, ensEMBL2id, by = "ensembl_gene_id")

res.menstrual.sig <- subset(res.menstrual.mer, abs(log2FoldChange) >= 1)
res.menstrual.sig <- subset(res.menstrual.sig, padj < 0.05)  ## 101

print(dim(res.menstrual.sig))

save(dds.menstrual, rld.menstrual, res.menstrual, file = paste0(out.dir, "/", Project, "-Control_Menstrual_dds_rld_res.RData"))

write.csv(res.menstrual.sig, file = paste0(out.dir, "/", Project, "-Menstrual_Control_sigDEGs_N101_List.csv"), row.names=F, quote=T)

message("+-------------------------------------------------------------------------------+")
message("+  Some individual genes counts plots identified use PCA gene                   +")
message("+-------------------------------------------------------------------------------+")

rv     <- rowVars(RLD)
select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
pca    <- prcomp(t(RLD[select, ]))
loadings                 <- as.data.frame(pca$rotation)
loadings$ensembl_gene_id <- rownames(loadings)
loadings                 <- merge(loadings, ensEMBL2id, by="ensembl_gene_id")
pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.5.gene      <- pca.1[c(1:5),c(1,16)]
pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.5.gene      <- pca.2[c(1:5),c(1,16)]

genestoplot <- c(pca.1.5.gene[,1], pca.2.5.gene[,1])
genesnames <- c(pca.1.5.gene[,2], pca.2.5.gene[,2])

makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, TITLE, gene2plot, outdir) {
  #
  # Plot the normalised read counts for a specified gene
  #
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep=""))
  
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2) <- c("count", "condition")
  
  
   plt <- ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
      geom_boxplot(width = 0.5, color=c("purple4","green"), outlier.shape=NA) + 
      geom_point(position=position_jitter(w=0.1,h=0), alpha=0.5) +
      ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot)) + 
      xlab("") + ylab("Normalised count") +
      scale_fill_manual(name="Genotype", values = c("purple4",  "green")) +
      theme(text = element_text(size=elementTextSize), legend.position="none")+
      theme_classic() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(size=elementTextSize))
  
  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "condition", "samples")
  t2           <- t2[order(t2$condition),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))
  
  plt.grp <- ggplot(t2, aes(x=samples2, y=count, fill=condition, group=condition)) + 
      geom_bar(stat="identity", alpha=.5) + 
      scale_fill_manual(name="Comparison", values = c("purple4", "green")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") +
      ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot))+
      theme_classic() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(size=elementTextSize))
  print(paste("Created plot for",  genename2plot), sep=" ")
  
  return(list(plt, plt.grp))
}
CONDIRION <- "condition"
GROUP <- "Condition"
out.dir <- outdir <- out.dir


for(i in 1:10){
  plts<- makeGeneCountPlot(DDS, ensEMBL2id, CONDITION, GROUP, genestoplot[i], out.dir)
  pdf(paste0(out.dir, "/", Project, "-GeneCount_", genesnames[i],".pdf"), height=6, width=8)
  plot_grid(plts[[1]], plts[[2]], labels = c('A', 'B'), label_size = 12, ncol = 1)
  dev.off() 
}



message("+-------------------------------------------------------------------------------+")
message("+                        MA plot, Volcano plot                                  +")
message("+-------------------------------------------------------------------------------+")

functionCustomMAPlot <- function(results, Project, FigureID, Title, significance, log2FC, topN, col_gender, ylim, yrange) {
  
  results$log2FoldChange[results$log2FoldChange > ylim]    <- ylim
  results$log2FoldChange[results$log2FoldChange < -(ylim)] <- -(ylim)
  
  res_ord_red  <- results[order(results$log2FoldChange,decreasing=TRUE),] 
  res_ord_red  <- subset(res_ord_red, log2FoldChange >= log2FC & padj <= significance)
  
  res_ord_blue <- results[order(results$log2FoldChange,decreasing=FALSE),] 
  res_ord_blue <- subset(res_ord_blue, log2FoldChange <= -(log2FC) & padj <= significance)
  
  plt <- ggplot(data = results, aes(x=baseMean, y=log2FoldChange )) +
    geom_abline(intercept = 0,       slope = 0, colour='black', alpha=0.5) +
    geom_abline(intercept = log2FC,  slope = 0, colour='black', alpha=0.5, linetype="dashed") +
    geom_abline(intercept = -log2FC, slope = 0, colour='black', alpha=0.5, linetype="dashed") +
    geom_point(size=0.5, alpha=0.5, col="black") +
    geom_point(data=subset(results, padj <= significance & log2FoldChange > l2fc), size=1, alpha=0.95,  col="red") +
    geom_point(data=subset(results, padj <= significance & log2FoldChange < -l2fc),  size=1, alpha=0.95,  col="blue") +
    
    geom_point(data=subset(results, results$chromosome_name == "Y" & padj <= significance & col_gender == 1), alpha=0.99, size=1, colour="cyan") +
    geom_point(data=subset(results, results$chromosome_name == "X" & padj <= significance & col_gender == 1), alpha=0.99, size=1, colour="pink") +
    
    geom_label_repel(data=res_ord_red[c(1:topN),], aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                     size=3, segment.size = 0.25, segment.color = 'darkred', force=5, nudge_x = 0, nudge_y=0) +
    geom_label_repel(data=res_ord_blue[c(1:topN),],
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                     size=3, segment.size = 0.25, segment.color = 'darkblue', force=5, nudge_x = 0, nudge_y=0) +
    scale_x_log10( ) +
    scale_y_continuous(breaks=seq(yrange[1],yrange[2],yrange[3])) + 
    xlab("Mean Normalised Read Count") + ylab("log2 Fold Change") + ggtitle(paste(Project, " ", FigureID, "\n", Title, " [fc ", log2FC, ", sig ", significance, "]", sep="")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  return (plt)
}

functionPlotDEVolcano <- function(results, sig_cut, logfc_cut, title,  xrange, yrange, topN, xlabel, ylabel, col_gender) {
  
  results       <- as.data.frame(results)
  results$genes <- results$external_gene_name
  results <- subset(results, !is.na(padj))
  results <- results[order(-results$log2FoldChange),]
  
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log10(padj), label=genes)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    
    geom_point(data=subset(results, abs(log2FoldChange) < logfc_cut | padj > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange >= logfc_cut),      alpha=0.75, size=0.8, colour="red") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=0.8, colour="blue") +
    geom_point(data=subset(results, chromosome_name == "Y" & padj<=sig_cut & col_gender == 1), alpha=0.99, size=0.8, colour="blue") +
    geom_point(data=subset(results, chromosome_name == "X" & padj<=sig_cut & col_gender == 1), alpha=0.99, size=0.8, colour="pink") +
    geom_text_repel( data= subset(results, log2FoldChange > logfc_cut & padj<= sig_cut & external_gene_name!="")[1:topN,],
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= tail(subset(results, log2FoldChange < (-logfc_cut) & padj<= sig_cut & external_gene_name!=""), n=topN),
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) +
    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    theme(aspect.ratio=1) +
    ggtitle(title) +
    theme_update(plot.title = element_text(size=16, face="bold", hjust=0.5),
                 axis.title.x = element_text(size=12, face= "bold"),
                 axis.text.x = element_text(size=12, face="bold"),
                 axis.title.y.left = element_text(size=12, face= "bold"),
                 axis.text.y = element_text(size=12, face="bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  
  
  return(volc.plt)
  
}

topN=10
sig_cut = 0.5
logfc_cut = 1
xrange=c(-12,12,4);
yrange=c(0,40,10)
col_gender=1
title = "Menstrual vs Control"
results = res.menstrual.mer
xlabel <- "log2FC (Menstrual/Control)"
ylabel <- bquote("-log"[10]~"(adj.p.value)")
res.control.menstrual.Vol <- functionPlotDEVolcano (results, sig_cut, logfc_cut, title,  xrange, yrange, topN, xlabel, ylabel, col_gender)
res.control.menstrual.MA <- functionCustomMAPlot(res.menstrual.mer, Project, "MA", "res.Control_Menstrual", significance, l2fc, topN, col_gender=1,  ylim=10, yrange=c(-12,8,4))

## labels the genes with abs(log2fold) >=2 and padj < 0.05
## also the pink genes are the chrX, cyan for the genes from Y. We had pink only. 
pdf(paste0(out.dir, "/CTR_gjb2_0010-Menstrual_Volcano_Plot.pdf"))
#res.control.menstrual.MA
res.control.menstrual.Vol
dev.off()


message("+-------------------------------------------------------------------------------+")
message("+               Heatmap for the 101 significant                                 +")
message("+-------------------------------------------------------------------------------+")


heatmap_mat_fn <- function(rlds, genes2plot, ensE, samTable, del.cols, sel.cols, outfile, width, height){
  mat.rows <- match(genes2plot, row.names(rlds))
  mat.rows <- mat.rows[!is.na(mat.rows)]
  mat1 <- assay(rlds)[mat.rows,]
  mat1.df <- data.frame(ensembl_gene_id=rownames(mat1), mat1)
  mat1.ann <- merge(mat1.df, ensE, by ="ensembl_gene_id")
  mat1.ann$dups <- duplicated(mat1.ann$ensembl_gene_id)
  mat1.ann.dedup <- subset(mat1.ann, mat1.ann$dups == F)
  rownames(mat1.ann.dedup) <-  mat1.ann.dedup$external_gene_name
  mat2 <- mat1.ann.dedup[,-c(del.cols)]
  samTable <- as.data.frame(samTable)
  samTable <- samTable[order(samTable$Labels),]
  rownames(samTable) <- samTable$Labels
  colnames(mat2) <- samTable$Labels
  condition <- samTable[, sel.cols]
  annotation_col <- as.data.frame(condition)
  names(annotation_col) <- "Treatment"
  rownames(annotation_col) <- samTable$Labels
  Treatment <- annotation_col$condition
  ann_colors <- list(Treatment = c(Scratch="lightseagreen", Menstrual="coral1" ))
  pdf(outfile, onefile=FALSE, width=width, height=height)
  par(bg=NA)
  pheatmap(mat2,  annotation_col = annotation_col, annotation_colors = ann_colors,
           fontsize=9, fontsize_row=6.5, show_rownames=T, cluster_cols = T,
           cluster_rows=T,treeheight_row = 0,
           fontsize_col = 9)
  dev.off()
  
}

outHfile <- paste0(out.dir, "/", Project, "-Heatmap_Control_Menstrual_padj0.05_lfc1_top101.pdf")
rlds <- rld.menstrual
genes2plot <- res.menstrual.sig$ensembl_gene_id
ensE <- ensEMBL2id
samTable <- sampleTable
outfile <- outHfile
del.cols <- c(1,16:19)
sel.cols <- 5

width <- 7; height = 13

heatmap_mat_fn(rlds, genes2plot, ensE, samTable, del.cols, sel.cols, outfile, width, height)


message("+--- Select stress markers to show that endometrial and menstrual organoids of the same patient in same level ----+")

smarkers.dat <- res.menstrual.mer[grep("stress", res.menstrual.mer$description),]
hmarkers.dat <- res.menstrual.mer[grep("hypoxia", res.menstrual.mer$description),]
stressmarkers <- unique(c("GDF15", "TNF", "HSPA5", "HSPA9", "AGER", "HSP90AA1", "HSP90AB1", "HSP90B1", "TRAP1", 
                   smarkers.dat$external_gene_name, hmarkers.dat$external_gene_name,
                   "APOE", "SELENOP", "SOD3", "FOXO1", "LEP", "SLC6A4", "PMAIP1", "EGR1", "CYP1A1"))
stressmarkers.dat <- res.menstrual.mer[res.menstrual.mer$external_gene_name%in%stressmarkers==T, ]
stressmarkers.dat$description <- gsub("..Source.*", "", stressmarkers.dat$description)

test <- stressmarkers[stressmarkers%in%res.menstrual.sig$external_gene_name==T]
testall <- stressmarkers[stressmarkers%in%res.menstrual.mer$external_gene_name==T]
## none of the above are in the significant DEGs

stressmarkers.dat <- stressmarkers.dat[order(stressmarkers.dat$description), ]

stressmarkers.dat <- stressmarkers.dat[-c(2,10,51,52,59),]
stressmarkers.dat$stress <- c("Inflammation", "ER stress", rep("Oxidative stress", 6), rep("Inflammation",2), "ER stress",
                              rep("Oxidative stress", 41), "ER stress", "Oxidative stress", rep("ER stress", 4), 
                              "Oxidative stress", rep("Inflammation",2))

stressmarkers.dat <- stressmarkers.dat[order(stressmarkers.dat$stress), ]

message("+---        Heatmap to all samples with these stress markers                        ----+")

functionPlotHeatmap_DataSort <- function(rlds, sel.cols, ensEMBL2id, selGenes, cont.cols){
  
  mat1 <- as.data.frame(assay(rlds)[,sel.cols])
  mat1$ensembl_gene_id <- rownames(mat1)
  mat1.merE <- merge(mat1, ensEMBL2id, by="ensembl_gene_id")
  mat1.sel <- mat1.merE[mat1.merE$external_gene_name%in%selGenes==T,]
  mat1.sHeatmat <- mat1.sel[,cont.cols]
  rownames(mat1.sHeatmat) <- mat1.sel$external_gene_name
  mat1.sHeatmat <- mat1.sHeatmat[order(rownames(mat1.sHeatmat)),]
  mat1.sHeatmat
}

selGenes.stress <- list(stressmarkers.dat[stressmarkers.dat$stress=="ER stress", 7], 
                     stressmarkers.dat[stressmarkers.dat$stress=="Inflammation", 7],
                     stressmarkers.dat[stressmarkers.dat$stress=="Oxidative stress", 7])
sel.cols <- c(1:14)
cont.cols <- c(2:15)
rlds <- rld.menstrual
stress.gomat.all <- list()
for (i in 1:3){
  selGenes <- selGenes.stress[[i]]
  stress.gomat.all[[i]] <- functionPlotHeatmap_DataSort(rlds, sel.cols, ensEMBL2id, selGenes, cont.cols)
  
  stress.gomat.all
}

stress.goHeatmat  <- rbind(stress.gomat.all[[1]], stress.gomat.all[[3]], stress.gomat.all[[2]])

minval <- floor(min(stress.goHeatmat ))-1
maxval <- round(max(stress.goHeatmat ))+1
breaksList.pro = seq(minval, maxval, by = 1)
ScaleCols.pro <- colorRampPalette(colors = c("purple4","white","green"))(length(breaksList.pro))

stress.goHeatmap <- Heatmap(as.matrix(stress.goHeatmat),
                               col = ScaleCols.pro, 
                               name = "Stress Relating", show_row_names=T,
                               show_column_names = T, width = unit(4, "cm"),
                               heatmap_legend_param = list(title = "Expression"),
                               cluster_rows = F,show_row_dend = F,
                               row_order=c(1:61),
                               column_title="Stress Relating",
                               column_order=c(1:14),
                               column_names_gp = gpar( fontsize = 6),
                               row_title_rot = 0,
                               row_gap = unit(8, "mm"),
                               row_names_gp = gpar( fontsize = 8),
                               row_title_gp = gpar(fontsize =10, fontface = "bold"),
                               row_split = rep(c("ER stress", 
                                                 "Oxidative stress",
                                                 "Inflammation"), 
                                               c(7,49,5)) )
pdf(paste0(out.dir, "/", Project, "-Heatmap_endometrial_Menstrual_stress_N61.pdf"))
print(stress.goHeatmap)
dev.off()




message("+---        Heatmap to B75 (EPL)  and B70 (BP) samples with these stress markers      ----+")

stress.goHeatmap.sub <- Heatmap(as.matrix(stress.goHeatmat[,c(7:8,13:14)]),
                                col = ScaleCols.pro, 
                                name = "Stress Relating", show_row_names=T,
                                show_column_names = T, width = unit(4, "cm"),
                                heatmap_legend_param = list(title = "Expression"),
                                cluster_rows = F,show_row_dend = F,
                                row_order=c(1:61),
                                column_title="Stress Relating",
                                column_order=c(1:4),
                                column_names_gp = gpar( fontsize = 6),
                                row_title_rot = 0,
                                row_gap = unit(8, "mm"),
                                row_names_gp = gpar( fontsize = 8),
                                row_title_gp = gpar(fontsize =10, fontface = "bold"),
                                row_split = rep(c("ER stress", 
                                                  "Oxidative stress",
                                                  "Inflammation"), 
                                                c(7,49,5)) )
pdf(paste0(out.dir, "/", Project, "-Heatmap_endometrial_Menstrual_stress_N61_B70_B75.pdf"))
print(stress.goHeatmap.sub)
dev.off()

message("+---- Add heatmap for stem cells as Tereza selected  ----+")

stemcell <- read_excel(paste0(out.dir,"/Menstrual_StemCell_markers.xls"))
colnames(stemcell) <- c("Category", "external_gene_name")
stemcell.mer <- merge(stemcell, res.menstrual.mer, by = "external_gene_name", all.x=T)

categories <- unique(stemcell.mer$Category)

selGenes.stem <- list(stemcell.mer[stemcell.mer$Category==categories[1], 1], 
                        stemcell.mer[stemcell.mer$Category==categories[2], 1],
                        stemcell.mer[stemcell.mer$Category==categories[3], 1],
                        stemcell.mer[stemcell.mer$Category==categories[4], 1],
                        stemcell.mer[stemcell.mer$Category==categories[5], 1],
                        stemcell.mer[stemcell.mer$Category==categories[6], 1],
                        stemcell.mer[stemcell.mer$Category==categories[7], 1])

sel.cols <- c(1:14)
cont.cols <- c(2:15)
rlds <- rld.menstrual
stem.cat.all <- list()
for (i in 1:7){
  selGenes <- selGenes.stem[[i]]
  stem.cat.all[[i]] <- functionPlotHeatmap_DataSort(rlds, sel.cols, ensEMBL2id, selGenes, cont.cols)
  
  stem.cat.all
}

stem.catHeatmat  <- rbind(stem.cat.all[[1]], stem.cat.all[[2]], stem.cat.all[[3]],
                          stem.cat.all[[4]], stem.cat.all[[5]], stem.cat.all[[6]],stem.cat.all[[7]])
stem.len <- unlist(lapply(stem.cat.all, function(x) dim(x)[1]))
stem.rowTot <- sum(stem.len)

minval <- floor(min(stem.catHeatmat ))-1
maxval <- round(max(stem.catHeatmat ))+1
breaksList.stem = seq(minval, maxval, by = 1)
ScaleCols.stem <- colorRampPalette(colors = c("purple4","white","green"))(length(breaksList.stem))

stem.catHeatmap.ori <- Heatmap(as.matrix(stem.catHeatmat),
                            col = ScaleCols.stem, 
                            name = "Stem Cell", show_row_names=T,
                            show_column_names = T, width = unit(4, "cm"),
                            heatmap_legend_param = list(title = "Expression"),
                            cluster_rows = F,show_row_dend = F,
                            cluster_columns=F,
                            row_order=c(1:stem.rowTot),
                            column_title="Stem Cell",
                            column_order=c(1:14),
                            column_names_gp = gpar( fontsize = 6),
                            row_title_rot = 0,
                            row_gap = unit(8, "mm"),
                            row_names_gp = gpar( fontsize = 7, fontface="bold"),
                            row_title_gp = gpar(fontsize =10, fontface = "bold"),
                            row_split = rep(categories, stem.len)
                            )
stem.catHeatmap.ord <- Heatmap(as.matrix(stem.catHeatmat[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]),
                           col = ScaleCols.stem, 
                           show_row_names=T,
                           show_column_names = T, 
                           width = unit(4, "cm"),
                           heatmap_legend_param = list(title = "Expression (rld)", 
                                                       legend_height = unit(3, "cm"), 
                                                       title_position = "leftcenter-rot"),
                           cluster_rows = F,show_row_dend = F,
                           cluster_columns=F,
                           row_order=c(1:stem.rowTot),
                           column_order=c(1:14),
                           column_split = rep(c("Scratch", "Menstrual"), c(7,7)),
                           column_names_gp = gpar(fontsize = 7,fontface = "bold"),
                           row_title_rot = 0,
                           column_title_rot = 0,
                           row_title = c("Cilia markers", "Endometrial gland \ndevelopment \nand function",
                                         "Epithelial cell \nmarkers", "Epithelial cell \nsecretion",
                                         "Hormonal receptors", "Progenitor cell \nmarkers",
                                         "Proliferative markers"),
                           row_gap = unit(8, "mm"),
                           row_names_gp = gpar( fontsize = 7, fontface = "italic"),
                           row_title_gp = gpar(fontsize =8, fontface = "bold"),
                           row_split = rep(categories, stem.len)
                            
)
pdf(paste0(out.dir, "/", Project, "-Heatmap_endometrial_Menstrual_stemcell_N38_Cat7_ord.pdf"), width=)
#print(stem.catHeatmap.ori)
#print(stem.catHeatmap.ord)
draw(stem.catHeatmap.ord, row_title = " ", row_title_gp = gpar(col = "red"),  
     column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

stem.catHeatmap.sub.ori <- Heatmap(as.matrix(stem.catHeatmat[,c(7:8,13:14)]),
                               col = ScaleCols.stem, 
                               name = "Stem Cell", show_row_names=T,
                               show_column_names = T, width = unit(4, "cm"),
                               heatmap_legend_param = list(title = "Expression"),
                               cluster_rows = F,show_row_dend = F,
                               cluster_columns=F,
                               row_order=c(1:stem.rowTot),
                               column_title="Stem Cell",
                               column_order=c(1:4),
                               column_names_gp = gpar( fontsize = 6),
                               column_split = rep(c("BP", "EPL"), c(2,2)),
                               row_title_rot = 0,
                               row_gap = unit(8, "mm"),
                               row_names_gp = gpar( fontsize = 8),
                               row_title_gp = gpar(fontsize =10, fontface = "bold"),
                               row_split = rep(categories, stem.len)
)
stem.catHeatmap.sub.ord <- Heatmap(as.matrix(stem.catHeatmat[,c(7,13,8,14)]),
                               col = ScaleCols.stem, 
                               name = "Stem Cell", 
                               show_row_names=T,
                               show_column_names = T, 
                               width = unit(4, "cm"),
                               heatmap_legend_param = list(title = "Expression (rld)", 
                                                           legend_height = unit(4, "cm"), 
                                                           title_position = "leftcenter-rot"),
                               cluster_rows = F,show_row_dend = F,
                               cluster_columns=F,
                               row_order=c(1:stem.rowTot),
                               column_title="Stem Cell",
                               column_order=c(1:4),
                               column_split = rep(c("Endometrial", "Menstrual"), c(2,2)),
                               column_names_gp = gpar( fontsize = 6),
                               row_title_rot = 0,
                               row_gap = unit(8, "mm"),
                               row_names_gp = gpar( fontsize = 8),
                               row_title_gp = gpar(fontsize =10, fontface = "bold"),
                               row_split = rep(categories, stem.len)
                              )


pdf(paste0(out.dir, "/", Project, "-Heatmap_endometrial_Menstrual_stemcell_N38_Cat7_B70_B75.pdf"))
#print(stem.catHeatmap.sub.ori)
draw(stem.catHeatmap.sub.ord, row_title = " ", row_title_gp = gpar(col = "red"),  
     column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()






message("+ ----   END OF SCRIPT                                                         ------ +")
