# Author: Cassie Litchfield
library(edgeR)
library(matrixStats)
library(RColorBrewer)
library(openxlsx)
library(ruv)
library(EDASeq)

indir <- "/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/analyses/"
outdir <- "/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/analyses/"
htseq_path <- "/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/ht_seq_counts/"

# -------------------------------------------------------------------------

# Read sample annotation
annot <- read.xlsx(paste0(indir,"metadata_v2.xlsx"), cols = 1:6)
annot$sampleID <- paste0("SampleID.",annot$SampleNo)
annot$sampleID <- gsub(" \\(",".",annot$sampleID)
annot$sampleID <- gsub("\\)","",annot$sampleID)

# Read htseq counts
files <- list.files(htseq_path,".txt", full.names = T)
raw.counts <- readDGE(files, columns=c(1,2),header=FALSE)
metaTags <- grep("^__", rownames(raw.counts))
count_stats <- raw.counts[metaTags, ] #rm metatags
raw.counts <- raw.counts[-metaTags, ] #rm metatags
dim(raw.counts)
colnames(raw.counts) <- basename(colnames(raw.counts))
colnames(raw.counts) <- gsub("-",".",colnames(raw.counts))

# Order the same
raw.counts <-  raw.counts[,annot$sampleID]  
annot$sampleID == colnames(raw.counts)

# Filter lowly expressed genes
libSizes <- colSums(raw.counts$counts)
summary(libSizes)
L <- min(libSizes) # min library size
P <- round(ncol(raw.counts)*0.50) #20% population
dge <- DGEList(counts = raw.counts)
keep <- rowSums(cpm(dge) > 5/L*1e6) > P
dge <- DGEList(dge[keep,,keep.lib.sizes=F]) 

# TMM method -------------------------------------------------------------------------
logCPM_unnorm <- cpm(dge,prior.count = 3, log = T)
dge2 <- dge
dge2 <- calcNormFactors(dge2)
logCPM <- cpm(dge2,prior.count = 3, log = T)
write.csv(logCPM,paste0(outdir,"/Norm.Log.Counts.csv"))

# RLE plots by group
medians <- rowMedians(as.matrix(logCPM_unnorm)) 
rle <- sweep(as.matrix(logCPM_unnorm), MARGIN = 1, medians, `-`)
mediansNorm <- rowMedians(as.matrix(logCPM))
rle.norm <- sweep(as.matrix(logCPM), MARGIN = 1, mediansNorm, `-`)

pdf(paste0(outdir,"/rle_plots_tmm.pdf"))
group <- factor(annot$GenotypeGroup_AgeGroup)
col <- brewer.pal(8,"Set2")[group]
boxplot(rle, outline=FALSE, ylab = "Relative log expression",col=col,main="Raw data",las=2,cex.axis = 0.5)
abline(h = 0)

group <- factor(annot$GenotypeGroup_AgeGroup)
col <- brewer.pal(8,"Set2")[group]
boxplot(rle.norm, outline=FALSE, ylab = "Relative log expression",col=col,main="Normalized data",las=2,cex.axis = 0.5)
abline(h = 0)
dev.off()

# MDS plot
pdf(paste0(outdir,"/MDS_tmm_norm_counts.pdf"))
group <- factor(annot$GenotypeGroup)
col <- brewer.pal(4,"Set1")[group]
plotMDS(logCPM, labels=group, col=col)
group <- factor(annot$GenotypeGroup_AgeGroup)
col <- brewer.pal(8,"Set2")[group]
plotMDS(logCPM, labels=group, col=col)
dev.off()


# RUV method -------------------------------------------------------------------------
k = 10
counts <- dge$counts
annot$sampleID == colnames(counts)

mouse_house_keep <- readRDS("~/Data/gene_sets/scHK_mouse.rds")
mouse_human_map <- read.csv("~/Data/gene_sets/Mouse_human_map.csv")

g <- factor(annot$GenotypeGroup_AgeGroup)
rep_matrix = replicate.matrix(g)
rownames(rep_matrix) = annot$sampleID

mouse_genes_mgi <- mouse_human_map[match(mouse_house_keep,mouse_human_map$hgnc_symbol),"mgi_symbol"]
negative.control.genes <- which(rownames(counts) %in% mouse_genes_mgi)

UQ = betweenLaneNormalization(counts, which = "upper")
ruvsl = RUVIII(Y = t(log2(UQ+0.5)), M = rep_matrix, k=k,ctl = negative.control.genes)
ruvsl <- t(ruvsl)
write.csv(ruvsl,paste0(outdir,"/Norm.Log.Counts.RUV_k",k,".csv"))

# RLE plots by group
medians <- rowMedians(as.matrix(logCPM_unnorm))
rle <- sweep(as.matrix(logCPM_unnorm), MARGIN = 1, medians, `-`)
mediansNorm <- rowMedians(as.matrix(ruvsl))
rle.norm <- sweep(as.matrix(ruvsl), MARGIN = 1, mediansNorm, `-`)

pdf(paste0(outdir,"/rle_plots_ruv_k",k,".pdf"))
group <- factor(annot$GenotypeGroup_AgeGroup)
col <- brewer.pal(8,"Set2")[group]
boxplot(rle, outline=FALSE, ylab = "Relative log expression",col=col,main="Raw data",las=2,cex.axis = 0.5)
abline(h = 0)

group <- factor(annot$GenotypeGroup_AgeGroup)
col <- brewer.pal(8,"Set2")[group]
boxplot(rle.norm, outline=FALSE, ylab = "Relative log expression",col=col,main="Normalized data",las=2,cex.axis = 0.5)
abline(h = 0)
dev.off()

# MDS plot
pdf(paste0(outdir,"/MDS_ruv_norm_counts_k",k,".pdf"))
group <- factor(annot$GenotypeGroup)
col <- brewer.pal(4,"Set1")[group]
plotMDS(ruvsl, labels=group, col=col)
group <- factor(annot$GenotypeGroup_AgeGroup)
col <- brewer.pal(8,"Set2")[group]
plotMDS(ruvsl, labels=group, col=col)
dev.off()
