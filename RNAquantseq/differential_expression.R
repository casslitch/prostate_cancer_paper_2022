# Author: Cassie Litchfield
# Description: DE analyses (limma-trend) using normalized log cpm

library(dplyr)
library(limma)
library(Mus.musculus)
library(openxlsx)

indir <- "/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/analyses/"
outdir <- "/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/analyses/"

# -------------------------------------------------------------------------
doDEanalyses <- function(group1,group2,fit,design){
  cont <- makeContrasts(contrasts = paste0(group1,"-",group2),levels = design)
  fit2 <- contrasts.fit(fit,cont)
  fit2 <- eBayes(fit2, trend = T, robust = T)
  
  top.t <- topTable(fit2, sort.by = 'p', number = Inf, p.value = 1,adjust.method = 'BH')
  top.t$SYMBOL <-  rownames(top.t)
  top.t$IS_SIGNIFICANT = top.t$adj.P.Val < 0.05
  return(top.t)
}

# -------------------------------------------------------------------------
# Load 
normCounts <-  read.csv(paste0(outdir,"/","Norm.Log.Counts.RUV_k10.csv"), row.names = 1)
# genes <- select(Mus.musculus, keys=rownames(normCounts), columns=c("TXCHROM","ENTREZID"),
#                 keytype="SYMBOL")
# genes <- genes[match(rownames(normCounts),genes$SYMBOL),] #keep first occurrence

# Read sample annotation
annot <- read.xlsx(paste0(indir,"metadata_v2.xlsx"), cols = 1:6)
annot$sampleID <- paste0("SampleID.",annot$SampleNo)
annot$sampleID <- gsub(" \\(",".",annot$sampleID)
annot$sampleID <- gsub("\\)","",annot$sampleID)

# Order the same
normCounts <-  normCounts[,annot$sampleID]  
annot$sampleID == colnames(normCounts)

# -------------------------------------------------------------------------
# DE analyses 
design <- model.matrix(~0 + factor(annot$GenotypeGroup_AgeGroup))
colnames(design) <- gsub('factor(annot$GenotypeGroup_AgeGroup)',"",colnames(design), fixed = T)

fit <- lmFit(normCounts, design = design)

group1 <- "TP53null_Early"
group2 <- "Control_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R172_Early"
group2 <- "Control_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R245_Early"
group2 <- "Control_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)


group1 <- "TP53null_Late"
group2 <- "Control_Late"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)


group1 <- "TP53_R172_Late"
group2 <- "Control_Late"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R245_Late"
group2 <- "Control_Late"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)


group1 <- "TP53_R245_Late"
group2 <- "TP53_R245_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)


group1 <- "TP53_R172_Late"
group2 <- "TP53_R172_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53null_Late"
group2 <- "TP53null_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "Control_Late"
group2 <- "Control_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R172_Early"
group2 <- "TP53_R245_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R172_Early"
group2 <- "TP53null_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R245_Early"
group2 <- "TP53null_Early"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R172_Late"
group2 <- "TP53_R245_Late"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R172_Late"
group2 <- "TP53null_Late"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)

group1 <- "TP53_R245_Late"
group2 <- "TP53null_Late"
top.t <- doDEanalyses(group1,group2,fit,design)
write.csv(top.t, paste0(outdir,"/",group1,"vs",group2,".csv"), row.names = F)
