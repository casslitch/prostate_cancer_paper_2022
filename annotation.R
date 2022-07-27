indir <- "/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/analyses/"
outdir <- "/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/analyses/"
fastqpath <- "~/ocavtio_data/" # copied the fastqs here

# Read sample annotation
annot <- read.xlsx(paste0(indir,"metadata_v2.xlsx"), cols = 1:6)
annot$sampleID <- paste0("SampleID.",annot$SampleNo)
annot$sampleID <- gsub(" \\(",".",annot$sampleID)
annot$sampleID <- gsub("\\)","",annot$sampleID)

# get the fastq file name
get_file <- function(x){
  tmp <- list.files(fastqpath)
  return(tmp[grepl(x,tmp)])
}

annot$fastq <- unlist(sapply(gsub("\\.","-",annot$sampleID),get_file))

write.csv(annot,paste0(outdir,"annot.csv"),row.names=F)
