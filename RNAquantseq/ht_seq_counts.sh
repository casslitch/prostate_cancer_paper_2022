bam_path="/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/"
gene_gtf="/data/reference/indexes/mouse/ensembl_GRCm38.73/bowtie2/Tophat2_transcriptome_data/Mus_musculus.GRCm38.73.gtf"
outdir="/researchers/simon.keam/RNA-QuantSeq/220120_NB501056_0747_AHFFW7BGXL/ht_seq_counts/"
logdir="${outdir}/logs/"

for bam in $bam_path/*/Bam/*.bam; do
  sampID=$(basename $bam _sorted.bam)
  echo $sampID
  outlog="${logdir}/${sampID}.out"
  sbatch --output=${outlog} ht_seq_count_RNA_quantseq.sbatch $bam $gene_gtf $outdir/${sampID}.txt
done
