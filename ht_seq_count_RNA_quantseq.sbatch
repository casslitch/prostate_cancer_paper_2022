#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --job-name=ht_seq_count_%j
#SBATCH --time=00-2:00:00
#SBATCH --mem=16GB
#SBATCH --partition=prod_short

module purge
module load htseq/0.9.1

bam=$1
gene_annot_gtf=$2
outfile=$3

# For 3' RNA quant seq

htseq-count \
-f bam \
-m intersection-nonempty \
-i gene_name \
-r pos \
-s yes \
$bam \
$gene_annot_gtf > $outfile