#!/usr/bin/env Rscript
library(argparse)
parser <- ArgumentParser(description="Run permutation analysis on input cell type (requires data tables in wkDir)")
parser$add_argument("-c", help="cell type (NPC, Neu, or Glia)", type="character", dest="c", required=TRUE) 
args <- parser$parse_args()
source("co-regulation_significance.R")

############################################################
# Read data tables
cell_type <- args$c
cell_intxns <- read.table(paste(cell_type, "-specific_PGC_intrxns.txt",sep=""), sep="\t", header=TRUE, comment.char="")
chromSizes <- read.table("chrom.sizes", sep="\t", header=FALSE)
RPKM_df <- read.table("COSgenesPositionsMatrix_remove_zero_var.txt", sep="\t", header=FALSE, comment.char="")
colnames(RPKM_df)[1:4] <- c("chrom", "start", "end", "gene")
# Separate RPKM table by chrom
RPKM_list <- split(RPKM_df, f=RPKM_df$chrom)

# Remove M chrom
chromSizes <- chromSizes[1:nrow(chromSizes) -1, ]
colnames(chromSizes) <- c("chrom", "bp")
# generate cumulative genome table
chrom <- as.character(chromSizes$chrom)
bp <- cumsum(as.numeric(chromSizes$bp))
cmChromSizes <- data.frame(chrom, bp)
########################################################################################
# Distance equivalent null distribution

permute_test <- function(intxns, cSizes, cmSizes, RPKM_l, RPKM_d, num_samples, perms, tstat, cell) {
  # Run permutation testing
  # Args:
  # intxns = dataframe cell type specific significant interactions 
  # cSizes = dataframe of bp chromosome sizes
  # cmSizes = dataframe of bp cumulative chromosome sizes
  # RPKM_l = list of chromosome dataframes of genes with RPKMs to sample from
  # RPKM_d = full dataframe of genes with RPKMs to sample from
  # sample_size = number of PGC interacting genes
  # perms = number of permutations to run
  # tstat = organization score test-statistic
  # cell = string of cell type
  # Returns:
  # org_scores.txt = organization scores for random gene samples
  # p = p-value for permutation analysis
  
  count = 0
  rndmeans = c()
  for (i in 1:perms)
  {
    rndcos <- sample_null_distribution(intxns, cSizes, cmSizes, RPKM_l, RPKM_d, num_samples)
    rndcos <- rndcos[5:ncol(rndcos)]
    rndcos=t(rndcos)
    n=length(colnames(rndcos))
    crnd <- cor(rndcos[,1:n],use="complete.obs")
    crnd1=abs(crnd)
    rndmean=mean(crnd1)
    rndmeans <- c(rndmeans, rndmean)
    if (rndmean > tstat) count = count+1 
  }
  print(rndcos[1:5,])
  p = count/perms
  print(cell)
  print("count:")
  print(count)
  print("p-value:")
  print(p)
  write(rndmeans, paste(cell, "_org_scores.txt", sep=""), sep="\t")

}
###########################################################################################
if (cell_type == "NPC") {
  permute_test(cell_intxns, chromSizes, cmChromSizes, RPKM_list, RPKM_df, 290, 100, 0.2455, "NPC")
} else if (cell_type == "Neu") {
  permute_test(cell_intxns, chromSizes, cmChromSizes, RPKM_list, RPKM_df, 314, 100, 0.2329, "Neu")
} else if (cell_type == "Glia") {
  permute_test(cell_intxns, chromSizes, cmChromSizes, RPKM_list, RPKM_df, 151, 100, 0.2066, "Glia")
} else {
  print("ERROR: wrong cell type")
  stop()
}
########################################################################################

