#!/usr/bin/env Rscript
library(argparse)
library(pheatmap)
library(RColorBrewer)
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
  # print(rndcos[1:5,])
  # p = count/perms
  # print(cell)
  # print("count:")
  # print(count)
  # print("p-value:")
  # print(p)
  # #write(rndmeans, paste(cell, "_org_scores.txt", sep=""), sep="\t")
  print("Observation Score: ")
  print(rndmean)
  pdf("NPC_random_heatmap.pdf", width=7.25, height=7)
    myBreaks <- c(seq(-1, 0, length.out=ceiling(100/2) + 1),seq(1/100, 1, length.out=floor(100/2)))
    pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA,
             breaks=myBreaks, show_rownames=FALSE, show_colnames=FALSE)
  dev.off()
}
###########################################################################################
if (cell_type == "NPC") {
  permute_test(cell_intxns, chromSizes, cmChromSizes, RPKM_list, RPKM_df, 290, 1, 0.2455, "NPC")
}
########################################################################################

