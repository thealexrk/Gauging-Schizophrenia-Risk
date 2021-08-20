#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Organization scores
df <- read.table(args, sep="\t", header=FALSE)
colnames(df)[1:4] <- c("chrom", "start", "end", "gene")
df <- df[5:ncol(df)]
row.names(df) <- 1:length(df[,1])
df=t(df)
n=length(colnames(df))
crnd <- cor(df[,1:n],use="complete.obs")
crnd1=abs(crnd)
rndmean=mean(crnd1)
print(rndmean)