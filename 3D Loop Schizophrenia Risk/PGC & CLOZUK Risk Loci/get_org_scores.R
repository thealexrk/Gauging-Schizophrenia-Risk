#!/usr/bin/env Rscript
library(pheatmap)
library(RColorBrewer)
args = commandArgs(trailingOnly = TRUE)

# Organization scores
df <- read.table(args, sep="\t", header=FALSE)
colnames(df)[1:4] <- c("chrom", "start", "end", "gene")
genes <- as.character(df[,4])
df <- df[5:ncol(df)]
row.names(df) <- 1:length(df[,1])
df=t(df)
n=length(colnames(df))
crnd <- cor(df[,1:n],use="complete.obs")
crnd1=abs(crnd)
rndmean=mean(crnd1)
print(rndmean)

heatmap_file <- paste(substr(args, 1, nchar(args) - 4), "_heatmap.pdf", sep="")
pdf(heatmap_file, width=7.25, height=7, onefile=FALSE)
myBreaks <- c(seq(-1, 0, length.out=ceiling(100/2) + 1),seq(1/100, 1, length.out=floor(100/2)))
pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA,
         breaks=myBreaks, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

# Write correlation matrix to file
outdf <- rbind(genes, crnd)
rownames(outdf) <- c("", genes)

 write.table(outdf, paste(substr(args, 1, nchar(args) - 4), "_cor_matrix.txt", sep=""), 
             quote=FALSE, sep="\t", row.names = TRUE, col.names=FALSE)
 
 