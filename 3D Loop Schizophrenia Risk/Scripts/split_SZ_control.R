#!/usr/bin/env Rscript
library(pheatmap)
library(RColorBrewer)
#args = commandArgs(trailingOnly = TRUE)

args= "COS_CLOZUK_PGC_Neu_header.txt"

get_org_score <- function(df){
  # Return org_score and cormatrix
  df=t(df)
  n=length(colnames(df))
  cor_m <- cor(df[,1:n],use="complete.obs")
  abs_cor_m=abs(cor_m)
  org_score=mean(abs_cor_m)
  return (list(org_score, cor_m))
} 

plot_cor_heatmap <- function(cormat, f, cluster) {
  # Plot correlation heatmap
  pdf(f, width=7.25, height=7, onefile=FALSE)
    myBreaks <- c(seq(-1, 0, length.out=ceiling(100/2) + 1),seq(1/100, 1, length.out=floor(100/2)))
    p1 <- pheatmap(cormat, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA,
                   breaks=myBreaks, show_rownames=FALSE, show_colnames=FALSE,
                   cluster_rows = cluster, cluster_cols=cluster)
  dev.off()
  return(p1)
}

# Organization scores
df <- read.table(args, sep="\t", header=TRUE)
df <- df[6:ncol(df)]
row.names(df) <- 1:length(df[,1])
genes <- as.character(df[,1])
r <- get_org_score(df)
os <- r[[1]]
cormat <- r[[2]]
print(os)

p1 <- plot_cor_heatmap(cormat, "COS_CLOZUK_PGC_Neu_heatmap.pdf", TRUE) 

# Reordered matrix from clustering
roworder <- p1$tree_row$order
ordered_genes <- data.frame(genes[roworder])
colnames(ordered_genes) <- "genes"
ctree <- cutree(p1$tree_row, k=4)
clust1_idx <- as.numeric(which(ctree == 4))
clust1_ordered_idx <- is.element(roworder, clust1_idx)
clust1_genes <- ordered_genes[clust1_ordered_idx,]

clust2_idx <- as.numeric(which(ctree == 2))
clust2_ordered_idx <- is.element(roworder, clust2_idx)
clust2_genes <- ordered_genes[clust2_ordered_idx,]

df <- read.table(args, sep="\t", header=TRUE)
df <- df[6:ncol(df)]
ordered_df <- df[roworder,]
clust_df <- ordered_df[clust1_ordered_idx | clust2_ordered_idx,]

r <- get_org_score(clust_df)
os <- r[[1]]
cormat <- r[[2]]
print(os)

p1 <- plot_cor_heatmap(cormat, "COS_CLOZUK_PGC_Neu_clusters_heatmap.pdf", FALSE) 

# SCZ vs Control
meta_df  <- read.table("metadata.csv", sep=",", header=TRUE)
#SCZ
SZ_df <- meta_df[meta_df$Dx == "SZ",]
SZ_IDs <- paste('X', SZ_df$Sample.Name, sep="")
SZ_IDs <- gsub("-", ".", SZ_IDs)
SZ_clust_df <- clust_df[SZ_IDs]
r <- get_org_score(SZ_clust_df)
os <- r[[1]]
cormat <- r[[2]]
print(os)

p1 <- plot_cor_heatmap(cormat, "COS_CLOZUK_PGC_Neu_clusters_SCZ_heatmap.pdf", FALSE) 


#Control
CT_df <- meta_df[meta_df$Dx == "CT",]
CT_IDs <- paste('X', CT_df$Sample.Name, sep="")
CT_IDs <- gsub("-", ".", CT_IDs)
CT_IDs <- gsub("R", ".2", CT_IDs)
CT_clust_df <- clust_df[CT_IDs]
r <- get_org_score(CT_clust_df)
os <- r[[1]]
cormat <- r[[2]]
print(os)
p1 <- plot_cor_heatmap(cormat, "COS_CLOZUK_PGC_Neu_clusters_CT_heatmap.pdf", FALSE) 


