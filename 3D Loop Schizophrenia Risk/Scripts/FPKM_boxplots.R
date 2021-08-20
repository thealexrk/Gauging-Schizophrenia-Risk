df <- read.table("gene_table", sep="\t", header=TRUE, quote="", comment.char="")

Astro_FPKM <- read.table("553_v._primary_astro_RNA_seq.txt", sep="\t", header=TRUE, quote="", comment.char="")

NPC_df <- read.table("NPC_overlapping_genes.txt", sep="\t", header=FALSE, quote="", comment.char="")
NPC_tss_df <- read.table("NPC_overlapping_tss.txt", sep="\t", header=FALSE, quote="", comment.char="")
NPC_full_df <- read.table("NPC_overlapping_genes_full_loop.txt", sep="\t", header=FALSE, quote="", comment.char="")

Astro_df <- read.table("Astro_overlapping_genes.txt", sep="\t", header=FALSE, quote="", comment.char="")
Astro_tss_df <- read.table("Astro_overlapping_tss.txt", sep="\t", header=FALSE, quote="", comment.char="")
Astro_full_df <- read.table("Astro_overlapping_genes_full_loop.txt", sep="\t", header=FALSE, quote="", comment.char="")

Neu_df <- read.table("Neu_overlapping_genes.txt", sep="\t", header=FALSE, quote="", comment.char="")
Neu_tss_df <- read.table("Neu_overlapping_tss.txt", sep="\t", header=FALSE, quote="", comment.char="")
Neu_full_df <- read.table("Neu_overlapping_genes_full_loop.txt", sep="\t", header=FALSE, quote="", comment.char="")

colnames(NPC_df) <- colnames(df)
colnames(NPC_tss_df) <- colnames(df)
colnames(NPC_full_df) <- colnames(df)
colnames(Astro_df) <- colnames(df)
colnames(Astro_tss_df) <- colnames(df)
colnames(Astro_full_df) <- colnames(df)
colnames(Neu_df) <- colnames(df)
colnames(Neu_tss_df) <- colnames(df)
colnames(Neu_full_df) <- colnames(df)

# Astrocyte data
Astro_FPKM <- Astro_FPKM[c(4,10)]
colnames(Astro_FPKM)[1] <- "Ensembl.Gene"
Astro_FPKM <- Astro_FPKM[!duplicated(Astro_FPKM$Ensembl.Gene),]

Astro_full_df <- merge(Astro_full_df, Astro_FPKM, by="Ensembl.Gene")
Astro_df <- merge(Astro_df, Astro_FPKM, by="Ensembl.Gene")
Astro_tss_df <- merge(Astro_tss_df, Astro_FPKM, by="Ensembl.Gene")

FPKM_list <- list(log10(df$X2607.1.AN.2 + 1), log10(NPC_full_df$X2607.1.AN.2 + 1),
                 log10(NPC_df$X2607.1.AN.2 + 1), log10(NPC_tss_df$X2607.1.AN.2 + 1),
                 log10(df$X2607.1.4.ngn2.2 + 1), log10(Neu_full_df$X2607.1.4.ngn2.2 + 1),
                 log10(Neu_df$X2607.1.4.ngn2.2 + 1), log10(Neu_tss_df$X2607.1.4.ngn2.2 + 1),
                 log10(Astro_FPKM$FPKM + 1), log10(Astro_full_df$FPKM + 1), 
                 log10(Astro_df$FPKM + 1), log10(Astro_tss_df$FPKM + 1))


z.test = function(x,mu,popvar){
  z.score <- (mean(x)-mu)/(popvar/sqrt(length(x)))
  print(z.score)
  one.tail.p <- pnorm(abs(z.score),lower.tail = FALSE)
  return(one.tail.p)
}

# pop_expr = log10(df$X2607.1.AN.2 + 1)
# loop_expr = log10(NPC_df$X2607.1.AN.2 + 1)
# 
# p <- z.test(loop_expr, mean(pop_expr), sd(pop_expr))
# print(p)

print("*****NPC p-vals*****")
print("full:")
print(z.test(FPKM_list[[2]], mean(FPKM_list[[1]]), sd(FPKM_list[[1]])))
print("anchor:")
print(z.test(FPKM_list[[3]], mean(FPKM_list[[1]]), sd(FPKM_list[[1]])))
print("tss:")
print(z.test(FPKM_list[[4]], mean(FPKM_list[[1]]), sd(FPKM_list[[1]])))
print("*****Neu p-vals*****")
print("full:")
print(z.test(FPKM_list[[6]], mean(FPKM_list[[5]]), sd(FPKM_list[[5]])))
print("anchor:")
print(z.test(FPKM_list[[7]], mean(FPKM_list[[5]]), sd(FPKM_list[[5]])))
print("tss:")
print(z.test(FPKM_list[[8]], mean(FPKM_list[[5]]), sd(FPKM_list[[5]])))
print("****Astro p-vals****")
print("full:")
print(z.test(FPKM_list[[10]], mean(FPKM_list[[9]]), sd(FPKM_list[[9]])))
print("anchor:")
print(z.test(FPKM_list[[11]], mean(FPKM_list[[9]]), sd(FPKM_list[[9]])))
print("tss:")
print(z.test(FPKM_list[[12]], mean(FPKM_list[[9]]), sd(FPKM_list[[9]])))


png('FPKM_boxplots.png', height=1800, width=2600, res=300)
boxplot(FPKM_list, names=rep(c("All", "Loop", "Anchor", "TSS"),3), 
        ylab=expression('log'[10]*'(FPKM + 1)'), xlab='Gene sets',
        col=c(rep("#d95f02", 4), rep("#1b9e77",4), rep("#7570b3",4)),
        ylim=c(0,3), outline = FALSE, cex.axis=0.8)
legend("topright", c("NPC", "Neu", "Astro"), 
       fill=c("#d95f02", "#1b9e77", "#7570b3"), bty="n")
dev.off()
