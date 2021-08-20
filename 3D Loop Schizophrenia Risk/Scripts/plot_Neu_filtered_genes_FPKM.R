gene <- read.table("gene_table", sep="\t", header=TRUE, quote="", comment.char="")

df <- read.table("Neu_overlapping_genes_cell_type_specific_filter.txt", sep="\t", header=FALSE, comment.char="")
colnames(df) <- colnames(gene)


names = c("cut-like homeobox 2",
          "family with sequence similarity 120B",
          "transmembrane protease, serine 13",
          "ketohexokinase",
          "neural cell adhesion molecule 1",
          "Rho GTPase activating protein 35",
          "neurofascin",
          "beta-transducin repeat containing E3 ubiquitin protein ligase",
          "long intergenic non-protein coding RNA 242",
          "long intergenic non-protein coding RNA 574",
          "","","")


# Astrocyte data
Astro_FPKM <- read.table("553_v._primary_astro_RNA_seq.txt", sep="\t", header=TRUE, quote="", comment.char="")
Astro_FPKM <- Astro_FPKM[c(4,10)]
colnames(Astro_FPKM) <- c("Ensembl.Gene", "Astro_FPKM_553")
Astro_FPKM <- Astro_FPKM[!duplicated(Astro_FPKM$Ensembl.Gene),]
df <- merge(df, Astro_FPKM, by="Ensembl.Gene")



png("plot_Neu_filtered_genes_FPKM.png",height=1500, width=2500, res=300)
par(las=1,mar=c(5,25,4,2))
barplot(log10(df$X2607.1.4.ngn2.2 + 1), horiz=TRUE, names.arg = names,
        col="#1b9e77", xlab=expression('log'[10]*'(FPKM + 1)'), xlim=c(0,1.65))
dev.off()

png("plot_NPC_filtered_genes_FPKM.png",height=1500, width=2500, res=300)
par(las=1,mar=c(5,25,4,2))
barplot(log10(df$X2607.1.AN.2 + 1), horiz=TRUE, names.arg = names,
        col="#d95f02", xlab=expression('log'[10]*'(FPKM + 1)'), xlim=c(0,1.65))
dev.off()

png("plot_Astro_filtered_genes_FPKM.png",height=1500, width=2500, res=300)
par(las=1,mar=c(5,25,4,2))
barplot(log10(df$Astro_FPKM_553 + 1), horiz=TRUE, names.arg = names,
        col="#7570b3", xlab=expression('log'[10]*'(FPKM + 1)'), xlim=c(0,1.65))
dev.off()
