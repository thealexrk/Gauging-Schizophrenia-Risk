z.test = function(x,mu,popvar){
  z.score <- (mean(x)-mu)/(popvar/sqrt(length(x)))
  print(z.score)
  one.tail.p <- pnorm(abs(z.score),lower.tail = FALSE)
  return(one.tail.p)
}

df <- read.table("gene_table", sep="\t", header=TRUE, quote="", comment.char="")
Astro_FPKM <- read.table("553_v._primary_astro_RNA_seq.txt", sep="\t", header=TRUE, quote="", comment.char="")

##################################################################################################
# NPC specific
NPC_specific <- read.table("NPC_overlapping_genes_cell_type_specific.txt", sep="\t", header=FALSE, quote="", comment.char="")
colnames(NPC_specific) <- colnames(df)

# Astrocyte data
Astro_FPKM <- Astro_FPKM[c(4,10)]
colnames(Astro_FPKM) <- c("Ensembl.Gene", "Astro_FPKM_553")
Astro_FPKM <- Astro_FPKM[!duplicated(Astro_FPKM$Ensembl.Gene),]

NPC_specific <- merge(NPC_specific, Astro_FPKM, by="Ensembl.Gene")

png('FPKM_NPC_specific_boxplots.png', height=1800, width=1600, res=300)
boxplot(log10(df$X2607.1.AN.2 + 1), log10(NPC_specific$X2607.1.AN.2 + 1), 
        log10(df$X2607.1.4.ngn2.2 + 1), log10(NPC_specific$X2607.1.4.ngn2.2 + 1),
        log10(Astro_FPKM$Astro_FPKM_553 + 1), log10(NPC_specific$Astro_FPKM_553 + 1),
        names=c("All", "Anchor", "All", "Anchor", "All", "Anchor"), 
        main= "NPC specific loops",
        ylab=expression('log'[10]*'(FPKM + 1)'), xlab='Gene sets',
        col=c(rep("#d95f02", 2), rep("#1b9e77",2), rep("#7570b3",2)),
        ylim=c(0,3), outline = FALSE, cex.axis=0.8)
legend("topleft", c("NPC", "Neu", "Astro"), 
       fill=c("#d95f02", "#1b9e77", "#7570b3"), bty="n")

p1 = z.test(log10(NPC_specific$X2607.1.AN.2 + 1), mean(log10(df$X2607.1.AN.2 + 1)), 
           sd(log10(df$X2607.1.AN.2 + 1)))
p2 = z.test(log10(NPC_specific$X2607.1.4.ngn2.2 + 1), mean(log10(df$X2607.1.4.ngn2.2 + 1)), 
            sd(log10(df$X2607.1.4.ngn2.2 + 1)))
p3 = z.test(log10(NPC_specific$Astro_FPKM_553 + 1), mean(log10(Astro_FPKM$Astro_FPKM_553 + 1)), 
            sd(log10(Astro_FPKM$Astro_FPKM_553 + 1)))

text(c(1, 3, 5) + .25, rep(1.75, 3), cex=0.75,
     c(paste("p =", sprintf("%.0e",p1, sep="")),
       paste("p = ", sprintf("%.0e",p2, sep="")),
       paste("p = ", sprintf("%.0e",p3, sep=""))))

dev.off()

print("*****NPC specific pvals*****")
print("NPC")
print(p1)
print("Neu")
print(p2)
print("Astro")
print(p3)

##################################################################################################
# Neu specific

Neu_specific <- read.table("Neu_overlapping_genes_cell_type_specific.txt", sep="\t", header=FALSE, quote="", comment.char="")
colnames(Neu_specific) <- colnames(df)

Neu_specific <- merge(Neu_specific, Astro_FPKM, by="Ensembl.Gene")

png('FPKM_Neu_specific_boxplots.png', height=1800, width=1600, res=300)
boxplot(log10(df$X2607.1.AN.2 + 1), log10(Neu_specific$X2607.1.AN.2 + 1), 
        log10(df$X2607.1.4.ngn2.2 + 1), log10(Neu_specific$X2607.1.4.ngn2.2 + 1),
        log10(Astro_FPKM$Astro_FPKM_553 + 1), log10(Neu_specific$Astro_FPKM_553 + 1),
        names=c("All", "Anchor", "All", "Anchor", "All", "Anchor"), 
        main= "Neu specific loops",
        ylab=expression('log'[10]*'(FPKM + 1)'), xlab='Gene sets',
        col=c(rep("#d95f02", 2), rep("#1b9e77",2), rep("#7570b3",2)),
        ylim=c(0,3), outline = FALSE, cex.axis=0.8)
legend("topleft", c("NPC", "Neu", "Astro"), 
       fill=c("#d95f02", "#1b9e77", "#7570b3"), bty="n")


p1 = z.test(log10(Neu_specific$X2607.1.AN.2 + 1), mean(log10(df$X2607.1.AN.2 + 1)), 
            sd(log10(df$X2607.1.AN.2 + 1)))
p2 = z.test(log10(Neu_specific$X2607.1.4.ngn2.2 + 1), mean(log10(df$X2607.1.4.ngn2.2 + 1)), 
            sd(log10(df$X2607.1.4.ngn2.2 + 1)))
p3 = z.test(log10(Neu_specific$Astro_FPKM_553 + 1), mean(log10(Astro_FPKM$Astro_FPKM_553 + 1)), 
                  sd(log10(Astro_FPKM$Astro_FPKM_553 + 1)))

text(c(1, 3, 5) + .25, rep(1.75, 3), cex=0.75,
     c(paste("p =", sprintf("%.0e",p1, sep="")),
       paste("p = ", sprintf("%.0e",p2, sep="")),
       paste("p = ", sprintf("%.0e",p3, sep=""))))


dev.off()

print("*****Neu specific pvals*****")
print("NPC")
print(p1)
print("Neu")
print(p2)
print("Astro")
print(p3)


##################################################################################################
# Astro specific

Astro_specific <- read.table("Astro_overlapping_genes_cell_type_specific.txt", sep="\t", header=FALSE, quote="", comment.char="")
colnames(Astro_specific) <- colnames(df)

Astro_specific <- merge(Astro_specific, Astro_FPKM, by="Ensembl.Gene")

png('FPKM_Astro_specific_boxplots.png', height=1800, width=1600, res=300)
boxplot(log10(df$X2607.1.AN.2 + 1), log10(Astro_specific$X2607.1.AN.2 + 1), 
        log10(df$X2607.1.4.ngn2.2 + 1), log10(Astro_specific$X2607.1.4.ngn2.2 + 1),
        log10(Astro_FPKM$Astro_FPKM_553 + 1), log10(Astro_specific$Astro_FPKM_553 + 1),
        names=c("All", "Anchor", "All", "Anchor", "All", "Anchor"), 
        main= "Astro specific loops",
        ylab=expression('log'[10]*'(FPKM + 1)'), xlab='Gene sets',
        col=c(rep("#d95f02", 2), rep("#1b9e77",2), rep("#7570b3",2)),
        ylim=c(0,3), outline = FALSE, cex.axis=0.8)
legend("topleft", c("NPC", "Neu", "Astro"), 
       fill=c("#d95f02", "#1b9e77", "#7570b3"), bty="n")


p1 = z.test(log10(Astro_specific$X2607.1.AN.2 + 1), mean(log10(df$X2607.1.AN.2 + 1)), 
            sd(log10(df$X2607.1.AN.2 + 1)))
p2 = z.test(log10(Astro_specific$X2607.1.4.ngn2.2 + 1), mean(log10(df$X2607.1.4.ngn2.2 + 1)), 
            sd(log10(df$X2607.1.4.ngn2.2 + 1)))
p3 = z.test(log10(Astro_specific$Astro_FPKM_553 + 1), mean(log10(Astro_FPKM$Astro_FPKM_553 + 1)), 
            sd(log10(Astro_FPKM$Astro_FPKM_553 + 1)))


text(c(1, 3, 5) + .25, rep(1.75, 3), cex=0.75,
     c(paste("p =", sprintf("%.0e",p1, sep="")),
       paste("p = ", sprintf("%.0e",p2, sep="")),
       paste("p = ", sprintf("%.0e",p3, sep=""))))

dev.off()

print("*****Astro specific pvals*****")
print("NPC")
print(p1)
print("Neu")
print(p2)
print("Astro")
print(p3)


##################################################################################################
# Neu specific filtered > 50 reads at loop pixel

Neu_specific <- read.table("Neu_overlapping_genes_cell_type_specific_filter.txt", sep="\t", header=FALSE, quote="", comment.char="")
colnames(Neu_specific) <- colnames(df)

Neu_specific <- merge(Neu_specific, Astro_FPKM, by="Ensembl.Gene")

png('FPKM_Neu_specific_filter_boxplots.png', height=1800, width=1600, res=300)
boxplot(log10(df$X2607.1.AN.2 + 1), log10(Neu_specific$X2607.1.AN.2 + 1), 
        log10(df$X2607.1.4.ngn2.2 + 1), log10(Neu_specific$X2607.1.4.ngn2.2 + 1),
        log10(Astro_FPKM$Astro_FPKM_553 + 1), log10(Neu_specific$Astro_FPKM_553 + 1),
        names=c("All", "Anchor", "All", "Anchor", "All", "Anchor"), 
        main= "Neu specific loops",
        ylab=expression('log'[10]*'(FPKM + 1)'), xlab='Gene sets',
        col=c(rep("#d95f02", 2), rep("#1b9e77",2), rep("#7570b3",2)),
        ylim=c(0,3), outline = FALSE, cex.axis=0.8)
legend("topleft", c("NPC", "Neu", "Astro"), 
       fill=c("#d95f02", "#1b9e77", "#7570b3"), bty="n")


p1 = z.test(log10(Neu_specific$X2607.1.AN.2 + 1), mean(log10(df$X2607.1.AN.2 + 1)), 
            sd(log10(df$X2607.1.AN.2 + 1)))
p2 = z.test(log10(Neu_specific$X2607.1.4.ngn2.2 + 1), mean(log10(df$X2607.1.4.ngn2.2 + 1)), 
            sd(log10(df$X2607.1.4.ngn2.2 + 1)))
p3 = z.test(log10(Neu_specific$Astro_FPKM_553 + 1), mean(log10(Astro_FPKM$Astro_FPKM_553 + 1)), 
            sd(log10(Astro_FPKM$Astro_FPKM_553 + 1)))


text(c(1, 3, 5) + .25, rep(1.75, 3), cex=0.75,
     c(paste("p =", sprintf("%.0e",p1, sep="")),
       paste("p = ", sprintf("%.0e",p2, sep="")),
       paste("p = ", sprintf("%.0e",p3, sep=""))))

dev.off()

print("*****Neu specific pvals*****")
print("NPC")
print(p1)
print("Neu")
print(p2)
print("Astro")
print(p3)


