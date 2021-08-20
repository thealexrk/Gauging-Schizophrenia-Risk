
lo <- read.table("overlap_loop_compartments.txt", header=TRUE, sep="\t", row.names=1)
# Remove total and U
rlo <- lo[,1:3]
tlo <- as.table(t(rlo))
colnames(tlo) <- c("Glia", "NPC", "Neuron")
rownames(tlo) <- c("AA", "BB", "AB")


pdf("downsampled_compartment_loop_overlaps.pdf", width=4, height=7)
par(mar=c(5, 5, 4, 2) + 0.1)
barplot(tlo, col=c("blue", "red", "violet"), legend.text=rownames(tlo), 
        ylab= "Number of loops", cex.lab=1.3, cex.names=1.3, cex.axis = 1.15, 
        args.legend=list(x=4.3, y=13000, bty="n", cex=1.3), ylim=c(0,12000))
dev.off()


# Percentages
print(paste("Astro_AA: ", round(lo["Astro", "AA"]/lo["Astro", "total"], 2)*100,
            "%", sep=""))
print(paste("Astro_BB: ", round(lo["Astro", "BB"]/lo["Astro", "total"], 2)*100,
            "%", sep=""))
print(paste("Astro_AB: ", round(lo["Astro", "AB"]/lo["Astro", "total"], 2)*100,
            "%", sep=""))
print(paste("NPC_AA: ", round(lo["NPC", "AA"]/lo["NPC", "total"], 2)*100,
            "%", sep=""))
print(paste("NPC_BB: ", round(lo["NPC", "BB"]/lo["NPC", "total"], 2)*100,
            "%", sep=""))
print(paste("NPC_AB: ", round(lo["NPC", "AB"]/lo["NPC", "total"], 2)*100,
            "%", sep=""))
print(paste("Neu_AA: ", round(lo["Neu", "AA"]/lo["Neu", "total"], 2)*100,
            "%", sep=""))
print(paste("Neu_BB: ", round(lo["Neu", "BB"]/lo["Neu", "total"], 2)*100,
            "%", sep=""))
print(paste("Neu_AB: ", round(lo["Neu", "AB"]/lo["Neu", "total"], 2)*100,
            "%", sep=""))

# Compartment size correct
Neu <- read.table("Neu.553x2_2607.11654.500kb.pc1_removed_chr.bedGraph", header=FALSE, sep="\t", comment.char="")
colnames(Neu) <- c("chrom", "start", "stop", "eigen")
Neu_A <- Neu[Neu$eigen > 0, ]
Neu_A_size <- sum(Neu_A$stop - Neu_A$start, na.rm=TRUE)
Neu_B <- Neu[Neu$eigen <= 0, ]
Neu_B_size <- sum(Neu_B$stop - Neu_B$start, na.rm=TRUE)

Astro <- read.table("Astro.553x2_2607.11654.500kb.pc1_removed_chr.bedGraph", header=FALSE, sep="\t", comment.char="")
colnames(Astro) <- c("chrom", "start", "stop", "eigen")
Astro_A <- Astro[Astro$eigen > 0, ]
Astro_A_size <- sum(Astro_A$stop - Astro_A$start, na.rm=TRUE)
Astro_B <- Astro[Astro$eigen <= 0, ]
Astro_B_size <- sum(Astro_B$stop - Astro_B$start, na.rm=TRUE)

NPC <- read.table("NPC.553x2_2607.11654.500kb.pc1_removed_chr.bedGraph", header=FALSE, sep="\t", comment.char="")
colnames(NPC) <- c("chrom", "start", "stop", "eigen")
NPC_A <- NPC[NPC$eigen > 0, ]
NPC_A_size <- sum(NPC_A$stop - NPC_A$start, na.rm=TRUE)
NPC_B <- NPC[NPC$eigen <= 0, ]
NPC_B_size <- sum(NPC_B$stop - NPC_B$start, na.rm=TRUE)

# Correct for size of compartments
Neu_A_correct <- (lo["Neu", "AA"]*1000000)/(Neu_A_size)
Neu_B_correct <- (lo["Neu", "BB"]*1000000)/(Neu_B_size)
NPC_A_correct <- (lo["NPC", "AA"]*1000000)/(NPC_A_size)
NPC_B_correct <- (lo["NPC", "BB"]*1000000)/(NPC_B_size)
Astro_A_correct <- (lo["Astro", "AA"]*1000000)/(Astro_A_size)
Astro_B_correct <- (lo["Astro", "BB"]*1000000)/(Astro_B_size)

Neu_correct <- c(Neu_A_correct, Neu_B_correct)
NPC_correct <- c(NPC_A_correct, NPC_B_correct)
Glia_correct <- c(Astro_A_correct, Astro_B_correct)

correct_mat <- as.matrix(data.frame(Glia_correct, NPC_correct, Neu_correct, row.names=c("A", "B")))
colnames(correct_mat) <- c("Glia", "NPC", "Neuron")

pdf("downsampled_compartment_loop_overlaps_size_correct.pdf", width=5, height=5)
par(mar=c(5.1, 5.1, 4.1, 2.1))
correct_tab <- as.table(correct_mat)
barplot(correct_tab, beside = TRUE, col=c("blue", "red"), 
        ylab= "Number of loops per 1 Mb", ylim= c(0,7))
legend("topright", rownames(correct_tab), fill=c("blue", "red"), bty="n")
dev.off()

