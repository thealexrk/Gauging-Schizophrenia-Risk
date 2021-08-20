Neu_10kb <- read.table("PM_NeuN_plus/354Neu_downsample_arrowhead_blocks_10kb", sep="\t", header=TRUE)
Neu_25kb <- read.table("PM_NeuN_plus/354Neu_downsample_arrowhead_blocks_25kb", sep="\t", header=TRUE)
Neu_50kb <- read.table("PM_NeuN_plus/354Neu_downsample_arrowhead_blocks_50kb", sep="\t", header=TRUE)
Neu_100kb <- read.table("PM_NeuN_plus/354Neu_downsample_arrowhead_blocks_100kb", sep="\t", header=TRUE)

Glia_10kb <- read.table("PM_NeuN_minus/354Glia_downsample_arrowhead_blocks_10kb", sep="\t", header=TRUE)
Glia_25kb <- read.table("PM_NeuN_minus/354Glia_downsample_arrowhead_blocks_25kb", sep="\t", header=TRUE)
Glia_50kb <- read.table("PM_NeuN_minus/354Glia_downsample_arrowhead_blocks_50kb", sep="\t", header=TRUE)
Glia_100kb <- read.table("PM_NeuN_minus/354Glia_downsample_arrowhead_blocks_100kb", sep="\t", header=TRUE)

count_table <- matrix(c(nrow(Neu_10kb), nrow(Neu_25kb), nrow(Neu_50kb), nrow(Neu_100kb),
                        nrow(Glia_10kb), nrow(Glia_25kb), nrow(Glia_50kb), nrow(Glia_100kb)), 
                        ncol=4, byrow=TRUE)
colnames(count_table) <- c("10kb", "25kb", "50kb", "100kb")
rownames(count_table) <- c("PM NeuN+", "PM NeuN-")

count_table <- as.table(count_table)

png("PM_TAD_counts.png", width=2500, height=2500, res=300)
par(lwd=2, mar=c(5.1, 5.1, 4.1, 2.1))
barplot(count_table, beside=TRUE, col = c("#1b9e77", "#7570b3"),
        ylab = "Number of TADs", cex.names=1.5, ylim=c(0,4200),cex.axis=1.25, cex.lab=1.5)
legend("topright", rownames(count_table), fill= c("#1b9e77", "#7570b3"), cex=1.5)
dev.off()

png("PM_TAD_size.png", width=2500, height=3500, res=300)
par(lwd=2, mar=c(5.1, 5.1, 4.1, 2.1))
boxplot((Neu_10kb$x2 - Neu_10kb$x1)/1000,
        (Glia_10kb$x2 - Glia_10kb$x1)/1000,
        (Neu_25kb$x2 - Neu_25kb$x1)/1000,
        (Glia_25kb$x2 - Glia_25kb$x1)/1000,
        (Neu_50kb$x2 - Neu_50kb$x1)/1000,
        (Glia_50kb$x2 - Glia_50kb$x1)/1000,
        (Neu_100kb$x2 - Neu_100kb$x1)/1000,
        (Glia_100kb$x2 - Glia_100kb$x1)/1000,
        names= c("10kb", "10kb", "25kb", "25kb", "50kb", "50kb", "100kb", "100kb"),
        col=c("#1b9e77", "#7570b3", "#1b9e77", "#7570b3", "#1b9e77", "#7570b3", "#1b9e77", "#7570b3"),
        ylab="Size of TADs (kb)", ylim=c(0, 7000), cex=1.5, cex.lab=2, cex.axis=1.25)
legend("topleft", rownames(count_table), fill= c("#1b9e77", "#7570b3"), cex=1.5)
dev.off()

