#!/usr/bin/env Rscript
library(argparse)
library(VennDiagram)

parser <- ArgumentParser(description= "Count overlap of called loops for 2 merged_loops files output by HICCUPS")
parser$add_argument("-f1", help= "merged_loops file 1", type="character", required=TRUE)
parser$add_argument("-f2", help= "merged_loops file 2", type="character", required=TRUE)
parser$add_argument("-r", help= "radius for merging in kb", type="integer", default=25)
args <- parser$parse_args()


count_overlaps <- function(F1_row, F2, radius) {
	# Get data
	chr_1 <- F1_row$chr1
	chr_2 <- F1_row$chr2
	x1 <- F1_row$x1
	y1 <- F1_row$y1
	# Count overlaps
	overlaps <- sum(any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & ((abs(x1 - F2$x1) <= radius)) & (abs(y1 - F2$y1) <= radius)))
	return(overlaps)
}

f1_overlap_loop <- function(F1_row, F2, radius) {
  # Get data
  chr_1 <- F1_row$chr1
  chr_2 <- F1_row$chr2
  x1 <- F1_row$x1
  y1 <- F1_row$y1
  # Does this loop overlap with loops from F2?
  overlap <- any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & ((abs(x1 - F2$x1) <= radius)) & (abs(y1 - F2$y1) <= radius))
  return(overlap)
}

F1 <- read.table(args$f1, header=TRUE, sep="\t")
F2 <- read.table(args$f2, header=TRUE, sep="\t")
radii <- seq(0, 200000, 5000)

total_F1 <- length(F1$x1)
total_F2 <- length(F2$x1)

counts <- c()

for (radius in radii) {
	count <-  0
	for (row_idx in 1:nrow(F1)) {
		count <- count + count_overlaps(F1[row_idx,], F2, radius)
	}
	counts <- c(counts, count)
}

png("overlap_loops_radius.png", width=2500, height=2500, res=300)
	plot(radii/1000, counts, xlab="radius in (kb)", ylab= "Number of overlapping loops", type="o", pch=20)
dev.off()

count <- 0
for (row_idx in 1:nrow(F1)) {
		count <- count + count_overlaps(F1[row_idx,], F2, args$r * 1000)
	}

png("overlap_loops_venn_25kb.png", width=2500, height=2500, res=300)
draw.pairwise.venn(area1 = total_F1, area2 = total_F2, cross.area = count, category = c("NPC2", 
    "Neu2"))
dev.off()

print("Done with Venn")

########################################################################################
# Write F1 loops that don't overlap
F1_no_overlap <- data.frame()

for (row_idx in 1:nrow(F1)) {
  if (! f1_overlap_loop(F1[row_idx,], F2, args$r * 1000)){
    F1_no_overlap <- rbind(F1_no_overlap, F1[row_idx,])
  }
}
write.table(F1_no_overlap, "NPC_no_overlap.txt", row.names=FALSE, sep="\t")
