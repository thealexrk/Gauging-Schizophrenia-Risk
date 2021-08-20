#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description= "Make overlapping qvalue histogram plots for complement and intersections of cell types")
parser$add_argument("-r", help= "radius for merging in kb", type="integer", default=25)
args <- parser$parse_args()


group_qvals <- function(F1_row, F2, F3, radius, nbd) {
  # Get data
  chr_1 <- F1_row$chr1
  chr_2 <- F1_row$chr2
  qval <- as.numeric(F1_row[nbd])
  x1 <- F1_row$x1
  y1 <- F1_row$y1
  F2_pair = FALSE
  F3_pair = FALSE
  # Check if loop intersects with F2 and F3 loops
  if (any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & (abs(x1 - F2$x1) <= radius) & (abs(y1 - F2$y1) <= radius))) {
    F2_pair = TRUE
  }
  if (any((chr_1 == F3$chr1) & (chr_2 == F3$chr2) & (abs(x1 - F3$x1) <= radius) & (abs(y1 - F3$y1) <= radius))) {
    F3_pair = TRUE
  }
  
  if(F2_pair & F3_pair) {
    return(list("group" = "all", "qval" = qval)) 
  }
  else if(xor(F2_pair,F3_pair)) {
    return(list("group" = "pair", "qval" = qval))
  }
  else {
    if(F2_pair | F3_pair) {
      print("ERROR: grouping qvals")
      stop()
    }
    else {
      return(list("group" = "alone", "qval" = qval))
    }
  }
}

Astro <- read.table('../Astro2/merged_loops_30', header=TRUE, sep="\t")
Neu2 <- read.table('../Neu2/merged_loops_30', header=TRUE, sep="\t")
NPC2 <- read.table('../NPC2/merged_loops_30', header=TRUE, sep="\t")

total_Astro <- length(Astro$x1)
total_Neu2 <- length(Neu2$x1)
total_NPC2 <- length(NPC2$x1)


######################################################################################################
# Histograms of qvalues for different neighborhoods in different cell types
# Fix zeros
Astro[Astro$fdrDonut == 0, "fdrDonut"] = 1.0e-42
Neu2[Neu2$fdrDonut == 0, "fdrDonut"] = 1.0e-42
NPC2[NPC2$fdrDonut == 0, "fdrDonut"] = 1.0e-42
Astro[Astro$fdrBL == 0, "fdrBL"] = 1.0e-42
Neu2[Neu2$fdrBL == 0, "fdrBL"] = 1.0e-42
NPC2[NPC2$fdrBL == 0, "fdrBL"] = 1.0e-42
Astro[Astro$fdrH == 0, "fdrH"] = 1.0e-42
Neu2[Neu2$fdrH == 0, "fdrH"] = 1.0e-42
NPC2[NPC2$fdrH == 0, "fdrH"] = 1.0e-42
Astro[Astro$fdrV == 0, "fdrV"] = 1.0e-42
Neu2[Neu2$fdrV == 0, "fdrV"] = 1.0e-42
NPC2[NPC2$fdrV == 0, "fdrV"] = 1.0e-42

png("qvalue_celltypes_hist_donut.png", height=2500, width=1000, res=300)
  par(mfrow=c(3,1))
  hist(log10(Astro$fdrDonut), breaks=seq(log10(min(Astro$fdrDonut)), 0, 1), ylim=c(0,1300), 
       main = "Astro", xlab="log10(q-value)", col="seagreen2")
  hist(log10(Neu2$fdrDonut), breaks=seq(log10(min(Neu2$fdrDonut)), 0, 1), ylim=c(0,1300),
       main = "Neu2", xlab="log10(q-value)", col= "orangered")
  hist(log10(NPC2$fdrDonut), breaks=seq(log10(min(NPC2$fdrDonut)), 0, 1), ylim=c(0,1300), 
       main = "NPC2", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()

png("qvalue_celltypes_hist_BL.png", height=2500, width=1000, res=300)
par(mfrow=c(3,1))
hist(log10(Astro$fdrBL), breaks=seq(log10(min(Astro$fdrBL)), 0, 1), ylim=c(0,1300), 
     main = "Astro", xlab="log10(q-value)", col="seagreen2")
hist(log10(Neu2$fdrBL), breaks=seq(log10(min(Neu2$fdrBL)), 0, 1), ylim=c(0,1300),
     main = "Neu2", xlab="log10(q-value)", col= "orangered")
hist(log10(NPC2$fdrBL), breaks=seq(log10(min(NPC2$fdrBL)), 0, 1), ylim=c(0,1300), 
     main = "NPC2", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()

png("qvalue_celltypes_hist_H.png", height=2500, width=1000, res=300)
par(mfrow=c(3,1))
hist(log10(Astro$fdrH), breaks=seq(log10(min(Astro$fdrH)), 0, 1), ylim=c(0,1700), 
     main = "Astro", xlab="log10(q-value)", col="seagreen2")
hist(log10(Neu2$fdrH), breaks=seq(log10(min(Neu2$fdrH)), 0, 1), ylim=c(0,1700),
     main = "Neu2", xlab="log10(q-value)", col= "orangered")
hist(log10(NPC2$fdrH), breaks=seq(log10(min(NPC2$fdrH)), 0, 1), ylim=c(0,1700), 
     main = "NPC2", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()

png("qvalue_celltypes_hist_V.png", height=2500, width=1000, res=300)
par(mfrow=c(3,1))
hist(log10(Astro$fdrV), breaks=seq(log10(min(Astro$fdrV)), 0, 1), ylim=c(0,1700), 
     main = "Astro", xlab="log10(q-value)", col="seagreen2")
hist(log10(Neu2$fdrV), breaks=seq(log10(min(Neu2$fdrV)), 0, 1), ylim=c(0,1700),
     main = "Neu2", xlab="log10(q-value)", col= "orangered")
hist(log10(NPC2$fdrV), breaks=seq(log10(min(NPC2$fdrV)), 0, 1), ylim=c(0,1700), 
     main = "NPC2", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()


######################################################################################################
# Histograms of complement and intersecting loops for NPC2, Neu2, Astro
# DONUT

donut_list <- list("alone" = c(), "pair" = c(), "all" = c())

for (row_idx in 1:nrow(NPC2)) {
  group_results <- group_qvals(NPC2[row_idx,], Neu2, Astro, args$r * 1000, "fdrDonut")
  grp <- group_results$group
  q <- group_results$qval
  donut_list[[grp]] <- c(donut_list[[grp]], q)
}

for (row_idx in 1:nrow(Neu2)) {
  group_results <- group_qvals(Neu2[row_idx,], NPC2, Astro, args$r * 1000, "fdrDonut")
  grp <- group_results$group
  q <- group_results$qval
  donut_list[[grp]] <- c(donut_list[[grp]], q)
}

for (row_idx in 1:nrow(Astro)) {
  group_results <- group_qvals(Astro[row_idx,], NPC2, Neu2, args$r * 1000, "fdrDonut")
  grp <- group_results$group
  q <- group_results$qval
  donut_list[[grp]] <- c(donut_list[[grp]], q)
}

donut_alone = donut_list$alone 
donut_pair = donut_list$pair
donut_all = donut_list$all

#Fix zeros
donut_alone[donut_alone == 0] = 1.0e-42
donut_pair[donut_pair == 0] = 1.0e-42
donut_all[donut_all == 0] = 1.0e-42

png("qvalue_groups_hist_donut.png", height=2500, width=1000, res=300)
par(mfrow=c(3,1))
hist(log10(donut_alone), breaks=seq(log10(min(donut_alone)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1250), main = "Cell Type Specific", xlab="log10(q-value)", col="seagreen2")
hist(log10(donut_pair), breaks=seq(log10(min(donut_pair)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1250), main = "2 Cell Types", xlab="log10(q-value)", col= "orangered")
hist(log10(donut_all), breaks=seq(log10(min(donut_all)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1250), main = "3 Cell Types", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()

# BL
BL_list <- list("alone" = c(), "pair" = c(), "all" = c())

for (row_idx in 1:nrow(NPC2)) {
  group_results <- group_qvals(NPC2[row_idx,], Neu2, Astro, args$r * 1000, "fdrBL")
  grp <- group_results$group
  q <- group_results$qval
  BL_list[[grp]] <- c(BL_list[[grp]], q)
}

for (row_idx in 1:nrow(Neu2)) {
  group_results <- group_qvals(Neu2[row_idx,], NPC2, Astro, args$r * 1000, "fdrBL")
  grp <- group_results$group
  q <- group_results$qval
  BL_list[[grp]] <- c(BL_list[[grp]], q)
}

for (row_idx in 1:nrow(Astro)) {
  group_results <- group_qvals(Astro[row_idx,], NPC2, Neu2, args$r * 1000, "fdrBL")
  grp <- group_results$group
  q <- group_results$qval
  BL_list[[grp]] <- c(BL_list[[grp]], q)
}

BL_alone = BL_list$alone 
BL_pair = BL_list$pair
BL_all = BL_list$all

#Fix zeros
BL_alone[BL_alone == 0] = 1.0e-42
BL_pair[BL_pair == 0] = 1.0e-42
BL_all[BL_all == 0] = 1.0e-42

png("qvalue_groups_hist_BL.png", height=2500, width=1000, res=300)
par(mfrow=c(3,1))
hist(log10(BL_alone), breaks=seq(log10(min(BL_alone)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1250), main = "Cell Type Specific", xlab="log10(q-value)", col="seagreen2")
hist(log10(BL_pair), breaks=seq(log10(min(BL_pair)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1250), main = "2 Cell Types", xlab="log10(q-value)", col= "orangered")
hist(log10(BL_all), breaks=seq(log10(min(BL_all)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1250), main = "3 Cell Types", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()

# Horizontal
H_list <- list("alone" = c(), "pair" = c(), "all" = c())

for (row_idx in 1:nrow(NPC2)) {
  group_results <- group_qvals(NPC2[row_idx,], Neu2, Astro, args$r * 1000, "fdrH")
  grp <- group_results$group
  q <- group_results$qval
  H_list[[grp]] <- c(H_list[[grp]], q)
}

for (row_idx in 1:nrow(Neu2)) {
  group_results <- group_qvals(Neu2[row_idx,], NPC2, Astro, args$r * 1000, "fdrH")
  grp <- group_results$group
  q <- group_results$qval
  H_list[[grp]] <- c(H_list[[grp]], q)
}

for (row_idx in 1:nrow(Astro)) {
  group_results <- group_qvals(Astro[row_idx,], NPC2, Neu2, args$r * 1000, "fdrH")
  grp <- group_results$group
  q <- group_results$qval
  H_list[[grp]] <- c(H_list[[grp]], q)
}

H_alone = H_list$alone 
H_pair = H_list$pair
H_all = H_list$all

#Fix zeros
H_alone[H_alone == 0] = 1.0e-42
H_pair[H_pair == 0] = 1.0e-42
H_all[H_all == 0] = 1.0e-42

png("qvalue_groups_hist_H.png", height=2500, width=1000, res=300)
par(mfrow=c(3,1))
hist(log10(H_alone), breaks=seq(log10(min(H_alone)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1700), main = "Cell Type Specific", xlab="log10(q-value)", col="seagreen2")
hist(log10(H_pair), breaks=seq(log10(min(H_pair)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1700), main = "2 Cell Types", xlab="log10(q-value)", col= "orangered")
hist(log10(H_all), breaks=seq(log10(min(H_all)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1700), main = "3 Cell Types", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()


# Vertical
V_list <- list("alone" = c(), "pair" = c(), "all" = c())

for (row_idx in 1:nrow(NPC2)) {
  group_results <- group_qvals(NPC2[row_idx,], Neu2, Astro, args$r * 1000, "fdrV")
  grp <- group_results$group
  q <- group_results$qval
  V_list[[grp]] <- c(V_list[[grp]], q)
}

for (row_idx in 1:nrow(Neu2)) {
  group_results <- group_qvals(Neu2[row_idx,], NPC2, Astro, args$r * 1000, "fdrV")
  grp <- group_results$group
  q <- group_results$qval
  V_list[[grp]] <- c(V_list[[grp]], q)
}

for (row_idx in 1:nrow(Astro)) {
  group_results <- group_qvals(Astro[row_idx,], NPC2, Neu2, args$r * 1000, "fdrV")
  grp <- group_results$group
  q <- group_results$qval
  V_list[[grp]] <- c(V_list[[grp]], q)
}

V_alone = V_list$alone 
V_pair = V_list$pair
V_all = V_list$all

#Fix zeros
V_alone[V_alone == 0] = 1.0e-42
V_pair[V_pair == 0] = 1.0e-42
V_all[V_all == 0] = 1.0e-42

png("qvalue_groups_hist_V.png", height=2500, width=1000, res=300)
par(mfrow=c(3,1))
hist(log10(V_alone), breaks=seq(log10(min(V_alone)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1700), main = "Cell Type Specific", xlab="log10(q-value)", col="seagreen2")
hist(log10(V_pair), breaks=seq(log10(min(V_pair)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1700), main = "2 Cell Types", xlab="log10(q-value)", col= "orangered")
hist(log10(V_all), breaks=seq(log10(min(V_all)), 0, 1), xlim=c(-42,0),
     ylim=c(0,1700), main = "3 Cell Types", xlab="log10(q-value)", col = "mediumorchid4")
dev.off()

