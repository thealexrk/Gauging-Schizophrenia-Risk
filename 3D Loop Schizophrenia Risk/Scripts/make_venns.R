#!/usr/bin/env Rscript
library(argparse)
library(VennDiagram)

parser <- ArgumentParser(description= "Make 4 Venn Diagram plots for Astro, NPC2, Neu2, and GM")
parser$add_argument("-r", help= "radius for merging in kb", type="integer", default=25)
args <- parser$parse_args()


count_overlap_pair <- function(F1_row, F2, radius) {
  # Get data
  chr_1 <- F1_row$chr1
  chr_2 <- F1_row$chr2
  x1 <- F1_row$x1
  y1 <- F1_row$y1
  # Count overlaps
  F1_F2_overlap <- any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & ((abs(x1 - F2$x1) <= radius)) & (abs(y1 - F2$y1) <= radius))
  if (F1_F2_overlap) {
    return(1)
  }
  else {
    return(0)
  }
}

count_overlap_all <- function(F1_row, F2, F3, radius) {
  # Get data
  chr_1 <- F1_row$chr1
  chr_2 <- F1_row$chr2
  x1 <- F1_row$x1
  y1 <- F1_row$y1
  # Count overlaps
  F1_F2_overlap <- any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & ((abs(x1 - F2$x1) <= radius)) & (abs(y1 - F2$y1) <= radius))
  F1_F3_overlap <- any((chr_1 == F3$chr1) & (chr_2 == F3$chr2) & ((abs(x1 - F3$x1) <= radius)) & (abs(y1 - F3$y1) <= radius))
  if (F1_F2_overlap & F1_F3_overlap) {
    return(1)
  }
  else {
    return(0)
  }
}


Astro <- read.table('../Astro2/merged_loops_30', header=TRUE, sep="\t")
Neu2 <- read.table('../Neu2/merged_loops_30', header=TRUE, sep="\t")
NPC2 <- read.table('../NPC2/merged_loops_30', header=TRUE, sep="\t")
GM <- read.table ('../GM/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_add_chr.txt', header=TRUE, sep="\t")

total_Astro <- length(Astro$x1)
total_Neu2 <- length(Neu2$x1)
total_NPC2 <- length(NPC2$x1)
total_GM <- length(GM$x1)

######################################################################################################
#Fig 1. NPC2, Neu2, Astro

NPC2_Neu2 <- 0
NPC2_Astro <- 0
Astro_Neu2 <- 0
all_3 <- 0

for (row_idx in 1:nrow(NPC2)) {
  NPC2_Neu2 <- NPC2_Neu2 + count_overlap_pair(NPC2[row_idx,], Neu2, args$r * 1000)
}

for (row_idx in 1:nrow(NPC2)) {
  NPC2_Astro <- NPC2_Astro + count_overlap_pair(NPC2[row_idx,], Astro, args$r * 1000)
}

for (row_idx in 1:nrow(Astro)) {
  Astro_Neu2 <- Astro_Neu2 + count_overlap_pair(Astro[row_idx,], Neu2,  args$r * 1000)
}

for (row_idx in 1:nrow(NPC2)) {
  all_3 <- all_3 + count_overlap_all(NPC2[row_idx,], Neu2, Astro, args$r * 1000)
}
png("Astro2_Neu2_NPC2.png", width=2500, height=2500, res=300)
draw.triple.venn(area1 = total_Astro, area2 = total_Neu2, area3= total_NPC2, 
                 n12=Astro_Neu2, n23=NPC2_Neu2, n13 = NPC2_Astro, n123=all_3,
                 category = c("Astro2", "Neu2", "NPC2"), 
                 fill = c("seagreen2", "orangered", "mediumorchid4"))
dev.off()
############################################################################################################   
######################################################################################################
#Fig 2. NPC2, Neu2, GM

NPC2_Neu2 <- 0
NPC2_GM <- 0
GM_Neu2 <- 0
all_3 <- 0

for (row_idx in 1:nrow(NPC2)) {
  NPC2_Neu2 <- NPC2_Neu2 + count_overlap_pair(NPC2[row_idx,], Neu2, args$r * 1000)
}

for (row_idx in 1:nrow(NPC2)) {
  NPC2_GM <- NPC2_GM + count_overlap_pair(NPC2[row_idx,], GM, args$r * 1000)
}

for (row_idx in 1:nrow(GM)) {
  GM_Neu2 <- GM_Neu2 + count_overlap_pair(GM[row_idx,], Neu2,  args$r * 1000)
}

for (row_idx in 1:nrow(NPC2)) {
  all_3 <- all_3 + count_overlap_all(NPC2[row_idx,], Neu2, GM, args$r * 1000)
}
png("GM_Neu2_NPC2.png", width=2500, height=2500, res=300)
draw.triple.venn(area1 = total_GM, area2 = total_Neu2, area3= total_NPC2, 
                 n12=GM_Neu2, n23=NPC2_Neu2, n13 = NPC2_GM, n123=all_3,
                 category = c("GM", "Neu2", "NPC2"), 
                 fill = c("lightskyblue", "orangered", "mediumorchid4"))
dev.off()
############################################################################################################ 
######################################################################################################
#Fig 3. NPC2, Astro, GM

NPC2_Astro <- 0
NPC2_GM <- 0
GM_Astro <- 0
all_3 <- 0

for (row_idx in 1:nrow(NPC2)) {
  NPC2_Astro <- NPC2_Astro + count_overlap_pair(NPC2[row_idx,], Astro, args$r * 1000)
}

for (row_idx in 1:nrow(NPC2)) {
  NPC2_GM <- NPC2_GM + count_overlap_pair(NPC2[row_idx,], GM, args$r * 1000)
}

for (row_idx in 1:nrow(GM)) {
  GM_Astro <- GM_Astro + count_overlap_pair(GM[row_idx,], Astro,  args$r * 1000)
}

for (row_idx in 1:nrow(NPC2)) {
  all_3 <- all_3 + count_overlap_all(NPC2[row_idx,], Astro, GM, args$r * 1000)
}
png("GM_Astro2_NPC2.png", width=2500, height=2500, res=300)
draw.triple.venn(area1 = total_GM, area2 = total_Astro, area3= total_NPC2, 
                 n12=GM_Astro, n23=NPC2_Astro, n13 = NPC2_GM, n123=all_3,
                 category = c("GM", "Astro2", "NPC2"), 
                 fill = c("lightskyblue", "seagreen2", "mediumorchid4"))
dev.off()
############################################################################################################
############################################################################################################ 
######################################################################################################
#Fig 4. Neu2, Astro, GM

Neu2_Astro <- 0
Neu2_GM <- 0
GM_Astro <- 0
all_3 <- 0

for (row_idx in 1:nrow(Neu2)) {
  Neu2_Astro <- Neu2_Astro + count_overlap_pair(Neu2[row_idx,], Astro, args$r * 1000)
}

for (row_idx in 1:nrow(Neu2)) {
  Neu2_GM <- Neu2_GM + count_overlap_pair(Neu2[row_idx,], GM, args$r * 1000)
}

for (row_idx in 1:nrow(GM)) {
  GM_Astro <- GM_Astro + count_overlap_pair(GM[row_idx,], Astro,  args$r * 1000)
}

for (row_idx in 1:nrow(Neu2)) {
  all_3 <- all_3 + count_overlap_all(Neu2[row_idx,], Astro, GM, args$r * 1000)
}
png("GM_Astro2_Neu2.png", width=2500, height=2500, res=300)
draw.triple.venn(area1 = total_GM, area2 = total_Astro, area3= total_Neu2, 
                 n12=GM_Astro, n23=Neu2_Astro, n13 = Neu2_GM, n123=all_3,
                 category = c("GM", "Astro2", "Neu2"), 
                 fill = c("lightskyblue", "seagreen2", "orangered"))
dev.off()
############################################################################################################