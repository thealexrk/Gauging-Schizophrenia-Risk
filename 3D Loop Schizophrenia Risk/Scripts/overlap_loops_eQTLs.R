df <- read.table("overlap_loops_eQTLs.txt", header=TRUE, sep="\t", comment.char="")

all_eQTL <- sum(df$GENE!= "")
all <- c(all_eQTL, length(df$GENE) - all_eQTL)

GM_eQTL <- sum(df$GMfdrDonut < -1 & df$GENE != "")
GM <- c(GM_eQTL, sum(df$GMfdrDonut < -1) - GM_eQTL)

Astro_eQTL <- sum(df$AstrofdrDonut < -1 & df$GENE != "")
Astro <- c(Astro_eQTL, sum(df$AstrofdrDonut < -1) - Astro_eQTL)

NPC_eQTL <- sum(df$NPCfdrDonut < -1 & df$GENE != "")
NPC <- c(NPC_eQTL, sum(df$NPCfdrDonut < -1) - NPC_eQTL)

Neu_eQTL <- sum(df$NeufdrDonut < -1 & df$GENE != "")
Neu <- c(Neu_eQTL, sum(df$NeufdrDonut < -1) - Neu_eQTL)

x = data.frame(all, GM, Astro, NPC, Neu, row.names= c("eQTL_loops", "non-eQTL_loops"))
x = as.matrix(x)
png("overlap_loops_eQTLs.png", height=2200, width=1700, res=300)
barplot(x, col=c("mediumorchid4", "orange"), legend=rownames(x), args.legend =list(bty="n"))
text(0.75, 400, paste(round((x["eQTL_loops", "all"]/sum(x[,"all"]))*100), '%', sep=""), col="white")
text(0.75 + 1.2, 400, paste(round((x["eQTL_loops", "GM"]/sum(x[,"GM"]))*100), '%', sep=""), col="white")
text(0.75 + 1.2*2, 400, paste(round((x["eQTL_loops", "Astro"]/sum(x[,"Astro"]))*100), '%', sep=""), col="white")
text(0.75 + 1.2*3, 400, paste(round((x["eQTL_loops", "NPC"]/sum(x[,"NPC"]))*100), '%', sep=""), col="white")
text(0.75 + 1.2*4, 400, paste(round((x["eQTL_loops", "Neu"]/sum(x[,"Neu"]))*100), '%', sep=""), col="white")
dev.off()

# Cell type specific
get_GM_eQTL_s <- function(x, gene) {
  if (gene) {
    return(sum(x$GMfdrDonut < -1 & 
                       x$AstrofdrDonut > -1 & 
                       x$NPCfdrDonut > -1 & 
                       x$NeufdrDonut > -1 &
                       x$GENE != ""))
  }
  else {
    return(sum(x$GMfdrDonut < -1 & 
                 x$AstrofdrDonut > -1 & 
                 x$NPCfdrDonut > -1 & 
                 x$NeufdrDonut > -1))
  }
}
get_Astro_eQTL_s <- function(x, gene) {
  if (gene) {
    return(sum(x$GMfdrDonut > -1 & 
                 x$AstrofdrDonut < -1 & 
                 x$NPCfdrDonut > -1 & 
                 x$NeufdrDonut > -1 &
                 x$GENE != ""))
  }
  else {
    return(sum(x$GMfdrDonut > -1 & 
                 x$AstrofdrDonut < -1 & 
                 x$NPCfdrDonut > -1 & 
                 x$NeufdrDonut > -1))
  }
}
get_NPC_eQTL_s <- function(x, gene) {
  if (gene) {
    return(sum(x$GMfdrDonut > -1 & 
                 x$AstrofdrDonut > -1 & 
                 x$NPCfdrDonut < -1 & 
                 x$NeufdrDonut > -1 &
                 x$GENE != ""))
  }
  else {
    return(sum(x$GMfdrDonut > -1 & 
                 x$AstrofdrDonut > -1 & 
                 x$NPCfdrDonut < -1 & 
                 x$NeufdrDonut > -1))
  }
}
get_Neu_eQTL_s <- function(x, gene) {
  if (gene) {
    return(sum(x$GMfdrDonut > -1 & 
                 x$AstrofdrDonut > -1 & 
                 x$NPCfdrDonut > -1 & 
                 x$NeufdrDonut < -1 &
                 x$GENE != ""))
  }
  else {
    return(sum(x$GMfdrDonut > -1 & 
                 x$AstrofdrDonut > -1 & 
                 x$NPCfdrDonut > -1 & 
                 x$NeufdrDonut < -1))
  }
}

GM_eQTL_s <- get_GM_eQTL_s(df, TRUE)
GM_s <- c(GM_eQTL_s, get_GM_eQTL_s(df,FALSE) - GM_eQTL_s)
Astro_eQTL_s <- get_Astro_eQTL_s(df, TRUE)
Astro_s <- c(Astro_eQTL_s, get_Astro_eQTL_s(df,FALSE) - Astro_eQTL_s)
NPC_eQTL_s <- get_NPC_eQTL_s(df, TRUE)
NPC_s <- c(NPC_eQTL_s, get_NPC_eQTL_s(df,FALSE) - NPC_eQTL_s)
Neu_eQTL_s <- get_Neu_eQTL_s(df, TRUE)
Neu_s <- c(Neu_eQTL_s, get_Neu_eQTL_s(df,FALSE) - Neu_eQTL_s)

x = data.frame(GM_s, Astro_s, NPC_s, Neu_s, row.names= c("eQTL_loops", "non-eQTL_loops"))
colnames(x) = c("GM", "Astro", "NPC", "Neu")
x = as.matrix(x)

png("overlap_celltype_specific_loops_eQTLs.png", height=2200, width=1300, res=300)
par(mar=c(5,4,4,5))
barplot(x, col=c("dodgerblue4", "palegreen"), legend=rownames(x), args.legend =list(x=7, y=1800, bty="n"))
text(-.5 + 1.2, 20, paste(round((x["eQTL_loops", "GM"]/sum(x[,"GM"]))*100), '%', sep=""), col="white")
text(-.5 + 1.2*2, 20, paste(round((x["eQTL_loops", "Astro"]/sum(x[,"Astro"]))*100), '%', sep=""), col="white")
text(-.5 + 1.2*3, 20, paste(round((x["eQTL_loops", "NPC"]/sum(x[,"NPC"]))*100), '%', sep=""), col="white")
text(-.5 + 1.2*4, 20, paste(round((x["eQTL_loops", "Neu"]/sum(x[,"Neu"]))*100), '%', sep=""), col="white")
dev.off()

# Boxplot of Number of SNPs per loop for cell type specific loops
GM_snps <- df[df$GMfdrDonut < -1 & 
  df$AstrofdrDonut > -1 & 
  df$NPCfdrDonut > -1 & 
  df$NeufdrDonut > -1 &
  df$GENE != "", "SNP"]

Astro_snps <- df[df$GMfdrDonut > -1 & 
                df$AstrofdrDonut < -1 & 
                df$NPCfdrDonut > -1 & 
                df$NeufdrDonut > -1 &
                df$GENE != "", "SNP"]

NPC_snps <- df[df$GMfdrDonut > -1 & 
                   df$AstrofdrDonut > -1 & 
                   df$NPCfdrDonut < -1 & 
                   df$NeufdrDonut > -1 &
                   df$GENE != "", "SNP"]


Neu_snps <- df[df$GMfdrDonut > -1 & 
                   df$AstrofdrDonut > -1 & 
                   df$NPCfdrDonut > -1 & 
                   df$NeufdrDonut < -1 &
                   df$GENE != "", "SNP"]

get_snp_number <- function(x) {
  y <- strsplit(as.character(x), "|", fixed=TRUE)
  return(length(y[[1]]))
}

GM_snp_number <- sapply(GM_snps, function(x) get_snp_number(x))
Astro_snp_number <- sapply(Astro_snps, function(x) get_snp_number(x))
NPC_snp_number <- sapply(NPC_snps, function(x) get_snp_number(x))
Neu_snp_number <- sapply(Neu_snps, function(x) get_snp_number(x))

png("eQTL_SNP_number_per_loop.png", height = 2200, width=1300, res = 300)
boxplot(GM_snp_number, Astro_snp_number, NPC_snp_number, Neu_snp_number, names = c("GM", "Astro", "NPC", "Neu"), 
        col=c("lightskyblue","#7570b3", "#d95f02", "#1b9e77"), ylab = "Number of eQTL SNPs per eQTL loop")
dev.off()
