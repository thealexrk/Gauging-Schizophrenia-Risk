library(VennDiagram)
source('C:/cygwin64/home/Tyler/Research/Schahram/Schahram-project/myfunctions.R')

signif_df <- function(d, cells) {
  for (i in 1:length(cells)) {
    d <- subset(d, d[cells[i]] < -1)
  }
  return(d)
}

count_loops_at_res <- function(d, res) {
  loop_count <- sum(abs(d$x1 - d$x2) == res)
  return(loop_count)
  
}

# Full dataset
full_df <- read.table("C:/cygwin64/home/Tyler/Research/Schahram/master_loops_Schahram/full/master_requested_loops", header=TRUE, sep="\t")
full_df <- clean_qvals(full_df)
# Downsampled dataset
ds_df <- read.table("C:/cygwin64/home/Tyler/Research/Schahram/master_loops_Schahram/downsample/master_requested_loops", header=TRUE, sep="\t")
ds_df <- clean_qvals(ds_df)


cell = c('AstrofdrDonut', 'GMfdrDonut', 'NeufdrDonut', 'NPCfdrDonut')


out = data.frame(matrix(nrow=6, ncol=4, 
       dimnames = list(c("5kb_full", "10kb_full", "5kb_10kb_full",
      "5kb_downsample", "10kb_downsample", "5kb_10kb_downsample"),
      c("Astro", "NPC", "Neu", "GM"))))

#Astro
sig_Astro_full <- signif_df(full_df, "AstrofdrDonut")
sig_Astro_ds <- signif_df(ds_df, "AstrofdrDonut")
out["5kb_full", "Astro"] <- count_loops_at_res(sig_Astro_full, 5000)
out["10kb_full", "Astro"] <- count_loops_at_res(sig_Astro_full, 10000)
out["5kb_10kb_full", "Astro"] <- nrow(sig_Astro_full)
out["5kb_downsample", "Astro"] <- count_loops_at_res(sig_Astro_ds, 5000)
out["10kb_downsample", "Astro"] <- count_loops_at_res(sig_Astro_ds, 10000)
out["5kb_10kb_downsample", "Astro"] <- nrow(sig_Astro_ds)

#NPC
sig_NPC_full <- signif_df(full_df, "NPCfdrDonut")
sig_NPC_ds <- signif_df(ds_df, "NPCfdrDonut")
out["5kb_full", "NPC"] <- count_loops_at_res(sig_NPC_full, 5000)
out["10kb_full", "NPC"] <- count_loops_at_res(sig_NPC_full, 10000)
out["5kb_10kb_full", "NPC"] <- nrow(sig_NPC_full)
out["5kb_downsample", "NPC"] <- count_loops_at_res(sig_NPC_ds, 5000)
out["10kb_downsample", "NPC"] <- count_loops_at_res(sig_NPC_ds, 10000)
out["5kb_10kb_downsample", "NPC"] <- nrow(sig_NPC_ds)

#Neu
sig_Neu_full <- signif_df(full_df, "NeufdrDonut")
sig_Neu_ds <- signif_df(ds_df, "NeufdrDonut")
out["5kb_full", "Neu"] <- count_loops_at_res(sig_Neu_full, 5000)
out["10kb_full", "Neu"] <- count_loops_at_res(sig_Neu_full, 10000)
out["5kb_10kb_full", "Neu"] <- nrow(sig_Neu_full)
out["5kb_downsample", "Neu"] <- count_loops_at_res(sig_Neu_ds, 5000)
out["10kb_downsample", "Neu"] <- count_loops_at_res(sig_Neu_ds, 10000)
out["5kb_10kb_downsample", "Neu"] <- nrow(sig_Neu_ds)

#GM
sig_GM_full <- signif_df(full_df, "GMfdrDonut")
sig_GM_ds <- signif_df(ds_df, "GMfdrDonut")
out["5kb_full", "GM"] <- count_loops_at_res(sig_GM_full, 5000)
out["10kb_full", "GM"] <- count_loops_at_res(sig_GM_full, 10000)
out["5kb_10kb_full", "GM"] <- nrow(sig_GM_full)
out["5kb_downsample", "GM"] <- count_loops_at_res(sig_GM_ds, 5000)
out["10kb_downsample", "GM"] <- count_loops_at_res(sig_GM_ds, 10000)
out["5kb_10kb_downsample", "GM"] <- nrow(sig_GM_ds)

write.table(out, "count_loops.txt", quote=FALSE, sep="\t")
