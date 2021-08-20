df <- read.table("clean_master_requested_loops", sep="\t", header=TRUE, quote="", comment.char="")
Neu <- df$observed[df$NeufdrDonut < -1]

png("read_count_at_loop_Neu_hist.png", width=1800, height=1800, res=300)
  hist(Neu, xlab="Number of reads at loop pixel", main = "Neu Loops", breaks=250,
       col = "#1b9e77", xlim=c(0,400))
  abline(v=50, col="red", lty=2)
dev.off()