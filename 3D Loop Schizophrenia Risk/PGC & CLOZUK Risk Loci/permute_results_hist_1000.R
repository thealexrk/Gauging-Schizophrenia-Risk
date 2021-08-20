# Read results
NPC_org <- scan("NPC_PGC_CLOZUK_org_scores.txt")
Neu_org <- scan("Neu_PGC_CLOZUK_org_scores.txt")
Glia_org <- scan("Glia_PGC_CLOZUK_org_scores.txt")

plot_org_hist <- function(scores, cell, tstat) {
  
  pdf(paste(cell, "_PGC_CLOZUK_org_scores_distance_constraint_hist_1000perm.pdf",sep=""), height=5, width=6.75)
  par(mar=c(5, 5, 4, 2) + 0.1)
  hist(scores, breaks= seq(0.12, 0.26, by = 0.003), freq=FALSE, xlim=c(0.12, 0.26), xlab = "Organization Score",
       col="lightgrey", main=cell, cex.lab=1.5, cex.main=1.5, ylim=c(2.65, 70))
  abline(v=tstat, col="red", lwd=2)
  box(bty="l")
  dev.off()
}


plot_org_hist(NPC_org, "NPC", 0.2484)
plot_org_hist(Neu_org, "Neu", 0.2446)
plot_org_hist(Glia_org, "Glia", 0.2186)


# p-values
print("NPC p-value: ")
NPC_p <- (sum(NPC_org > 0.2484))/ 1000
print(NPC_p)
print("Neu p-value: ")
Neu_p <- (sum(Neu_org > 0.2446))/ 1000
print(Neu_p)
Glia_p <- (sum(Glia_org > 0.2186))/ 1000
print("Glia p-value: ")
print(Glia_p)
      