source('myfunctions.R')

df <- read.table('../master_loops_Schahram/full/master_requested_loops', header=TRUE, sep="\t")
df <- clean_qvals(df)

png("loop_distance.png", height=1500, width=2500, res=300)
par(mar=c(5, 5, 4, 2) + 0.1,lwd=2)
Astro_loops <- df[df$AstrofdrDonut < -1,]
Astro_distances <- Astro_loops$y2 - Astro_loops$x1
Astro_ecdf <- ecdf(log10(Astro_distances))
plot(Astro_ecdf, verticals=TRUE, do.points=FALSE, col="#7570b3", xlim=c(4.3, 8.2), 
     main="", xlab=expression(log[10]* '(bp distance)'),
     col.01line = NULL, cex.lab=1.5, ylab=expression(F[n]*'(x)'))

NPC_loops <- df[df$NPCfdrDonut < -1,]
NPC_distances <- NPC_loops$y2 - NPC_loops$x1
NPC_ecdf <- ecdf(log10(NPC_distances))
plot(NPC_ecdf, verticals=TRUE, do.points=FALSE, col="#d95f02", add=TRUE,
     col.01line = NULL)

Neu_loops <- df[df$NeufdrDonut < -1,]
Neu_distances <- Neu_loops$y2 - Neu_loops$x1
Neu_ecdf <- ecdf(log10(Neu_distances))
plot(Neu_ecdf, verticals=TRUE, do.points=FALSE, col="#1b9e77", add=TRUE,
     col.01line = NULL)

GM_loops <- df[df$GMfdrDonut < -1,]
GM_distances <- GM_loops$y2 - GM_loops$x1
GM_ecdf <- ecdf(log10(GM_distances))
plot(GM_ecdf, verticals=TRUE, do.points=FALSE, col="lightskyblue", add=TRUE,
     col.01line = NULL)

legend("bottomright", c("NPC", "Glia", "Neuron", "GM12878"), 
       col=c("#d95f02", "#7570b3", "#1b9e77", "lightskyblue"),
       lty=rep(1,4), bty="n")

dev.off()


#KS
print("K-S test results:")
print(paste("Two_samples", "p-value", sep="    "))
ksr <- ks.test(Astro_distances, NPC_distances)
print(paste("Astro_NPC", ksr$p.value, sep="    "))
ksr <- ks.test(Astro_distances, Neu_distances)
print(paste("Astro_Neu", ksr$p.value, sep="    "))
ksr <- ks.test(Astro_distances, GM_distances)
print(paste("Astro_GM", ksr$p.value, sep="    "))
ksr <- ks.test(NPC_distances, Neu_distances)
print(paste("NPC_Neu", ksr$p.value, sep="    "))
ksr <- ks.test(NPC_distances, GM_distances)
print(paste("NPC_GM", ksr$p.value, sep="    "))
ksr <- ks.test(Neu_distances, GM_distances)
print(paste("Neu_GM", ksr$p.value, sep="    "))

