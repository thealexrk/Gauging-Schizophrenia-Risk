library(VennDiagram)
source('myfunctions.R')

signif <- function(d, cells) {
  for (i in 1:length(cells)) {
    d <- subset(d, d[cells[i]] < -1)
  }
  return(nrow(d))
}


df <- read.table('../master_loops_Schahram/master_requested_loops', header=TRUE, sep="\t")

df <- clean_qvals(df)

cell = c('AstrofdrDonut', 'GMfdrDonut', 'NeufdrDonut', 'NPCfdrDonut')
png('4wayVenn.png', width=3000, height=2700, res=300)
draw.quad.venn(signif(df, cell[1]), signif(df, cell[2]), signif(df, cell[3]), signif(df, cell[4]),
               signif(df, cell[1:2]), signif(df, cell[c(1,3)]), signif(df, cell[c(1,4)]), 
               signif(df, cell[2:3]), signif(df, cell[c(2,4)]), signif(df, cell[3:4]),
               signif(df, cell[1:3]), signif(df, cell[c(1,2,4)]), signif(df, cell[c(1,3,4)]),
               signif(df, cell[2:4]), signif(df, cell),
               category=c('Astro', 'GM', 'Neu', 'NPC'),
               fill = c('seagreen2', 'lightskyblue', 'orangered', 'mediumorchid4')
)
dev.off()

# Get brain specific loops
bs <- df
for (cell in c('AstrofdrDonut', 'NeufdrDonut', 'NPCfdrDonut')) {
  bs <- subset(bs, bs[cell] < -1)
}

bs <- subset(bs, bs['GMfdrDonut'] > -1)

write.table(bs, 'brain_specific_loops.txt', row.names=FALSE, sep='\t', quote=FALSE)

gm <- df
for (cell in c('AstrofdrDonut', 'NeufdrDonut', 'NPCfdrDonut')) {
  gm <- subset(gm, gm[cell] > -1)
}
gm <- subset(gm, gm['GMfdrDonut'] < -1)
write.table(gm, 'GM_specific_loops.txt', row.names=FALSE, sep='\t', quote=FALSE)