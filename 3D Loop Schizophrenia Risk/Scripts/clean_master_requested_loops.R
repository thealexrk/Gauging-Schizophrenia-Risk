source('myfunctions.R')

df <- read.table('../master_loops_Schahram/full/master_requested_loops', header=TRUE, sep="\t")
df <- clean_qvals(df)

write.table(df, '../master_loops_Schahram/full/clean_master_requested_loops', quote=FALSE, 
            sep="\t", row.names=FALSE)


df <- read.table('../master_loops_Schahram/downsample/master_requested_loops', header=TRUE, sep="\t")
df <- clean_qvals(df)

write.table(df, '../master_loops_Schahram/downsample/clean_master_requested_loops', quote=FALSE, 
            sep="\t", row.names=FALSE)
