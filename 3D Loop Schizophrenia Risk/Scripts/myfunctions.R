clean_qvals <- function(d) {
  print(nrow(d))
  # Remove NA cases
  d <- na.omit(d)
  print(nrow(d))
  print(head(d))
  # Fix zeros
  d[d$AstrofdrDonut == 0, "AstrofdrDonut"] = 1.0e-42
  d[d$GMfdrDonut == 0, "GMfdrDonut"] = 1.0e-42
  d[d$NeufdrDonut == 0, "NeufdrDonut"] = 1.0e-42
  d[d$NPCfdrDonut == 0, "NPCfdrDonut"] = 1.0e-42
  print(head(d))
  # Change q values to log10(qvals)
  d$AstrofdrDonut <- log10(d$AstrofdrDonut)
  d$GMfdrDonut <- log10(d$GMfdrDonut)
  d$NeufdrDonut <- log10(d$NeufdrDonut)
  d$NPCfdrDonut <- log10(d$NPCfdrDonut)
  print(head(d))
  
  # Remove duplicates
  d <- d[!duplicated(d[,c('chr1','x1','x2','chr2','y1','y2')]),]
  print(nrow(d))
  return(d)
  
  
}

