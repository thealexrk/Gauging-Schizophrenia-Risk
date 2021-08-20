
is_cis <- function(row) {
  # Determine if significant interaction is in
  # cis or trans
  # Returns: boolean; True if in cis"
  
  chr1 <- sub(":.*", "", row$anchor.bin.coord)
  chr2 <- sub(":.*", "", row$target.bin.coord)
  if(chr1 == chr2){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

get_distance <- function(row) {
  # Returns: bp distance between significant interaction
  
  start <- as.numeric(sub(".*:(.*)-.*", "\\1", row$anchor.bin.coord))
  end <- as.numeric(sub(".*:(.*)-.*", "\\1", row$target.bin.coord))
  d <- abs(start - end)
  return(d)
}

rand_genomic_pos <- function(d, c, cm) {
  # Args: '
  # d = bp distance
  # c = chromosome size table
  # cm = cumulative chromosome size table 
  # Returns: list coords
  #   chrom = chromosome
  #   g1 = random genomic coordinate
  #   g2 = paired genomic coordinate d bp from g1
  
  total <- cm$bp[length(cm$bp)]
  chrStarts <- c(1,cm$bp[1:length(cm$bp) -1] + 1)
  rand_pos <- sample(total, 1)
  chr_idx <- max(which(chrStarts <= rand_pos))
  chrom <- as.character(c$chrom[chr_idx])
  g1 <- rand_pos - chrStarts[chr_idx] + 1
  g2 <- g1 + d
  coords <- list(chrom, g1, g2)
  if (g2 > c$bp[chr_idx]) {
    coords<- rand_genomic_pos(d, c, cm)
  }
  return(coords)
}

overlap <- function(a, b){
  # Return boolean for whether interval b
  # overlaps interval a
  
  if (((a[1] <= b[1]) & (b[1]<= a[2])) | ((b[1] <= a[1]) & (a[1]<= b[2]))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

get_overlap_gene_idxs <- function(RPKMl, rand_coord) {
  # Return row numbers of RPKM_df for genes
  # overlapping input random coordinates
  # Args:
  # RPKMl = list of chromosome dataframes of genes with RPKMs to sample from
  # rand_coord = list of [chrom, coord1, coord2]
  # Returns: 
  # vector idx = row indexes for genes that overlap 
  # random genome coordinates
  
  rand_chrom <- rand_coord[[1]]
  # 10kb bins
  rand_g1 <- c(rand_coord[[2]], rand_coord[[2]] + 10000)
  rand_g2 <- c(rand_coord[[3]], rand_coord[[3]] + 10000)
  idx = c()
  chrom_gene_df <-  RPKMl[[rand_chrom]]
  for (i in 1:nrow(chrom_gene_df)) {
    gene_pos <- c(chrom_gene_df[i,]$start, chrom_gene_df[i,]$end)
    if (overlap(gene_pos, rand_g1) | overlap(gene_pos, rand_g2)) {
      idx <- c(idx, as.numeric(rownames(chrom_gene_df[i,])))
    }
  }
  return(idx)
}

sample_null_distribution <- function(intxns, cSizes, cmcSizes, RPKM_l, RPKM_d, sample_size) {
  # sample genes from RPKM_df that have distances with an equivalent 
  # distribution as distances between genes with significant 
  # interactions with PGC loci
  # Args:
  # intxns = dataframe of significant PGC interactions
  # cSizes = dataframe of bp chromosome sizes
  # cmSizes = dataframe of bp cumulative chromosome sizes
  # RPKM_l = list of chromosome dataframes of genes with RPKMs to sample from
  # RPKM_d = full dataframe of genes with RPKMs to sample from
  # sample_size = number of PGC interacting genes
  # Returns:
  # dataframe of sample genes for null distribution
  
  rand_idxs <- c()
  while(length(rand_idxs) < sample_size) { 
    # Get random PGC interaction
    r <- sample(1:nrow(intxns), 1)
    randRow = intxns[r,]
    if (is_cis(randRow)) {
      # Get distance of interaction
      d <- get_distance(randRow)
      # Select random coordinates with same distance
      rand_coord <- rand_genomic_pos(d, cSizes, cmcSizes)
      #print(rand_coord)
      # Get genes overlapping random coordinate
      rand_idx <- get_overlap_gene_idxs(RPKM_l, rand_coord)
      rand_idxs <- c(rand_idxs, rand_idx)
    } else {
      print("ERROR: trans interaction")
      stop()
    }
    # Check for and remove duplicates
    if(anyDuplicated(rand_idxs) != 0) {
      rand_idxs <- rand_idxs[!duplicated(rand_idxs)]
    }
  }
  # Trim to sample_size in cases it went over
  rand_idxs <- rand_idxs[1:sample_size]
  return(RPKM_d[rand_idxs,])
}




