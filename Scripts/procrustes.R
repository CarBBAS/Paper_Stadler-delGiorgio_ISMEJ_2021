

dna <- subset_samples(pb, DnaType == 'DNA')
rna <- subset_samples(pb, DnaType == 'cDNA')


dna.mat <- t(otu_mat(dna))
rna.mat <- t(otu_mat(rna))

dna.df <- setDT(as.data.frame(dna.mat), keep.rownames = 'Sample')
rna.df <- setDT(as.data.frame(rna.mat), keep.rownames = 'Sample')

# Get meta data to rename DNA and RNA data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(DnaType), 
                   stringsAsFactors = F)

dna.df <- merge(dna.df, meta, by =  "Sample")
rna.df <- merge(rna.df, meta, by = "Sample")

# correct a few wrong sample names for matching DNA and RNA
uni.name <- llply(list(dna.df, rna.df), function(x){
  setDT(x)
  x$Sample <- str_replace(x$Sample, "RO2R52R","RO2.52R")
  x$Sample <- str_replace(x$Sample, "SWR34R","SW34R")
  x$Sample <- str_replace(x$Sample, "RO2.36pD","RO2.36D")
  x$Sample <- str_replace(x$Sample, "RO2.36pR","RO2.36R")
  x$Sample <- str_replace(x$Sample, "RO2111.60mD", "RO2111.90mD")
  x$Sample <- str_replace(x$Sample, "RO2.30DPR", "RO2.30R") # two DNA
  x$Sample <- str_replace(x$Sample, "RO301.HypoR", "RO31.HypoR")
  x$Sample <- str_replace(x$Sample, "RO301R", "RO31R")
  x$Sample <- str_replace(x$Sample, "RO304R", "RO34R")
  x$Sample <- str_replace(x$Sample, "RO307R", "RO37R")
  x$Sample <- str_replace(x$Sample, "L230R", "L330R")
  
  x[x$DnaType == "DNA", ID := str_replace(x$Sample[x$DnaType == "DNA"], "D$", "")]
  x[x$DnaType == "cDNA", ID := str_replace(x$Sample[x$DnaType == "cDNA"], "R$", "")]
  x[,DnaType := NULL]
  
  return(x)
})

# keep only samples that are in both datasets
uni.name[[1]] <- uni.name[[1]][uni.name[[1]]$ID %in% uni.name[[2]]$ID,]
uni.name[[2]] <- uni.name[[2]][uni.name[[2]]$ID %in% uni.name[[1]]$ID,]

# summarise duplicates
sum <- llply(uni.name, function(x){
  mel.df <- melt.data.table(
    x,
    id.vars = c("ID"),
    measure.vars = patterns("^ASV_"),
    variable.name = "ASV",
    value.name = "reads"
  )
  
  mel.df <- mel.df[, .(reads = mean(reads, na.rm = T)), by = .(ASV, ID)]
  out <- dcast(mel.df, ID ~ ASV, value.var = "reads")
  return(out)
})

pcoa.ls <- llply(sum, function(x){
  pb.mat <- as.matrix(x[,c(2:ncol(x))], rownames = x$ID)
  pb.mat <- log2(pb.mat + 1)
  bray <- vegdist(pb.mat, method = "bray")
  bray <- sqrt(bray) # make euclidean
  
  pcoa <- ape::pcoa(bray)
  
  return(list(Bray = bray, PCOA = pcoa))
})




# convert distance matrix into long format
bray.melt <- melt.dist(bray) %>% mutate(Metric = "Bray")
jacc.melt <- melt.dist(jacc) %>% mutate(Metric = "Jaccard")

# merge both distances
dist.df <- rbind(bray.melt, jacc.melt)

# Get meta data to rename DNA and RNA data
meta <- data.frame(Sample = as.character(row.names(sample_df(physeq))),
                   sample_df(physeq) %>% dplyr::select(DnaType, Year, Season, sample.type.year), 
                   stringsAsFactors = F)

# add meta data for both sample x and sample y
dist.df <- merge(dist.df, meta, by.x =  "Sample.x", by.y = "Sample")
dist.df <- merge(dist.df, meta %>% select(DnaType, Sample, Year), by.x =  "Sample.y", by.y = "Sample")

# omit all samples of 2015 (no RNA was taken, and sample name strategy changed -> creates duplicates)
dist.df <- dist.df[dist.df$Year.x != 2015,]
dist.df <- dist.df[dist.df$Year.y != 2015,]

# correct a few wrong sample names for matching DNA and RNA
dist.df <- as.data.frame(cbind(apply(dist.df[1:2], 2, function(x){
  x <- str_replace(x, "RO2R52R","RO2.52R")
  x <- str_replace(x, "SWR34R","SW34R")
  x <- str_replace(x, "RO2.36pD","RO2.36D")
  x <- str_replace(x, "RO2.36pR","RO2.36R")
  x <- str_replace(x, "RO2111.60mD", "RO2111.90mD")
  x <- str_replace(x, "RO2.30DPR", "RO2.30R") # two DNA
  x <- str_replace(x, "RO301.HypoR", "RO31.HypoR")
  x <- str_replace(x, "RO301R", "RO31R")
  x <- str_replace(x, "RO304R", "RO34R")
  x <- str_replace(x, "RO307R", "RO37R")
  x <- str_replace(x, "L230R", "L330R")
}), dist.df[3:ncol(dist.df)]))

# remove Ds and Rs to match counterpart samples
dist.df$ID.x[dist.df$DnaType.x == "DNA"] <- str_replace(dist.df$Sample.x[dist.df$DnaType.x == "DNA"], "D$", "")
dist.df$ID.x[dist.df$DnaType.x == "cDNA"] <- str_replace(dist.df$Sample.x[dist.df$DnaType.x == "cDNA"], "R$", "")
dist.df$ID.y[dist.df$DnaType.y == "DNA"] <- str_replace(dist.df$Sample.y[dist.df$DnaType.y == "DNA"], "D$", "")
dist.df$ID.y[dist.df$DnaType.y == "cDNA"] <- str_replace(dist.df$Sample.y[dist.df$DnaType.y == "cDNA"], "R$", "")
