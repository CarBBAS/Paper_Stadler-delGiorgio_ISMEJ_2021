

dna <- subset_samples(pb, DnaType == 'DNA')
rna <- subset_samples(pb, DnaType == 'RNA')


dna.mat <- otu_mat(dna)
rna.mat <- otu_mat(rna)

dna.df <- setDT(as.data.frame(dna.mat), keep.rownames = 'Sample')
rna.df <- setDT(as.data.frame(rna.mat), keep.rownames = 'Sample')

# Get meta data to rename DNA and RNA data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(DR.names), 
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
d <- dna.df[dna.df$DR.names %in% rna.df$DR.names,]
r <- rna.df[rna.df$DR.names %in% d$DR.names,]

pcoa.ls <- llply(list(d,r), function(x){
  x <- as.data.frame(x)
  rownames(x) <- x$DR.names
  x$DR.names <- NULL
  pb.mat <- as.matrix(x[,-1])
  
  pb.mat <- log2(pb.mat + 1)
  bray <- vegdist(pb.mat, method = "bray")
  bray <- sqrt(bray) # make euclidean
  
  pcoa <- cmdscale(bray, k = nrow(pb.mat)-1, eig = T)
  #PCoA will return the number of axis given as k, usually df
  # eigenvalues should be returned
  
  return(list(Bray = bray, Pcoa = pcoa))
}, .parallel = T)

# Map RNA on DNA ordination, only first two axes
proc <- procrustes(pcoa.ls[[1]]$Pcoa, pcoa.ls[[2]]$Pcoa, scores = "sites", choices = c(1,2))
summary(proc)

# The sum of squared residuals between scaled and rotated configurations 
# of each matrix is used as a metric of association (m2) (Peres-Neto and Jackson, 2001).
# The m2 metric varies between 0 and 1, and smaller values of m2 indicate stronger concordance between
# dissimilarity matrices.

protest(pcoa.ls[[1]]$Pcoa, pcoa.ls[[2]]$Pcoa, scores = "sites", permutations = 999)
# report 
# procrustes sum of squares m12 squared
# correlation r
# significance p

# get residuals of procrustes
resid <- data.frame(ID = names(residuals(proc)), 
                    Resid = as.vector(residuals(proc)))
setDT(resid)
#merge with meta data
resid[uni.name[[1]], c("sample.type.year",
                        "Season" ) := list(i.sample.type.year,
                                                   i.Season), on = .(ID)]

# add panels for plotting
resid[, panels := "main"]
resid[sample.type.year == "Tributary" |
          sample.type.year == "Lake" |
          sample.type.year == "Riverine \nLakes" |
          sample.type.year == "Sediment", panels := "side"]

# calculate confidence interval and means of sample type and season combinations
resid <- resid[, .(mean =  mean(Resid, na.rm = T),
                       conf.int = conf.int(Resid),
                   std = sd(Resid)),
               by = .(sample.type.year, Season, panels)]

(
  main.b <-
    ggplot(resid[panels == "main", ], aes(
      x = sample.type.year, y = mean, fill = Season
    )) +
    theme_cust(base_theme = "pubr") +
    geom_errorbar(aes(ymin = mean - std, ymax = mean + std, colour = Season),
                  position = position_dodge(0.7), width = 0) +
    geom_jitter(aes(fill = Season), shape = 21, 
                position = position_dodge(0.7), size = 2.5) +
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
    labs(x = "Sample type", 
         y = paste0("Procrustes residuals")) +
    lims(y = c(min(resid$mean - resid$std, na.rm = T),
               max(resid$mean + resid$std, na.rm = T))) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title = element_text(size = 10)
    )
)

# side panel
(
  side.b <-
    ggplot(resid[panels == "side", ], aes(
      x = sample.type.year, y = mean, fill = Season
    )) +
    theme_cust(base_theme = "pubr") +
    geom_errorbar(aes(ymin = mean - std, ymax = mean + std, colour = Season),
                  position = position_dodge(0.7), width = 0) +
    geom_jitter(aes(fill = Season), shape = 21, 
                position = position_dodge(0.7), size = 2.5) +
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
    labs(x = "Sample type", 
         y = paste0("Procrustes residuals")) +
    lims(y = c(min(resid$mean - resid$std, na.rm = T),
               max(resid$mean + resid$std, na.rm = T))) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 8),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
)

ggarrange(main.b, 
          side.b,
          widths = c(3,1),
          ncol = 2, nrow = 1, 
          common.legend = T,
          legend = "top",
          align = "h",
          font.label = list(size = 10))

# calculate means and confidence intervals by habitat
resid[, .(mean)]

resid <- merge(resid, uni.name[[1]] %>% select(ID, sample.type.year, Season))

setDT(resid)

plot(proc, kind = 1)
# residuals indicate the disagreement of points after rotation and adjustment
plot(proc, kind = 2)

biplot(pcoa.ls[[2]]$PCOA)



t <- data.frame(fit.1 = fitted(proc)[,1], fit.2 = fitted(proc)[,2],
           dna.1 = pcoa.ls[[1]]$Pcoa$points[,1], pcoa.ls[[1]]$Pcoa$points[,2],
           rna.1 = pcoa.ls[[2]]$Pcoa$points[,1], pcoa.ls[[2]]$Pcoa$points[,2])

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
