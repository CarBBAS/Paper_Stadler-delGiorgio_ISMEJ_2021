###-----------------------------###
#-   Multivariate exploration   - #
###-----------------------------###

# This script will execute an overall multivariate analysis to get a first grasp of the data set
#setwd("/media/shared/Documents/University/PhD/Analyses/Molecular/lr.chapter1")

#----------#
# PACKAGES #
#----------#
library(phyloseq)
library(tidyverse)
library(data.table)
library(plyr)
library(ggpubr) # arrange ggplots
library(ggnewscale) # two aesthetic scales
library(plotly)
library(doMC) # parallel computing
library(vegan) #vegdist, diversity
library(ape) #pcoa
library(ade4) #is.euclid

#-----------#
# FUNCTIONS #
#-----------#
source("./Functions/custom_fun.R")

#-----------------#
# PARALLEL SET-UP #
#-----------------#
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# prepare for parallel processing (with 'parallel') for 'vegan' functions
numCores <- detectCores()
cl <- makeCluster(numCores, type = "FORK") # using forking

#-----------#
# 2015-2017 #
#-----------#
#################################################################################################################
#################################################################################################################
#-----------#
# Read data #
#-----------#

# do we have several files per object? -> take newest version
# ASV CSS transformed table
asv.tab <- select_newest("./Output", "201520162017_CSS_asvtab")
asv.tab <- read.csv(
  paste0("./Output/", asv.tab),
  sep = "\t",
  dec = ".",
  stringsAsFactors = F
)

# transpose back to ASV in cols, samples in rows
row.names(asv.tab) <- asv.tab$Taxa.and.Samples
asv.tab[, "Taxa.and.Samples"] <- NULL
asv.tab <- as.matrix(asv.tab)
# row orders need to match between tax.tab and asv.tab
asv.tab <- asv.tab[order(row.names(asv.tab)),]

# Taxonomy table
tax.tab <- select_newest("./Output", "201520162017_tax_table")
tax.tab <-
  as.matrix(read.csv(
    paste0("./Output/", tax.tab),
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ))
# orders need to match between tax.tab and asv.tab
tax.tab <- tax.tab[order(row.names(tax.tab)),]

# Meta data
met.df <-
  select_newest(path = "./Output", file.pattern = "201520162017_meta_data")
met.df <-
  read.csv(
    paste0("./Output/", met.df),
    sep = ";",
    dec = ".",
    stringsAsFactors = F
  )

# phyloseq needs the sample names of the meta data to be the same as the microbial data
met.df <- sample_data(met.df)

# Assign rownames to be Sample ID's
rownames(met.df) <- met.df$DadaNames

# Construct phyloseq object
pb <- phyloseq(otu_table(asv.tab, taxa_are_rows = T),
               sample_data(met.df),
               tax_table(tax.tab))




# create colour vector for later plotting
# ensure consistent colours for all sample types
met.df$sample.type.year <- factor(met.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheicwater", 
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver","RO3", "RO2", "RO1","Deep",
                                                                      "Downriver",
                                                                      "Marine"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Hyporheicwater", 
                                             "Groundwater","Stream", "Tributary",
                                             "Headwater \nLakes", "Upstream \nPonds", "Lake", "Lake",
                                             "Upriver","RO3","RO2", "RO1","Hypolimnion","Downriver",
                                             "Estuary"))

sample.factors <- levels(met.df$sample.type.year)
# create colour vector for plotting
colvec <- c("red4","chocolate3","orangered2","orange3",
            "khaki","cadetblue","darksalmon",
            "darkolivegreen","darkolivegreen3","gold",
            "royalblue","mediumorchid4", "violet","palevioletred2","navy","skyblue",
            "seagreen3")
names(colvec) <- as.character(sample.factors)



############
# Analysis #
############
####################################################################################
#--------------------------#
#- Mulativariate analysis -#
#--------------------------#
# We are dealing with large environmental gradients, thus we expect a high proportion of zeros
# Thus we need to select an asymmetrical similarity distance coefficient
# If we're dealing with samples from failry homogeneous environmental conditions (short envir. gradients)
# and we expect few zeros and symmetric association coefficients we can use the Euclidean distance.
# Logarithm transformation log2(x + 1) (in case of microbial data) can be used to make asymmetric species
# distributions less asymmetric
# Transformations can be used to make asymmetric species distributions more symmetric so that e.g.
# Euclidean distances can be used

# Several methods were compared:
# PCA
# CA
# PCoA
# NMDS - no convergence

# PCoA was chosen for all ordinations for consistency and smallest horse-shoe effect

rel.df <- select_newest("./Objects", "201520162017_css")
rel.df <- readRDS(
  paste0("./Objects/", rel.df))

rel.df$sample.type.year <- factor(rel.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheicwater", 
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver","RO3", "RO2", "RO1","Deep",
                                                                      "Downriver",
                                                                      "Marine"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Hyporheicwater", 
                                             "Groundwater","Stream", "Tributary",
                                             "Headwater \nLakes", "Upstream \nPonds", "Lake", "Lake",
                                             "Upriver","RO3","RO2", "RO1","Hypolimnion","Downriver",
                                             "Estuary"))

setDT(rel.df)
# calculate means by sample type
means <- rel.df[, .(mean.css = mean(css.reads, na.rm = T),
                    sd.css = sd(css.reads, na.rm = T)), by = .(sample.type.year, DnaType, ASV)]

dnarna <- dcast(means, ASV ~ sample.type.year + DnaType, value.var = "mean.css")
setDF(dnarna)
row.names(dnarna) <- dnarna$ASV
dnarna <- as.matrix(dnarna[,-1])

# make mini meta data
meta <- data.frame(ID = colnames(dnarna)) %>%
  separate(ID, into = c("sample.type", "DnaType"), sep = "_", remove = F)
row.names(meta) <- meta$ID

# Construct phyloseq object
pb <- phyloseq(otu_table(dnarna, taxa_are_rows = T),
               sample_data(meta),
               tax_table(tax.tab))

pb <- prune_taxa(names(sort(taxa_sums(pb), TRUE)[1:5000]), pb)
plot_heatmap(pb)

# try to make gradient heatmap
# our main gradient is sample types
sub.pb <- subset_samples(pb, DnaType == "DNA")
sub.pb <- prune_taxa(names(sort(taxa_sums(sub.pb), TRUE)[1:1000]), sub.pb)
plot_heatmap(sub.pb)

rel.df <- select_newest("./Objects", "201520162017_css")
rel.df <- readRDS(
  paste0("./Objects/", rel.df))

rel.df$sample.type.year <- factor(rel.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheicwater", 
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver","RO3", "RO2", "RO1","Deep",
                                                                      "Downriver",
                                                                      "Marine"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Hyporheicwater", 
                                             "Groundwater","Stream", "Tributary",
                                             "Headwater \nLakes", "Upstream \nPonds", "Lake", "Lake",
                                             "Upriver","RO3","RO2", "RO1","Hypolimnion","Downriver",
                                             "Estuary"))

setDT(rel.df)
# calculate means by sample type
means <- rel.df[, .(mean.css = mean(css.reads, na.rm = T),
                    sd.css = sd(css.reads, na.rm = T)), by = .(sample.type.year, DnaType, ASV)]
real.mean <- rel.df[, .(mean.css = mean(css.reads, na.rm = T)), by = .(DnaType, ASV)]
real.mean <- real.mean[DnaType == "DNA",]

mean.order <- real.mean[order(mean.css, decreasing = T)]$ASV

# only do this exercise with DNA first
means <- means[DnaType == "DNA",]

t <- means[order(match(ASV, mean.order))]

sub <- levels(factor(means$ASV))[1:10000][order(match(levels(factor(means$ASV))[1:10000], mean.order))]

p <- ggplot(means[ASV %in% sub,], aes(x = sample.type.year, y = ASV)) +
  geom_tile(aes(fill = log2(mean.css + 1))) +
  scale_fill_viridis_c()

# transpose back to community matrix
com.mat <- dcast(means, sample.type.year ~ ASV, value.var = "mean.css")
rown <- as.character(com.mat$sample.type.year)
com.mat <- setDF(com.mat[,-1])
row.names(com.mat) <- rown

com.mat <- as.matrix(com.mat)
library(spaa)

niches <- niche.width(com.mat, method = "levins")
niches <- niches[1,]

############
#---------------------------------------------------------------------------------------------------#
############
# Only DNA #
############
# subset only DNA samples
dna <- subset_samples(pb, DnaType == "DNA")

# extract ASV matrix
pb.mat <- t(otu_mat(dna))

# 1. PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # FALSE

pb.bray.pcoa <- ape::pcoa(pb.bray)
  
# 2. PCoA with Jaccard (presence-absence)
# less obvious horse-shoe effect
pb.jac <- vegdist(pb.mat, method = "jaccard", binary = T)
is.euclid(pb.jac) # TRUE
pb.jac.pcoa <- ape::pcoa(pb.jac)

# plot
dna.bray <- plot_bray_n_jacc(bray = pb.bray.pcoa,
                 jacc = pb.jac.pcoa,
                 colours = colvec,
                 plot.name = "DNA_bray_jac",
                 output = T)
dna.bray[["bray"]]
# no major difference in bray curtis and jaccard, go only with Bray Curtis
# Jaccard weaker horse-shoe, keep Bray Curtis for consistency

# Try 3D plot
plot.df <- dna.bray[["df"]]
plot.df <- plot.df[plot.df$Metric == "Bray",]
p <- plot_ly(type = "scatter", mode = "markers")

p <- plot_ly(plot.df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, color = ~sample.type.year,
             colors = colvec, size = 5, symbol = ~DnaType, symbols = c(21,22))
p <- p %>% add_markers()
p <- p %>% layout(scene = list(xaxis = 
                                 list(title = paste("PC1 [", unique(plot.df$x), "%]")),
                               yaxis = 
                                 list(title = paste("PC2 [", unique(plot.df$y), "%]")),
                               zaxis = 
                                 list(title = paste("PC3 [", unique(plot.df$z), "%]"))))
p

htmlwidgets::saveWidget(as_widget(p), "PCoA_DNA_3D.html")

# Test for significant difference between factors
ord.df <- dna.bray[["df"]]
adonis(pb.mat ~ sample.type.year + Season, data = ord.df[ord.df$Metric == "Bray",], 
       sqrt.dist = T, method = "bray")
#Permutation: free
#Number of permutations: 999
#Terms added sequentially (first to last)
#                 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#sample.type.year  16    50.914  3.1821  14.219 0.36423  0.001 ***
#  Season             2     4.724  2.3621  10.555 0.03380  0.001 ***
#  Residuals        376    84.147  0.2238         0.60198           
#Total            394   139.784                 1.00000           

############
#---------------------------------------------------------------------------------------------------#
############
# Only RNA #
############
rna <- subset_samples(pb, DnaType == "cDNA")

# extract ASV matrix
pb.mat <- t(otu_mat(rna))

# 1. PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray)
is.euclid(pb.bray) # TRUE
pb.bray.pcoa <- ape::pcoa(pb.bray)

# 2. PCoA with Jaccard (presence-absence)
# no horse-shoe effect
pb.jac <- vegdist(pb.mat, method = "jaccard", binary = T)
is.euclid(pb.jac) # TRUE
pb.jac.pcoa <- ape::pcoa(pb.jac)

# plot
rna.bray <- plot_bray_n_jacc(bray = pb.bray.pcoa,
                 jacc = pb.jac.pcoa,
                 colours = colvec,
                 plot.name = "RNA_brayjac",
                 output = T)
rna.bray[["bray"]]
# no major difference in bray curtis and jaccard, both horse-shoe and arch effect
# splits in sample types and seasons, with summer and autumn cluster

####################
#---------------------------------------------------------------------------------------------------#
####################
# Both DNA and RNA #
####################
# extract species table with species in columns
pb.mat <- t(otu_mat(pb))

# 1. PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE

pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray)
pb.bray.pcoa <- ape::pcoa(pb.bray)

# 2. PCoA with Jaccard (presence-absence)
pb.jac <- vegdist(pb.mat, method = "jaccard", binary = T)
is.euclid(pb.jac) # TRUE
pb.jac.pcoa <- ape::pcoa(pb.jac)

# plot
dnarna.bray <- plot_bray_n_jacc(bray = pb.bray.pcoa,
                 jacc = pb.jac.pcoa,
                 colours = colvec,
                 plot.name = "DNARNA_brayjac",
                 output = T)
dnarna.bray[["bray"]]

# Try 3D plot
plot.df <- dnarna.bray[["df"]]
plot.df <- plot.df[plot.df$Metric == "Bray",]
p <- plot_ly(type = "scatter", mode = "markers")

p <- plot_ly(plot.df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, color = ~Season,
             size = 5, symbol = ~DnaType, symbols = c(21,22))
p <- p %>% add_markers()
p <- p %>% layout(scene = list(xaxis = 
                                 list(title = paste("PC1 [", unique(plot.df$x), "%]")),
                               yaxis = 
                                 list(title = paste("PC2 [", unique(plot.df$y), "%]")),
                               zaxis = 
                                 list(title = paste("PC3 [", unique(plot.df$z), "%]"))))
p

htmlwidgets::saveWidget(as_widget(p), "PCoA_DNARNA_Season_3D.html")

# statistically test if factors are different
# extract ordination coordinates and factors
ord.df <- dnarna.bray[["df"]]
# run PERMANOVA
dnarna.perm <- adonis(pb.mat ~ sample.type.year + Season + DnaType, data = ord.df[ord.df$Metric == "Bray",], 
       sqrt.dist = T, method = "bray", parallel = cl)
# Permutation: free
#Number of permutations: 999
#Terms added sequentially (first to last)
#                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  sample.type.year  16    64.462  4.0289  17.187 0.30202  0.001 ***
#  Season             2     7.436  3.7182  15.861 0.03484  0.001 ***
#  DnaType            1     6.278  6.2783  26.783 0.02942  0.001 ***
#  Residuals        577   135.257  0.2344         0.63372           
#Total            596   213.434                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

####################################################################################
#-------------------------------------------------------#
# Extract pair-wise dissimilarity among DNA-RNA samples #
#-------------------------------------------------------#
# UPDATE: Instead of using distance between DNA and RNA wihin PCoA space, we use dissimilarity
# This approach is favoured as PCoA reduces sample complexity into two dimensions and misses variation
# Pair-wise dissimilarity is a more direct approach to grasp how different DNA and RNA of the same sample are

dist.dr <- dissim.dnarna(pb, save.name = "All", output = T)

#-------------------------------------------------------------#
# Calculate distance between DNA and RNA points in PCoA space #
#-------------------------------------------------------------#
# use custom function to correct a few wrong sample names and match DNA-RNA counterpart samples
# calculating distance between points in two-dimensional space for both Bray-Curtis and Jaccard
# Use three distances, as variances explained by second and third axes are almost identical
dist.dr <- dist.dnarna(dnarna.bray[["df"]], save.name = "3D", dimensions = 3)
#dist.dr <- dist.dnarna(dnarna.bray[["df"]], save.name = "All_2D", dimensions = 2)

# dissimilarity and distance show different patterns....
# what is behind this difference?

####################################################################################
#-----------------------------------------------#
# Same exercise with different abundance groups #
#-----------------------------------------------#
# Define abundance groups
rel.df <- select_newest("./Objects", "201520162017_css")
rel.df <- readRDS(
  paste0("./Objects/", rel.df))

# Conventional grouping does not work with this data set
# Otherwise, there are no abundant ASVs as we cover too different ecosystem types

# Try to come up with a new classification
# 1. Calcaulte the mean abundance of each ASV for a sample type (e.g. reservoir, lake, stream, soil etc)
# 2. Define abundance groups based on rank abundance curves:
#-- * For each sample type, we create a rank abundance curve
#-- * Log-transform and take the derivative of the curve
#-- * Use second derivative of log-curve to define where the curve starts to bend.

# calculate the mean abundance of each ASV for all sample types
setDT(rel.df)
# calculate means by sample type
means <- rel.df[, .(mean.css = mean(css.reads, na.rm = T),
                    sd.css = sd(css.reads, na.rm = T)), by = .(sample.type.year, DnaType, ASV)]
# order the abundances to make ranks
means <- means[mean.css != 0,] # remove all 0 observations
means <- means[order(mean.css, decreasing = T)]
means[, rank.abun := 1:.N, by = .(DnaType, sample.type.year)]
means[, log.mean := log(mean.css), by = .(DnaType, sample.type.year)]

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
# Get derivatives for an example and plot for supplementary material
x <- means[sample.type.year == "Soilwater" & DnaType == "DNA"]

spl <- smooth.spline(x$rank.abun, x$log.mean, spar = 0.5)
pred <- predict(spl)
first <- predict(spl, deriv = 1) # first derivative
sec <- predict(spl, deriv = 2) # second derivative

options(scipen = 999) # avoid scientific annotations
raw <- ggplot() +
  theme_pubr() +
  annotate(xmax = x$rank.abun[localMinima(sec$y)[1]],
           xmin = x$rank.abun[localMaxima(sec$y)[1]],
           ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
  annotate(xmax = x$rank.abun[localMaxima(sec$y)[2]],
           xmin = x$rank.abun[localMaxima(sec$y)[1]],
           ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "forestgreen") +
  annotate(xmax = max(x$rank.abun),
           xmin = x$rank.abun[localMaxima(sec$y)[2]],
           ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
  geom_line(aes(x = x$rank.abun, y = x$mean.css)) +
  geom_point(aes(x = x$rank.abun[localMaxima(sec$y)[1:2]],
                 y = x$mean.css[localMaxima(sec$y)[1:2]]), colour = "tomato") +
  geom_point(aes(x = x$rank.abun[localMinima(sec$y)[1]],
                 y = x$mean.css[localMinima(sec$y)[1]]), colour = "royalblue") +
  labs(x = "", y = "CSS transformed mean read number") +
  theme(axis.title = element_text(size = 9))

logged <- ggplot() +
  theme_pubr() +
  annotate(xmax = pred$x[localMinima(sec$y)[1]],
           xmin = pred$x[localMaxima(sec$y)[1]],
           ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "tomato") +
  annotate(xmax = pred$x[localMaxima(sec$y)[2]],
           xmin = pred$x[localMaxima(sec$y)[1]],
           ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "forestgreen") +
  annotate(xmax = max(pred$x),
           xmin = pred$x[localMaxima(sec$y)[2]],
           ymin = -Inf, ymax = Inf, geom = "rect", alpha = 0.2, fill = "aquamarine") +
  geom_line(aes(x = pred$x, y = pred$y)) +
  geom_point(aes(x = pred$x[localMaxima(sec$y)[1:2]],
                 y = pred$y[localMaxima(sec$y)[1:2]]), colour = "tomato") +
  geom_point(aes(x = pred$x[localMinima(sec$y)[1]],
                 y = pred$y[localMinima(sec$y)[1]]), colour = "royalblue") +
  labs(x = "", y = "CSS transformed mean read number (log)") +
  theme(axis.title = element_text(size = 9))

deriv <- ggplot() +
  theme_pubr() +
  geom_line(aes(x = sec$x, y = sec$y * 1000)) +
  geom_point(aes(x = sec$x[localMaxima(sec$y)[1:2]],
                 y = sec$y[localMaxima(sec$y)[1:2]]* 1000), colour = "tomato") +
  geom_point(aes(x = sec$x[localMinima(sec$y)[1]],
                 y = sec$y[localMinima(sec$y)[1]]* 1000), colour = "royalblue") +
  labs(x = "", y = expression(paste("Acceleration x10"^3, " (2"^"nd", " derivative)"))) +
  theme(axis.title = element_text(size = 9))

p <- ggarrange(raw, logged, deriv, ncol = 3, labels = "auto")
# add x axis title to be in the middle of two panels
(p <- annotate_figure(p, bottom = text_grob("ASV Rank")))
ggsave("./Figures/Final/abundance_class_ex.png", p,
       width = 22, height = 9, unit = "cm")

rm(x, p, raw, logged, pred, sec, first, spl)

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#####################
# Apply to the data #
#####################

classif.thres <- ddply(means, .(DnaType, sample.type.year), function(x){
  spl <- smooth.spline(x$rank.abun, x$log.mean, spar = 0.5)
  #pred <- predict(spl)
  #first <- predict(spl, deriv = 1) # first derivative
  sec <- predict(spl, deriv = 2) # second derivative
  setDT(x)
  x[rank.abun <= x$rank.abun[localMaxima(sec$y)[1]], ab.group := "Abundant", by = .(ASV)]
  x[rank.abun > x$rank.abun[localMaxima(sec$y)[1]] &
      rank.abun <= x$rank.abun[localMaxima(sec$y)[2]], ab.group := "Medium", by = .(ASV)]
  x[rank.abun > x$rank.abun[localMaxima(sec$y)[2]], ab.group := "Rare", by = .(ASV)]
  
  out <- x[, .(max = max(mean.css),
               min = min(mean.css)), by = .(ab.group)]
  return(out)
}, .parallel = T)

setDT(classif.thres)
classif.thres[, .(mean.max = mean(max),
                  sd.max = sd(max),
                  mean.min = mean(min),
                  sd.min = sd(min)), by = .(DnaType, ab.group)]
#saveRDS(classif.thres, "./Objects/abundance.classification.threshold.rds")

####################################################################################################

# After consulting the means and deviations of the maximum and minimum thresholds across samples
# we settle with:
# Abundant >= 100 css reads
# Medium < 100 & >= 20 css reads
# Rare < 20 css reads

means[mean.css >= 100, ab.group := 1] # Abundant
means[mean.css < 100 & mean.css >= 20, ab.group := 2] # Medium
means[mean.css < 20, ab.group := 3] # Rare

# code to resusicate absent rows
temp <- dcast(means, DnaType + ASV ~ sample.type.year, value.var = c("ab.group")) # simple wide first
temp[is.na(temp)] <- 0 # much faster to overwrite NAs this way
ag.class <- melt(temp, id.vars = c("DnaType", "ASV"),
                 variable.name = "sample.type.year", value.name = "ab.group") # make long again

# Add ecosystem domain identifier
ag.class[, ecosys.domain := "Freshwater"]
ag.class[sample.type.year == "Soil" | sample.type.year == "Soilwater" | sample.type.year == "Hyporheicwater"
         | sample.type.year == "Wellwater" | sample.type.year == "Sediment",
         ecosys.domain := "Soily"]
ag.class[sample.type.year == "Marine", ecosys.domain := "Marine"]

# Custom function to identify ASVs for each group
# 1. Abundant in all sample types never absent = Universally abundant
# 2. Medium in all sample types never absent = Universally medium
# 3. Always rare = Universally rare / absent (Assumption: what is absent is so rare that it's below detection limit)
# 4. Abundant in all of one ecosystem domain (soily, fresh or marine) and never rare or absent = Cosmopolitan

abun.list <- dlply(ag.class, .(DnaType), function(x){
  setDT(x)
  # 1. Universally abundant
  univ.abun <- as.character(x[ab.group == 1, abundant.counts := .N, by = .(ASV)][abundant.counts == length(levels(factor(sample.type.year))),]$ASV)
  # 2. Universally medium
  univ.med <- as.character(x[ab.group == 2, abundant.counts := .N, by = .(ASV)][abundant.counts == length(levels(factor(sample.type.year))),]$ASV)
  # 3. Universally rare
  univ.rare <- as.character(unique(x[ab.group == 3 | ab.group == 0, 
                                     rare.counts :=  .N, 
                                     by = .(ASV)][rare.counts == length(levels(factor(x$sample.type.year))),]$ASV))
  # 4. Cosmopolitan
  x[, domain.count := .N, by = .(ASV, ecosys.domain)]
  abun.dom <- unique(as.character(x[ab.group == 1, abundant.counts.by.dom := .N, by = .(ecosys.domain, ASV)][abundant.counts.by.dom == domain.count,]$ASV))
  rare <- as.character(unique(x[ab.group == 3 | ab.group == 0,]$ASV))
  cosmo <- as.character(unique(x[ASV %in% abun.dom & !(ASV %in% rare),]$ASV))
  
  # 5. Abundant in all of one ecosystem domain, but never abundant in other domains == Specialist
  abun.dom <- x[ab.group == 1 | ab.group == 2,
                abun.counts.by.dom := .N, by = .(ecosys.domain, ASV)][abun.counts.by.dom == domain.count, abun.domain := ecosys.domain, by = .(ASV)]
  spec.asvs <- as.character(unique(abun.dom[!is.na(abun.domain),]$ASV))
  abun.dom[ASV %in% spec.asvs & is.na(abun.domain), abun.domain := "No" ]
  abun.dom <- as.character(unique(abun.dom[ASV %in% spec.asvs & 
                                             abun.domain != ecosys.domain &
                                             (ab.group == 1 | ab.group == 2),
                                           abun.counts := .N, by = .(ASV, ecosys.domain)][ASV %in% spec.asvs & abun.domain != ecosys.domain, .(sum.abun.counts = sum(abun.counts, na.rm = T)), by = .(ASV)][sum.abun.counts == 0,]$ASV))
  
  # 6. Shifters
  rest <- c(univ.abun, univ.med, univ.rare, cosmo, abun.dom)
  remain <- x[ASV %in% unique(x$ASV)[!(unique(x$ASV) %in% rest)]]
  
  remain[,
         abundant.counts := nrow(.SD[ab.group == 2 | ab.group == 3]), by = .(ASV)]
  remain[, rare.counts := nrow(.SD[ab.group == 1 | ab.group == 0]), by = .(ASV)]
  
  # Mostly abundant, sometimes rare
  abundant.shifter <- as.character(unique(remain[abundant.counts >= rare.counts,]$ASV))
  # Mostly rare, sometimes abundant
  rare.shifter <- as.character(unique(remain[abundant.counts < rare.counts,]$ASV))
  
  #length(univ.abun) + length(univ.med) + length(univ.rare) + length(cosmo) + length(rare.dom) + length(abundant.shifter) + length(rare.shifter) == length(unique(x$ASV))
  list(universal.abundant = univ.abun,
       universal.medium = univ.med,
       universal.rare = univ.rare,
       cosmopolitan = cosmo,
       specialist = abun.dom,
       abundant.shifter = abundant.shifter,
       rare.shifter = rare.shifter)
}, .parallel = T)

# Only extract abundance classification based on DNA
dna.ab.group <- abun.list[["DNA"]]

# only keep those abundance groups that have ASVs assigned
ls.asvs <- dna.ab.group[which(sapply(dna.ab.group, length) != 0L)]
grp.names <- names(ls.asvs)


# extract raw DNA - RNA relationship
rel.df <- rel.df[Year != 2015,]

# correct a few wrong sample names for matching DNA and RNA
# correct a few wrong sample names for matching DNA and RNA
rel.df[Sample == "RO2R52R", Sample := "RO2.52R"]
rel.df[Sample == "SWR34R", Sample := "SW34R"]
rel.df[Sample == "RO2.36pD", Sample :="RO2.36D"]
rel.df[Sample == "RO2.36pR", Sample :="RO2.36R"]
rel.df[Sample == "RO2111.60mD", Sample := "RO2111.90mD"]
rel.df[Sample == "RO2.30DPR", Sample := "RO2.30R"] # two DNA
rel.df[Sample == "RO301.HypoR", Sample := "RO31.HypoR"]
rel.df[Sample == "RO301R", Sample := "RO31R"] 
rel.df[Sample == "RO304R", Sample := "RO34R" ]
rel.df[Sample == "RO307R", Sample := "RO37R"] 
rel.df[Sample == "L230R", Sample := "L330R"] # L230 does not exist

# remove Ds and Rs to match counterpart samples
rel.df[DnaType == "DNA", ID := str_replace(rel.df$Sample[rel.df$DnaType == "DNA"], "D$", "")]
rel.df[DnaType == "cDNA", ID := str_replace(rel.df$Sample[rel.df$DnaType == "cDNA"], "R$", "")]

# calculate mean coordinates for duplicates
sum <- rel.df[, .(mean.css = mean(css.reads),
                  n = .N), by = .(sample.type.year, DnaType, ASV)]

sum <- sum[mean.css != 0,]

cast.sum <- dcast(sum, sample.type.year + ASV ~ DnaType, value.var = "mean.css")

cast.sum <- cast.sum[!(DNA == 0 & cDNA == 0),]
cast.sum <- cast.sum[!(DNA == 0 | cDNA == 0),]

# add into means dataframe
for(i in 1:length(grp.names)){
  cast.sum[ASV %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
  sum[ASV %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
}

ggplot(cast.sum, aes(x = log2(cDNA +1), y = log2(DNA+1))) +
  geom_point(aes(fill = sample.type.year), shape = 21) +
  geom_smooth(aes(group = ab.names, colour = ab.names), method = "lm", se = F) +
  facet_grid(.~ab.names)

ggplot(sum, aes(x = sample.type.year, y = log2(mean.css +1))) +
  geom_violin(aes(colour = DnaType)) +
  facet_wrap(.~ab.names)

sum$sample.type.year <- factor(sum$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheicwater", 
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver","RO3", "RO2", "RO1","Deep",
                                                                      "Downriver",
                                                                      "Marine"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Hyporheicwater", 
                                             "Groundwater","Stream", "Tributary",
                                             "Headwater \nLakes", "Upstream \nPonds", "Lake", "Lake",
                                             "Upriver","RO3","RO2", "RO1","Hypolimnion","Downriver",
                                             "Estuary"))

ggplot(means[ab.names == "rare.shifter",], 
       aes(x = sample.type.year, y = ASV)) +
  geom_raster(aes(fill = log2(mean.css + 1))) +
  scale_fill_viridis_c() +
  facet_wrap(.~DnaType)




# put names of each list bin into list itself for plotting
for(i in 1:length(ls.asvs)){
  ls.asvs[[i]] <- as.list(ls.asvs[i])
  names(ls.asvs[[i]]) <- "ASVs"
  ls.asvs[[i]]$name <- names(ls.asvs)[i]
}

# Run multivariate analysis
ls.pcoa <- llply(ls.asvs, function(x){
  if(length(x$ASVs) == 0L){
    warning("Abundance group is empty")
  } else {
    # extract species table with species in columns
    pb.mat <- t(otu_mat(pb))
    
    # take only ASVs part of one abundance group
    pb.mat <- pb.mat[, c(x$ASVs)]
    # remove samples (= rows) that have all 0 values
    pb.mat <- pb.mat[apply(pb.mat, 1, function(x) !all(x == 0)),]
    
    # 1. PCoA with Bray-Curtis
    pb.bray <- vegdist(pb.mat, method = "bray")
    pb.bray <- sqrt(pb.bray) # make Euclidean
    pb.bray.pcoa <- ape::pcoa(pb.bray)
    
    # 2. PCoA with Jaccard (presence-absence)
    pb.jac <- vegdist(pb.mat, method = "jaccard", binary = T)
    pb.jac.pcoa <- ape::pcoa(pb.jac)
    
    # merge outputs into a list
    list(asv.tab = pb.mat,
         pcoa.bray = pb.bray.pcoa,
         pcoa.jacc = pb.jac.pcoa,
         name = x$name)
  }
}, .parallel = T)

# Plot PCoAs
pcoa.df <- llply(ls.pcoa, function(x){
  out <- plot_bray_n_jacc(bray = x$pcoa.bray,
                   jacc = x$pcoa.jacc,
                   colours = colvec,
                   plot.name = x$name,
                   output = T)
  out$name <- x$name # add group name to list
  return(out)
}, .parallel = T)

# plot DNA-RNA distance
dnarna.dist <- llply(pcoa.df, function(x){
  dist.dr <- dist.dnarna(x$df, save.name = x$name, dimensions = 3)
}, .parallel = T)

# plot DNA-RNA dissimilarity
dnarna.abgroup <- llply(ls.asvs, function(x){
  physeq <- subset_asv(pb, x$ASVs)
  
  dist.dr <- dissim.dnarna(physeq, save.name = x$name, output = T)
  
}, .parallel = T)

set.seed(3) #plot heatmap works with NMDS
# plot heatmaps
heatmap.ab <- llply(ls.asvs, function(x){
  physeq <- subset_asv(pb, x$ASVs)
  
  if(any(taxa_sums(physeq) == 0) == T){
    physeq <- prune_taxa(!taxa_sums(physeq) == 0, physeq)
  }
  if(any(sample_sums(physeq) == 0) == T){
    physeq <- prune_samples(!sample_sums(physeq) == 0, physeq)
  }
  
  if(length(x$ASVs) > 1000){
    physeq <- prune_taxa(names(sort(taxa_sums(physeq), TRUE)[1:5000]), physeq)
  }
  
  p <- plot_heatmap(physeq)
  return(p)
}, .parallel = T)



out.groups <- ldply(dnarna.abgroup, function(x){
  x$wide
}, .parallel = T)


out.all <- dist.dr$wide %>% mutate(.id = "All")

out <- rbind(out.all, out.groups)

ggplot(out, aes(x = Bray, y = Jaccard)) +
  geom_point(aes(fill = sample.type.year), shape = 21) +
  geom_smooth(aes(group = sample.type.year, colour = sample.type.year), method = "lm", se = F) +
  annotate(geom = "segment", x = 0, y = 0, xend = 1, yend = 1)

# Bray Curtis
# 1 = don't share any species
# 2 = share all species








##########################################################################
#-----------------------#
# Correlate to richness #
#-----------------------#
# read in rarefied datasets
min_lib <- c(15000, 25000, 50000)

perm.rar <- select_newest("./Output",
              "perm.rar_lib", by = min_lib)

perm.rar <- c(perm.rar, select_newest("./Output", "201520162017_CSS_asvtab"))

# register number of cores for parallel computing
detectCores() #24, we do not have 24. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12

alpha <- llply(as.list(perm.rar), function(x){
  if(grepl(x, pattern = "asvtab") == T){
    rar <- read.csv(paste0("./Output/", x), sep = "\t", dec = ".", stringsAsFactors = F)
    colnames(rar)[1] <- "ASV"
    rar <- melt.data.table(setDT(rar),
                    id.vars = "ASV", # skip measure.var, takes all columns
                    variable.name = "Sample",
                    value.name = "iter.mean")
  } else {
    rar <- read.csv(paste0("./Output/", x), sep = ";", stringsAsFactors = F)
  }
  data.table::setDT(rar)
  # make mean column to count
  rar[, iter.mean := round(iter.mean, 0)]
  rar <- setDF(dcast(rar, Sample ~ ASV, value.var = "iter.mean"))
  
  # shannon-wiener index (H')
  # quantifies uncertainty associated with predicted identity of a new taxa
  # given number of taxa and evenness in abundances of individuals within each taxa
  # assumes the sample for site was collected randomly
  # more sensitive to rare species, can be positive and negative, typically range from 1.5 to 3.5
  
  # simpson's index (λ)
  # measure of dominance and as such weights towards the abundance of the most common taxa
  # higher values represent higher diversity, ranges from 0 to 1
  
  # pielou's evenness (J)
  # compares actual diversity value (e.g. Shannon-Wiener) to the maximum possible diversity value
  # when all species are equally common it is considered as the highest degree of evenness
  # ranges from 0 and 1
  # the more variation in abundances between different taxa within the community the lower is J
  # highly dependent on sample size and highly sensitive to rare taxa
  
  # Chao 1 richness esimator
  alpha <- plyr::ddply(rar, ~ Sample, function(z) {
    data.frame(Shannon = vegan::diversity(z[,-1], index = "shannon"),
               Simpson = vegan::diversity(z[,-1], index = "simpson"),
               Pielou = vegan::diversity(z[,-1], index = "shannon") / log(sum(z[,-1] > 0)),
               Chao1 = vegan::estimateR(z[,-1],)["S.chao1",])
  })
  
  return(alpha)
}, .parallel = T)

names(alpha) <- c(paste0("lib", min_lib), "css")
alpha.df <- bind_rows(alpha, .id = "Data")

write.table(alpha.df, "./Output/alpha_div_summary.csv", sep = ";", dec = ".", row.names = F)
alpha.df <- read.csv("./Output/alpha_div_summary.csv", sep = ";", dec = ".", stringsAsFactors = F)

#------------------------------------------------------------------------------------------------#
# Combine with bray distance data
# correct a few wrong sample names for matching DNA and RNA
alpha.df[alpha.df$Sample == "RO2R52R", "Sample"] <- "RO2.52R"
alpha.df[alpha.df$Sample == "SWR34R", "Sample"] <- "SW34R"
alpha.df[alpha.df$Sample == "RO2.36pD", "Sample"] <- "RO2.36D"
alpha.df[alpha.df$Sample == "RO2.36pR", "Sample"] <- "RO2.36R"
alpha.df[alpha.df$Sample == "RO2111.60mD", "Sample"] <- "RO2111.90mD"
alpha.df[alpha.df$Sample == "RO2.30DPR", "Sample"] <- "RO2.30R" # two DNA
alpha.df[alpha.df$Sample == "RO301.HypoR", "Sample"] <- "RO31.HypoR"
alpha.df[alpha.df$Sample == "RO301R", "Sample"] <- "RO31R" 
alpha.df[alpha.df$Sample == "RO304R", "Sample"] <- "RO34R" 
alpha.df[alpha.df$Sample == "RO307R", "Sample"] <- "RO37R" 
alpha.df[alpha.df$Sample == "L230R", "Sample"] <- "L330R" # L230 does not exist

# get ID
setDT(alpha.df); setDT(pdataframe)
alpha.df[pdataframe[Metric == "Bray"], c("ID", "DnaType") := 
           list(i.ID, i.DnaType), on = .(Sample)]

# split DNA and RNA
dna.alpha <- alpha.df[alpha.df$DnaType == "DNA",]
rna.alpha <- alpha.df[alpha.df$DnaType == "RNA",]

############################
# Explore DNA and RNA relations
# calculate mean coordinates for duplicates
sum <- alpha.df %>%
  dplyr::group_by(ID, Data, DnaType) %>%
  dplyr::summarise(Shannon = mean(Shannon, na.rm = T),
                   Simpson = mean(Simpson, na.rm = T),
                   Pielou = mean(Pielou, na.rm = T),
                   Chao1 = mean(Chao1, na.rm = T)) %>%
  ungroup()

setDT(sum)
temp <- dcast(sum, ID + Data ~ DnaType, value.var = c("Simpson","Shannon","Pielou","Chao1"))

temp <- melt(temp, id.vars = c("ID","Data"),
     variable.name = "Index",
     value.name = "Diversity") %>%
  separate(Index, into = c("Index","DnaType"), sep = "_")

temp <- dcast(temp, ID + Data + Index ~ DnaType, value.var = "Diversity") 


ggplot(temp[Index == "Pielou",], aes(x = DNA, y = RNA)) +
  geom_point()

# Correlate richness to distance in PCoA space between DNA and RNA

dna.alpha <- melt(setDT(dna.alpha), id.vars = c("ID","Data"),
                  measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
                  variable.name = "Index",
                  value.name = "Diversity")

# merge with distance
plot.df <- dna.alpha[dist.dr[Metric == "Bray"], c("distance.bray",
                                                  "sample.type.year") := 
                       list(i.distance, i.sample.type.year), on = .(ID)]

#plot.df <- dna.alpha[dist.dr[Metric == "Jaccard"], c("distance.jac",
#                                                  "sample.type.year") := 
#                       list(i.distance, i.sample.type.year), on = .(ID)]

# create colour vector for plotting
colvec <- c("red4","chocolate3","orangered2","orange3",
            "cadetblue", "darksalmon",
            "darkolivegreen","darkolivegreen3",
            "royalblue","mediumorchid4", "violet","palevioletred2","navy","skyblue",
            "seagreen3")

plot.df$Data <- factor(plot.df$Data, levels = c("css", "lib15000", "lib25000", "lib50000"),
                       labels = c("CSS", "Rarefied: Lib15000", "Rarefied: Lib25000","Rarefied: Lib50000"))

# make plot with ggpubr to include pearson's correlation outputs directly in the plot
p <- ggscatter(plot.df[!is.na(distance.bray),], x = "distance.bray", y = "Diversity", 
               color = "sample.type.year",
               palette = colvec,
               xlab = "Distance between DNA and RNA \nin PCoA space (Bray Curtis)",
               ylab = "DNA Alpha Diversity",
               legend.title = "Sample Type")
(pf <- facet(p, facet.by = c("Index","Data"), ncol = 3, nrow = 4, scales = "free")+
  stat_cor(method = "pearson", label.x = 0.1, cor.coef.name = "r"))

ggsave("./Figures/Final/All_DNARNA_distance_alphadiv_comparison.png", pf,
       width = 30, height = 18, unit = "cm")


# There is no substantial effect of data transformation on alpha diversity results (CSS vs Rarefaction)
# We continue with CSS to keep consistent with the underlying data structure

# We focus on two alpha diversity indices:
# Shannon-Wiener and Pielou

alpha.df <- plot.df[Data == "CSS" & (Index == "Shannon" | Index == "Pielou"),]
alpha.df <- alpha.df[!is.na(distance.bray),]
setorderv(alpha.df, c("Index","Diversity")) # rearrange dataframe

z <- alpha.df[Index == "Shannon"]
z <- z %>%
  dplyr::group_by(sample.type.year) %>%
  dplyr::summarise(Diversity = mean(Diversity, na.rm = T),
                   distance.bray = mean(distance.bray, na.rm = T))
plot(Diversity ~ distance.bray, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)
lm1 <- lm(z$Diversity ~ z$distance.bray)
lm2 <- lm(z$Diversity ~ poly(z$distance.bray, 2))
lm3 <- lm(z$Diversity ~ poly(z$distance.bray, 3))
anova(lm0,lm1) # preferred model is lm1
anova(lm1,lm2) 
anova(lm2,lm3) 
rm(z)


z <- alpha.df[Index == "Pielou"]
z <- z %>%
  dplyr::group_by(sample.type.year) %>%
  dplyr::summarise(Diversity = mean(Diversity, na.rm = T),
                   distance.bray = mean(distance.bray, na.rm = T))
plot(Diversity ~ distance.bray, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)
lm1 <- lm(z$Diversity ~ z$distance.bray)
lm2 <- lm(z$Diversity ~ poly(z$distance.bray, 2))
lm3 <- lm(z$Diversity ~ poly(z$distance.bray, 3))
anova(lm0,lm1)
anova(lm1,lm2) # preferred model is lm2
anova(lm2,lm3) 

# lm2 equation:
# y = beta[0] + beta[1]x + beta[2]x^2 + epsilon
rm(z)

# both models are best with a polynomial degree 2

lin.ls <- dlply(alpha.df, .(Index), function(z){
  
  means <- z %>%
    dplyr::group_by(sample.type.year) %>%
    dplyr::summarise(Diversity = mean(Diversity, na.rm = T),
                     distance.bray = mean(distance.bray, na.rm = T))
  if(unique(z$Index) == "Shannon"){
    lin <- lm(means$Diversity ~ means$distance.bray)
  } else if(unique(z$Index) == "Pielou") {
    lin <-  lm(means$Diversity ~ poly(means$distance.bray, 2))
  }

  # check linear assumptions
  #plot(lin) # normality not good
  # large sample sizes, normality does not affect results too much (central limit theorem)
  # homoscedasticity and independence important
  #summary(lin)
  #confint(lin, level = 0.95)
  
  # get data for plotting
  x <- data.frame(x = sort(means$distance.bray))
  pred <- predict(lin, newdata = x, se = T)
  ci <- pred$se.fit[order(means$distance.bray)] * qt(0.95 / 2 + 0.5, pred$df)
  y <- pred$fit[order(means$distance.bray)]
  ymin <- y - ci
  ymax <- y + ci
  
  plot.df <- data.frame(x = sort(means$distance.bray),
                        y = y,
                        ymin = ymin,
                        ymax = ymax,
                        se = pred$se.fit[order(means$distance.bray)])
  
  colvec <- c("red4","chocolate3","orangered2","orange3",
              "cadetblue", "darksalmon",
              "darkolivegreen","darkolivegreen3",
              "royalblue","mediumorchid4", "violet","palevioletred2","navy","skyblue",
              "seagreen3")
  
  p <- ggplot() +
    theme_pubr() +
    geom_point(data = z, aes(x = distance.bray, y = Diversity), 
               colour = "gray40", alpha = 0.3, size = 2) +
    geom_line(data = plot.df, aes(x = x, y = y), inherit.aes = F, size = 1) +
    geom_line(data = plot.df, aes(x = x, y = ymax), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_line(data = plot.df, aes(x = x, y = ymin), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_point(data = means,
               aes(x = distance.bray, y = Diversity, fill = sample.type.year), shape = 21, size = 3) +
    scale_fill_manual(values = colvec, name = "Sample Type") +
    labs(x = "",
         y = paste0(unique(z$Index)))

  # get model statistics
  options(scipen = 999) # avoid scientific annotations
  
  fnr <- substitute(italic(R)^2~"="~r2*","~~italic(F)[df]~"="~Fstat,
                    list(r2 = format(summary(lin)$r.squared, digits = 2),
                         Fstat = format(summary(lin)$fstatistic[[1]], digits = 4),
                         df = paste0(format(summary(lin)$fstatistic[[2]], digits = 0),
                                     ",", format(summary(lin)$fstatistic[[3]], digits = 0))))
  pv1 <- summary(lin)$coefficients[2,4]
  pv1 <- if(pv1 < 0.0001){
    "< 0.0001"} else if(pv1 < 0.001){
      "< 0.001"} else if(pv1 < 0.01){
        "< 0.01"} else if(pv1 < 0.05){
          "< 0.05"
        } else {
          paste("=",round(pv1, 2))
        }
  
  if(unique(z$Index) == "Pielou"){
    eq1 <- substitute(italic(y) == a - b %.% italic(x) + b2 %.% italic(x)^2,
                      list(a = format(as.vector(coef(lin)[1]), digits = 2),
                           b = format(as.vector(abs(coef(lin)[2])), digits = 2),
                           b2 = format(as.vector(coef(lin)[3]), digits = 2)))
    
    pv2 <- summary(lin)$coefficients[3,4]
    pv2 <- if(pv2 < 0.0001){
      "< 0.0001"} else if(pv2 < 0.001){
        "< 0.001"} else if(pv2 < 0.01){
          "< 0.01"} else if(pv2 < 0.05){
            "< 0.05"
          } else {
            paste("=",round(pv2, 2))
          }
    ps <- substitute(italic(p)[beta[1]]~pval1*","~italic(p)[beta[2]]~pval2,
                     list(pval1 = pv1,
                          pval2 = pv2))
    (p <- p + 
        annotate("text", x = 0.4, y = 0.9, 
                 label = as.character(as.expression(eq1)), parse = T, size = 2.5)+
        annotate("text", x = 0.4, y = 0.88, 
                 label = as.character(as.expression(fnr)), parse = T, size = 2.5) +
        annotate("text", x = 0.4, y = 0.86, 
                 label = as.character(as.expression(ps)), parse = T, size = 2.5) 
    )
  } else {
    eq1 <- substitute(italic(y) == a - b %.% italic(x),
                      list(a = format(as.vector(coef(lin)[1]), digits = 2),
                           b = format(as.vector(abs(coef(lin)[2])), digits = 2)))
    ps <- substitute(italic(p)~pval1,
                     list(pval1 = pv1))
    
    (p <- p + 
        annotate("text", x = 0.4, y = 6.5, 
                 label = as.character(as.expression(eq1)), parse = T, size = 2.5)+
        annotate("text", x = 0.4, y = 6.3, 
                 label = as.character(as.expression(fnr)), parse = T, size = 2.5) +
        annotate("text", x = 0.4, y = 6.1, 
                 label = as.character(as.expression(ps)), parse = T, size = 2.5) 
    )
  }
  
  #, abs(round(coef(lin)[2], 2)), "*x +",
  #round(coef(lin)[3], 2), "*x"^2*""
  
  list(original = z,
       binned = means,
       lin = lin,
       coef = coef(lin),
       fitted = lin.df,
       plot = p)

})

(p <- ggarrange(lin.ls[[1]]$plot, lin.ls[[2]]$plot,
                ncol = 2, common.legend = T, legend = "right"))

(p <- annotate_figure(p, 
                      bottom = text_grob("Distance between DNA and RNA in PCoA space (Bray Curtis)",
                                         just = "centre")))
ggsave("./Figures/Final/Richness_distance_nonlin_reg.png", p,
       width = 22, height = 11, unit = "cm")

ggplot() +
  theme_pubr() +
  facet_grid(.~Index) + 
  geom_point(data = alpha.df, aes(x = distance.bray, y = Diversity), colour = "gray40", alpha = 0.5) +
  geom_line(data = lin.df, aes(x = x, y = y), inherit.aes = F, size = 1) +
  geom_line(data = lin.df, aes(x = x, y = ymax), inherit.aes = F, size = 0.7, linetype = "dashed") +
  geom_line(data = lin.df, aes(x = x, y = ymin), inherit.aes = F, size = 0.7, linetype = "dashed") +
  geom_point(data = means,
             aes(x = mean.dist, y = mean.div, fill = sample.type.year), shape = 21, size = 3) +
  labs(x = "Distance between DNA and RNA \nin PCoA space (Bray Curtis)", y = "Shannon-Wiener Index")

ggplot(shan, aes(x = distance.bray, y = Diversity)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x,3)) +
  stat_smooth(aes(outfit = fit <<- ..y..), method = "lm", formula = y ~ poly(x,3))

lin.p <- lm(exp(pie$Diversity) ~ poly(pie$distance.bray,2))
plot(lin.p) # normality not good


ggplot(alpha.df, aes(x = distance.bray, y = Diversity)) +
  geom_point() +
  facet_grid(.~Index, scales = "free")


# Taxonomic composition
# We want to show the taxonomic composition of our samples plus the abundance
# As we have too many samples, best would probably be to calcualte the mean abundance for each group
# Groups are: sample.type.year + DnaType + Season

# melt ASV table
pb.df <- as.data.frame(otu_mat(pb)) # make data.frame
setDT(pb.df, keep.rownames = "ASV") # make data.table

# melt
pb.df <- melt.data.table(pb.df, 
                id.vars = "ASV", # skip measure.var, takes all columns
                variable.name = "Sample",
                value.name = "css")

# add taxonomy data to the abundance data
tax.df <- as.data.frame(tax.tab) # make data.frame
setDT(tax.df, keep.rownames = "ASV") # make data.table
pb.df <- pb.df[tax.df, on = .(ASV)] # merge

# add meta data
meta <- sample_df(met.df) %>%
  dplyr::select(DnaType, Year, Season, sample.type.year, soilorwater, LibrarySize)
setDT(meta, keep.rownames = "Sample") # make data.table

pb.df <- pb.df[meta, on = .(Sample)] # merge

# calculate mean abundance (css) for each category
mean.pb <- pb.df[, .(css.mean = mean(css, na.rm = T),
                     css.sd = sd(css, na.rm = T)),
          by = .(ASV, sample.type.year, Season, DnaType)][
            , css.sum := sum(css.mean, na.rm = T),
            by = .(sample.type.year, Season, DnaType)
          ][, css.rel := css.mean * 1 / css.sum]
mean.pb <- mean.pb[tax.df, on = .(ASV)] # merge with taxonomy data
mean.pb[, ID := paste(DnaType, Season, sample.type.year, sep = "_")] # add plot ID

# calculate mean library size per category
lib.size <- pb.df[, .(lib.mean = mean(LibrarySize, na.rm = T),
                      lib.sd = sd(LibrarySize, na.rm = T)),
                  by = .(sample.type.year, Season, DnaType)]
lib.size[, ID := paste(DnaType, Season, sample.type.year, sep = "_")] # add plot ID

# Overwrite IDs as factors to set a order for plotting
mean.pb$ID <- factor(mean.pb$ID, 
                                      levels = c("DNA_spring_Soil", "DNA_spring_Sediment",
                                                 "DNA_spring_Soilwater","DNA_spring_Hyporheicwater", 
                                                 "DNA_spring_Wellwater", "DNA_spring_Stream",
                                                 "DNA_spring_Tributary", "DNA_spring_HeadwaterLakes",
                                                 "DNA_spring_PRLake", "DNA_spring_Lake",
                                                 "DNA_spring_IslandLake", "DNA_spring_Upriver",
                                                 "DNA_spring_RO2", 
                                                 "DNA_spring_RO1",
                                                 "DNA_spring_Downriver", "DNA_spring_Marine",
                                                 "DNA_summer_Soil", "DNA_summer_Sediment",
                                                 "DNA_summer_Soilwater","DNA_summer_Hyporheicwater", 
                                                 "DNA_summer_Stream",
                                                 "DNA_summer_Tributary", "DNA_summer_HeadwaterLakes",
                                                 "DNA_summer_PRLake", "DNA_summer_Lake",
                                                 "DNA_summer_IslandLake", "DNA_summer_Upriver",
                                                 "DNA_summer_RO3", "DNA_summer_RO2", 
                                                 "DNA_summer_RO1", "DNA_summer_Deep",
                                                 "DNA_summer_Downriver", "DNA_summer_Marine",
                                                 "DNA_autumn_Tributary", "DNA_autumn_HeadwaterLakes",
                                                 "DNA_autumn_Lake", "DNA_autumn_Upriver",
                                                 "DNA_autumn_RO3", "DNA_autumn_RO2", 
                                                 "DNA_autumn_RO1", "DNA_autumn_Deep",
                                                 "DNA_autumn_Downriver",
                                                 "cDNA_spring_Soil",
                                                 "cDNA_spring_Soilwater","cDNA_spring_Hyporheicwater", 
                                                 "cDNA_spring_Stream",
                                                 "cDNA_spring_Tributary", "cDNA_spring_HeadwaterLakes",
                                                 "cDNA_spring_Lake",
                                                 "cDNA_spring_Upriver",
                                                 "cDNA_spring_RO2", 
                                                 "cDNA_spring_RO1",
                                                 "cDNA_spring_Downriver",
                                                 "cDNA_summer_Soil", "cDNA_summer_Sediment",
                                                 "cDNA_summer_Soilwater","cDNA_summer_Hyporheicwater", 
                                                 "cDNA_summer_Stream",
                                                 "cDNA_summer_Tributary", "cDNA_summer_HeadwaterLakes",
                                                 "cDNA_summer_Lake",
                                                 "cDNA_summer_Upriver",
                                                 "cDNA_summer_RO3", "cDNA_summer_RO2", 
                                                 "cDNA_summer_RO1", "cDNA_summer_Deep",
                                                 "cDNA_summer_Downriver", "cDNA_summer_Marine",
                                                 "cDNA_autumn_Tributary", "cDNA_autumn_HeadwaterLakes",
                                                 "cDNA_autumn_Lake", "cDNA_autumn_Upriver",
                                                 "cDNA_autumn_RO3", "cDNA_autumn_RO2", 
                                                 "cDNA_autumn_RO1", "cDNA_autumn_Deep",
                                                 "cDNA_autumn_Downriver"))

lib.size$ID <- factor(lib.size$ID, 
                     levels = c("DNA_spring_Soil", "DNA_spring_Sediment",
                                "DNA_spring_Soilwater","DNA_spring_Hyporheicwater", 
                                "DNA_spring_Wellwater", "DNA_spring_Stream",
                                "DNA_spring_Tributary", "DNA_spring_HeadwaterLakes",
                                "DNA_spring_PRLake", "DNA_spring_Lake",
                                "DNA_spring_IslandLake", "DNA_spring_Upriver",
                                "DNA_spring_RO2", 
                                "DNA_spring_RO1",
                                "DNA_spring_Downriver", "DNA_spring_Marine",
                                "DNA_summer_Soil", "DNA_summer_Sediment",
                                "DNA_summer_Soilwater","DNA_summer_Hyporheicwater", 
                                "DNA_summer_Stream",
                                "DNA_summer_Tributary", "DNA_summer_HeadwaterLakes",
                                "DNA_summer_PRLake", "DNA_summer_Lake",
                                "DNA_summer_IslandLake", "DNA_summer_Upriver",
                                "DNA_summer_RO3", "DNA_summer_RO2", 
                                "DNA_summer_RO1", "DNA_summer_Deep",
                                "DNA_summer_Downriver", "DNA_summer_Marine",
                                "DNA_autumn_Tributary", "DNA_autumn_HeadwaterLakes",
                                "DNA_autumn_Lake", "DNA_autumn_Upriver",
                                "DNA_autumn_RO3", "DNA_autumn_RO2", 
                                "DNA_autumn_RO1", "DNA_autumn_Deep",
                                "DNA_autumn_Downriver",
                                "cDNA_spring_Soil",
                                "cDNA_spring_Soilwater","cDNA_spring_Hyporheicwater", 
                                "cDNA_spring_Stream",
                                "cDNA_spring_Tributary", "cDNA_spring_HeadwaterLakes",
                                "cDNA_spring_Lake",
                                "cDNA_spring_Upriver",
                                "cDNA_spring_RO2", 
                                "cDNA_spring_RO1",
                                "cDNA_spring_Downriver",
                                "cDNA_summer_Soil", "cDNA_summer_Sediment",
                                "cDNA_summer_Soilwater","cDNA_summer_Hyporheicwater", 
                                "cDNA_summer_Stream",
                                "cDNA_summer_Tributary", "cDNA_summer_HeadwaterLakes",
                                "cDNA_summer_Lake",
                                "cDNA_summer_Upriver",
                                "cDNA_summer_RO3", "cDNA_summer_RO2", 
                                "cDNA_summer_RO1", "cDNA_summer_Deep",
                                "cDNA_summer_Downriver", "cDNA_summer_Marine",
                                "cDNA_autumn_Tributary", "cDNA_autumn_HeadwaterLakes",
                                "cDNA_autumn_Lake", "cDNA_autumn_Upriver",
                                "cDNA_autumn_RO3", "cDNA_autumn_RO2", 
                                "cDNA_autumn_RO1", "cDNA_autumn_Deep",
                                "cDNA_autumn_Downriver"))


# colour blind friendly, derived from https://medialab.github.io/iwanthue/
col_vector<-c("#dcd873","#350070","#b9ce40","#003499","#f5d249","#6b8aff","#06a644","#ec83f6","#7beb8b",
"#c947b1","#01b072","#df3587","#006f26","#ff83d0","#215a00","#99a1ff","#668200",
"#015db8","#e77e28","#019cf8","#b5221d","#67b8ff","#bd0a35","#b2e294","#840066",
"#314800","#ffb0ed","#954700","#3d0e52","#ff9b61","#59003b","#ff6e83","#aa7dbf","#620009")

# normal
#col_vector<-c("#6893ff","#bbce1a","#6b63ed","#18cc58","#b22bb2",
#"#71dd6e","#ff6fea","#428700","#005eca","#f1bf45","#36508f","#d78600","#4bb8ff",
#"#ce5300","#5bd5f6","#da0053","#01823f","#c80074","#97d68c","#a5006f","#bccf63",
#"#85307a","#848a00","#ff94d5","#006845","#ff6f91","#00b5be","#a8041e","#9ac89a",
#"#ff7f64","#707c48","#ff9b3a","#635b00","#9c5700")

x.types.labs<-c("Soil","Sediment","Soilwater","Hyporheicwater","Groundwater","Stream",
"Tributary","Headwater Lakes","Upstream Ponds","Lake",
"IslandLake","Upriver","RO2","RO1","Downriver","Marine",
"Soil","Sediment","Soilwater","Hyporheicwater","Stream","Tributary","Headwater Lakes",
"Upstream Ponds","Lake","IslandLake","Upriver","RO3","RO2","RO1","Hypolimnion","Downriver","Marine",
"Tributary","Headwater Lakes","Lake","Upriver","RO3","RO2","RO1","Hypolimnion","Downriver",
"Soil","Soilwater","Hyporheicwater","Stream","Tributary","Headwater Lakes",
"Lake","Upriver","RO2","RO1","Downriver",
"Soil","Sediment","Soilwater","Hyporheicwater",
"Stream","Tributary","Headwater Lakes","Lake","Upriver",
"RO3","RO2","RO1","Hypolimnion","Downriver","Marine",
"Tributary","Headwater Lakes","Lake","Upriver",
"RO3","RO2","RO1","Hypolimnion","Downriver")

lib.bar <- 
    ggplot(lib.size, aes(x = ID, y = lib.mean / 1000)) +
    theme_pubr() +
    geom_bar(stat = "identity", width = 0.7, fill = "white", colour = "grey40") +
    geom_errorbar(aes(ymin = (lib.mean - lib.sd) / 1000, 
                      ymax = (lib.mean + lib.sd) / 1000), width = 0.4) + 
    labs(y = expression(paste("Average library size [x10"^3,"]"))) +
    scale_x_discrete(labels = x.types.labs, expand = c(0.025,0.025)) +
    theme(axis.text.x = element_blank(),
         axis.ticks = element_blank(), axis.title.x = element_blank(),
         axis.title = element_text(size = 9)) +
  coord_cartesian(clip="off")

  tax.leg <-
    get_legend(
      ggplot(mean.pb, 
             aes(x = ID, y = css.rel, fill = phylum, colour = phylum)) +
        theme_bw() +
        geom_bar(stat = "identity", width = 0.7) +
        labs(y = "Average relative abundance") +
        scale_fill_manual(name = "Phylum", values = rev(col_vector)) +
        scale_colour_manual(name = "Phylum", values = rev(col_vector)) +
        theme(legend.position = "right")
    )


tax <- ggplot(mean.pb, 
       aes(x = ID, y = css.rel, fill = phylum, colour = phylum)) +
  theme_pubr() +
  geom_bar(stat = "identity", width = 0.8) +
  guides(fill = "none", colour = "none") +
  labs(y = "Average relative abundance") +
  scale_fill_manual(values = rev(col_vector)) +
    scale_colour_manual(values = rev(col_vector)) +
  scale_x_discrete(labels = x.types.labs, expand = c(0.025,0.025)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(2,2,15,2), "mm"))

# adding axis labels to the bottom of x axis
labelled.tax <- tax + 
  coord_cartesian(ylim=c(0,1), clip="off") +
  annotate("segment", x = 0, xend = 16.5, y = -0.4, yend = -0.4) +
  annotate("segment", x = 16.7, xend = 33.5, y = -0.4, yend = -0.4) +
  annotate("segment", x = 33.7, xend = 42.5, y = -0.4, yend = -0.4) +
  annotate("segment", x = 43, xend = 53.5, y = -0.4, yend = -0.4) +
  annotate("segment", x = 53.7, xend = 68.5, y = -0.4, yend = -0.4) +
  annotate("segment", x = 68.7, xend = 78.5, y = -0.4, yend = -0.4) +
    annotate("text", x = 8.5, y = -0.425, label = "Spring") +
    annotate("text", x = 25.1, y = -0.425, label = "Summer") +
    annotate("text", x = 38.1, y = -0.425, label = "Autumn") +
    annotate("text", x = 48.25, y = -0.425, label = "Spring") +
    annotate("text", x = 61.1, y = -0.425, label = "Summer") +
    annotate("text", x = 73.635, y = -0.425, label = "Autumn") +
    annotate("segment", x = 0, xend = 42.5, y = -0.48, yend = -0.48) +
    annotate("segment", x = 43, xend = 78.5, y = -0.48, yend = -0.48) +
    annotate("text", x = 21.25, y = -0.52, label = "DNA") +
    annotate("text", x = 60.75, y = -0.52, label = "RNA")

combo  <- ggarrange(lib.bar, labelled.tax, nrow = 2,
                    align = "v", heights = c(0.2, 0.8), legend.grob = tax.leg, legend = "right")
          
#tax.leg, nrow = 2, heights = c(0.9, 0.1), widths = c(0.8, 0.2)
ggsave("./Figures/Final/Tax_LibSiz_Phyla.png", combo,
       width = 360, height = 203, unit = "mm", dpi = 300)


#-----------#
# 2015-2016 #
#-----------#
#################################################################################################################
#################################################################################################################
#-----------#
# Read data #
#-----------#

# do we have several files per object? -> take newest version
# ASV CSS transformed table
asv.tab <- select_newest("./Output", "20152016_CSS_asvtab")
asv.tab <- read.csv(
  paste0("./Output/", asv.tab),
  sep = "\t",
  dec = ".",
  stringsAsFactors = F
)

# transpose back to ASV in cols, samples in rows
row.names(asv.tab) <- asv.tab$Taxa.and.Samples
asv.tab[, "Taxa.and.Samples"] <- NULL
asv.tab <- as.matrix(asv.tab)
# row orders need to match between tax.tab and asv.tab
asv.tab <- asv.tab[order(row.names(asv.tab)),]

# Taxonomy table
tax.tab <- select_newest("./Output", "20152016_tax_table")
tax.tab <-
  as.matrix(read.csv(
    paste0("./Output/", tax.tab),
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ))
# orders need to match between tax.tab and asv.tab
tax.tab <- tax.tab[order(row.names(tax.tab)),]

# Meta data
met.df <-
  select_newest(path = "./Output", file.pattern = "20152016_meta_data")
met.df <-
  read.csv(
    paste0("./Output/", met.df),
    sep = ";",
    dec = ".",
    stringsAsFactors = F
  )
# phyloseq needs the sample names of the meta data to be the same as the microbial data
met.df <- sample_data(met.df)

# Assign rownames to be Sample ID's
rownames(met.df) <- met.df$DadaNames

# Construct phyloseq object
pb <- phyloseq(otu_table(asv.tab, taxa_are_rows = T),
               sample_data(met.df),
               tax_table(tax.tab))

############
# Only DNA #
############

# filter only DNA
dna <- pb %>%
  subset_samples(
    DnaType == "DNA"
  )

############
# Analysis #
############
# make ordinations
# tried NMDS stress does not reach convergence
dist <- "bray"
ord_meths <- "PCoA"

####################
# Both DNA and RNA #
####################

plist <- llply(as.list(ord_meths), function(i, physeq, dist) {
  ordi = ordinate(physeq = physeq,
                  method = i ,
                  distance = dist)
  plot_ordination(physeq = physeq, ordi)
}, pb, dist)

names(plist) <- ord_meths

pdataframe <- ldply(plist, function(x){
  df <- x$data[, 1:2]
  return(cbind(df, x$data[,-c(1:2)]))
})

var.exp <- ldply(plist, function(z) {
  data.frame(x = str_extract(z$labels$x, "\\[.+?\\]"),
             y = str_extract(z$labels$y, "\\[.+?\\]"))
})

pdataframe <- merge(pdataframe, var.exp, by = ".id")

############
# Plotting #
############
# put factor order and rename for plotting
# overwrite factors for plotting
pdataframe$sample.type.year <- factor(pdataframe$sample.type.year, levels = c("Soil","Sediment",
                                                                            "Soilwater","Hyporheicwater", 
                                                                            "Wellwater","Stream", "Tributary",
                                                                            "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                            "Upriver","RO2", "RO1","Downriver",
                                                                            "Marine", "Bioassay", "Blank"),
                                     labels = c("Soil","Sediment",
                                                "Soilwater","Hyporheicwater", 
                                                "Groundwater","Stream", "Tributary",
                                                "Headwater \nLakes", "Upstream \nPonds", "Lake", "Island \nLakes",
                                                "Upriver","RO2", "RO1","Downriver",
                                                "Estuary", "Bioassay", "Blank"))
pdataframe$Season <- factor(pdataframe$Season, levels = c("spring", "summer"), labels = c("Spring", "Summer"))
pdataframe$DnaType <- factor(pdataframe$DnaType, levels = c("DNA", "cDNA"), labels = c("DNA", "RNA"))

# create colour vector for plotting
colvec <- c("red4","chocolate3","orangered2","orange3",
            "cadetblue", "khaki", "darksalmon",
            "darkolivegreen","darkolivegreen3", "chartreuse3","yellow3",
            "royalblue","mediumorchid4", "violet","navy",
            "seagreen3","black", "beige")

(pcoa <- ggplot(pdataframe, aes(x = Axis.1, y = Axis.2)) +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
               size = 3) +
    #scale_fill_viridis_d(name = "Sample Type") +
    scale_fill_manual(values = colvec, name = "Sample Type") +
    scale_shape_manual(values = c(21,23,25)) +
    scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
    labs(x = paste("PCoA 1 ", pdataframe$x), y = paste("PCoA 2 ", pdataframe$y)) +
    theme(legend.key.size = unit(1.5, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(override.aes=list(shape=21))))


ggsave("./Figures/General/All_DNARNA_SampleType_PCoA.png", pcoa.arr,
       width = 23, height = 18, unit = "cm")

####################################################################################
# All plots show horseshoe effect (for 2015-2016)
# Try different methods to overcome this bias

# extract species table with species in columns
pb.mat <- t(otu_mat(pb))

# 1. Hellinger's transformation
pb.h <- decostand(pb.mat, "hellinger")
pb.h.pcoa <- ape::pcoa(vegdist(pb.h, method = "euc"))
is.euclid(dist(pb.h, method = "euclidean")) # TRUE
summary(pb.h.pcoa)
head(pb.h.pcoa$values) # eigenvalues = proportions explained by each axis
sum(pb.h.pcoa$values$Eigenvalues[1:2]) # 53.18 % explained with first two axes
biplot(pb.h.pcoa) #horseshoe

# 2. Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray.pcoa <- ape::pcoa(pb.bray)
head(pb.bray.pcoa$values)
sum(pb.bray.pcoa$values$Eigenvalues[1:2]) # 27.87 % explained with first two axes
biplot(pb.bray.pcoa) # horseshoe

# 3. Bray-Curtis with log+1
pb.log <- log1p(pb.mat)
pb.bray.log <- vegdist(pb.log, method = "bray")
is.euclid(pb.bray.log) # FALSE
pb.bray.log.pcoa <- ape::pcoa(pb.bray.log)
sum(pb.bray.log.pcoa$values$Eigenvalues[1:2]) # 26.57 % explained with first two axes
biplot(pb.bray.log.pcoa) # no horseshoe

# 4. Bray-Curtis with square root transformation
is.euclid(sqrt(pb.bray)) # TRUE
pb.braysq.pcoa <- ape::pcoa(sqrt(pb.bray))
head(pb.braysq.pcoa$values)
sum(pb.braysq.pcoa$values$Eigenvalues[1:2]) # 18.99 % explained with first two axes
biplot(pb.braysq.pcoa) # less strong horseshoe

# 5. Bray-Curtis with Lingoes (corrects negative eigenvalues)
pb.brayl.pcoa <- ape::pcoa(pb.bray, correction = "lingoes")
head(pb.brayl.pcoa$values)
sum(pb.brayl.pcoa$values$Eigenvalues[1:2]) # 27.87 % explained with first two axes
biplot(pb.brayl.pcoa) #horseshoe

# 6. Bray-Curtis with Cailliez (corrects negative eigenvalues)
pb.brayc.pcoa <- ape::pcoa(pb.bray, correction = "cailliez")
head(pb.brayc.pcoa$values)
sum(pb.brayc.pcoa$values$Eigenvalues[1:2]) # 27.87 % explained with first two axes
biplot(pb.brayc.pcoa) # horseshoe

### go with log1p for all ###
rm(pb.bray, pb.bray.pcoa, pb.brayl.pcoa, pb.brayc.pcoa, pb.braysq.pcoa, pb.h, pb.h.pcoa)

# Extract PCoA output based on log transformed data
# for all samples
pb.scores <- data.frame(Sample = row.names(pb.bray.log.pcoa$vectors),
                        pb.bray.log.pcoa$vectors[,1:2]) # get first two axes
pb.var <- pb.bray.log.pcoa$values$Eigenvalues[1:2]
meta <- data.frame(Sample = row.names(sample_df(pb)),
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, DnaType, LibrarySize,
                                                   bact.abundance, bact.production,catchment.area))
pb.scores <- merge(pb.scores, meta, by = "Sample")
pb.scores$Sample <- as.character(pb.scores$Sample)

# overwrite factors for plotting
pb.scores$sample.type.year <- factor(pb.scores$sample.type.year, levels = c("Soil","Sediment",
                                                                              "Soilwater","Hyporheicwater", 
                                                                              "Wellwater","Stream", "Tributary",
                                                                              "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                              "Upriver","RO2", "RO1","Downriver",
                                                                              "Marine", "Bioassay", "Blank"),
                                      labels = c("Soil","Sediment",
                                                 "Soilwater","Hyporheicwater", 
                                                 "Groundwater","Stream", "Tributary",
                                                 "Headwater \nLakes", "Upstream \nPonds", "Lake", "Island \nLakes",
                                                 "Upriver","RO2", "RO1","Downriver",
                                                 "Estuary", "Bioassay", "Blank"))
pb.scores$Season <- factor(pb.scores$Season, levels = c("spring", "summer"), labels = c("Spring", "Summer"))
pb.scores$DnaType <- factor(pb.scores$DnaType, levels = c("DNA", "cDNA"), labels = c("DNA", "RNA"))

# create colour vector for plotting
colvec <- c("red4","chocolate3","orangered2","orange3",
            "cadetblue", "khaki", "darksalmon",
            "darkolivegreen","darkolivegreen3", "chartreuse3","yellow3",
            "royalblue","mediumorchid4", "violet","navy",
            "seagreen3","black", "beige")

## Sample Type
(pcoa <- ggplot(pb.scores, aes(x = Axis.1, y = Axis.2)) +
    theme_bw()+
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
               size = 3) +
    #scale_fill_viridis_d(name = "Sample Type") +
    scale_fill_manual(values = colvec, name = "Sample Type") +
    scale_shape_manual(values = c(21,23)) +
    scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
    labs(x = paste0("PCoA 1 (", round(pb.var[1],2), " %)"),
         y = paste0("PCoA 2 (", round(pb.var[2], 2), " %)")) +
    guides(fill = guide_legend(override.aes=list(shape=21)))+
    lims(x = c(-0.45, 0.45), y = c(-0.35, 0.3)) +
    theme(legend.key.size = unit(1.5, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

ggsave("./Figures/General/2016_log_DNARNA_SampleType_PCoA.png", pcoa, width = 23, height = 20, unit = "cm")


#-----------------------------------------------#
# Calculate distance between DNA and RNA points #
#-----------------------------------------------#
# correct a few wrong sample names
pb.scores[pb.scores$Sample == "RO2R52R", "Sample"] <- "RO2.52R"
pb.scores[pb.scores$Sample == "SWR34R", "Sample"] <- "SW34R"
pb.scores[pb.scores$Sample == "RO2.36pD", "Sample"] <- "RO2.36D"
pb.scores[pb.scores$Sample == "RO2.36pR", "Sample"] <- "RO2.36R"

pb.scores$ID[pb.scores$DnaType == "DNA"] <- str_replace(pb.scores$Sample[pb.scores$DnaType == "DNA"], "D$", "")
pb.scores$ID[pb.scores$DnaType == "RNA"] <- str_replace(pb.scores$Sample[pb.scores$DnaType == "RNA"], "R$", "")

# calculate mean coordinates for duplicates
sum <- pb.scores %>% 
  filter(Year == 2016) %>% 
  dplyr::group_by(ID, DnaType) %>%
  dplyr::summarise(x = mean(Axis.1), y = mean(Axis.2),
                   sample.type.year = unique(sample.type.year),
                   Season = unique(Season),
                   Year = unique(Year), 
                   bact.abundance = unique(bact.abundance),
                   bact.production = unique(bact.production),
                   catchment.area = unique(catchment.area),
                   n = n()) %>%
  ungroup()

# exctract those that have both
dnarna <- sum[sum$ID %in% sum$ID[duplicated(sum$ID)],]

# calculate distance
setDT(dnarna)
temp <- dcast(dnarna, ID ~ DnaType, value.var = c("x","y"))
temp[, distance := sqrt((x_DNA - x_RNA)^2 + (y_DNA - y_RNA)^2)]

# combine back with categories
dist.dr <- merge(unique(dnarna, by = "ID"), temp, by = "ID")

# overwrite Headwater Lakes as Lake
dist.dr[sample.type.year == "Headwater \nLakes", sample.type.year := "Lake"]
dist.dr$sample.type.year <- factor(dist.dr$sample.type.year, levels = c("Soilwater","Hyporheicwater", 
                                                                        "Stream", "River",
                                                                        "Lake", 
                                                                        "Reservoir"))

(dr.dist <- ggplot(dist.dr, aes(x = sample.type.year, y = distance, fill = Season))+
    geom_boxplot(outlier.alpha = 0, width = 0.5) +
    geom_point(position = position_jitterdodge(), colour = "gray20", alpha = 0.5) +
    scale_fill_manual(values = c("#2D708EFF","#FDE725FF")) +
    labs(x = "Sample type", y = "Distance between DNA and RNA \nin PCoA space (log+1 Bray-Curtis)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

ggsave("./Figures/General/2016_log_DNARNA_withinPCoA_distance.png", dr.dist, width = 20, height = 12, unit = "cm")


#---------------------#
#------- Done! -------#
# Move to next script #
#---------------------#
sessionInfo()


#########################################################################
#- Extra code

# try NMDS
set.seed(3)
# Stress values equal to or below 0.1
# Stress values equal to or below 0.05 indicate good fit
nmds <- metaMDS(pb.bray, k = 6, trymax = 1000) # k = 6, minimum dim fulfilling stress
# but 1000 iterations did not reach convergence...

# try PCA
#pb.mat <-decostand(pb.mat, method = "norm") # chord transformation
#pb.mat <-decostand(pb.mat, method = "hellinger") # hellinger transformation
# Transformations were explored, no major difference, go with simplest solution

pca<-prcomp(pb.mat,retx=T,center=T,scale.=F)
scores<-pca$x
loadings<-pca$rotation

screeplot(pca)
pca.pct<-100*round(summary(pca)$importance[2,],3)
barplot(pca.pct)

# extract % of variance explained by each PC
pca.res <- summary(pca)
# extract the proprtion of variance explained by each PCA axis
exp<-data.frame(pca.res$importance)
exp <- exp[,1:2]

# site coordinates
pca.sites <- data.frame(PC1 = scores[,1]/pca$sdev[1], PC2 = scores[,2]/pca$sdev[2]) 

# species coordinates
pca.species <- loadings*matrix(pca$sdev,nrow=nrow(loadings),ncol=ncol(loadings),byrow=TRUE)
#pca.species <- pca.species *2 # choose extension factor
pca.species <- data.frame(PC1 = pca.species[,1], PC2 = pca.species[,2], Species = row.names(pca.species),
                          stringsAsFactors = F)
# loadings are weighted by sqrt(eigenvalues) (multiplied by sqrt(eigenvalues))

# export PCA scores and merge with meta data
ord.df<-data.frame(PC1=pca.sites$PC1,PC2=pca.sites$PC2, Samples = row.names(pca.sites),
                   stringsAsFactors = F)
ord.df <- merge(ord.df, sample_df(dna) %>% 
                  mutate(Samples = row.names(sample_df(dna))) %>%
                  dplyr::select(Samples, Year, Season, sample.type.year, DnaType), by = "Samples")

# set factors for plotting
ord.df$sample.type.year <- factor(ord.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheicwater", 
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver","RO3", "RO2", "RO1","Deep",
                                                                      "Downriver",
                                                                      "Marine"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Hyporheicwater", 
                                             "Groundwater","Stream", "Tributary",
                                             "Headwater \nLakes", "Upstream \nPonds", "Lake", "Lake",
                                             "Upriver","RO3","RO2", "RO1","Hypolimnion","Downriver",
                                             "Estuary"))
ord.df$Season <- factor(ord.df$Season, levels = c("spring", "summer", "autumn"), 
                        labels = c("Spring", "Summer","Autumn"))
ord.df$DnaType <- factor(ord.df$DnaType, levels = c("DNA", "cDNA"), labels = c("DNA", "RNA"))

# extract sample types represented in subset
colvec <- colvec[names(colvec) %in% as.character(levels(ord.df$sample.type.year))]

# get legend of plot separately
(
  bot.leg <-
    get_legend(
      ggplot() +
        theme_bw() +
        geom_point(data = ord.df, aes(x = PC1, y = PC2,
                                      fill = sample.type.year, shape = Season),
                   size = 2.5) +
        scale_fill_manual(values = colvec, name = "Sample Type") +
        theme(legend.position = "bottom") +
        scale_shape_manual(values = c(21, 23, 25)) +
        scale_alpha_manual(values = c(1, 0.5), name = "Nucleic Acid Type") +
        guides(shape = guide_legend(order = 1),
               alpha = guide_legend(order = 2), fill = "none")
    )
)

(pca.bi <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(data = ord.df, aes(x = PC1, y = PC2,
                                  fill = sample.type.year, shape = Season),
               size = 2.5) +
    #geom_text(data = pca.species, 
    #          aes(x = PC1, y = PC2, label = Species),
    #          size = 4) +
    #coord_fixed(ratio = .75) +
    scale_fill_manual(values = colvec, name = "Ecosystem Type") +
    scale_shape_manual(values = c(21,23,25)) +
    labs(x = paste("PC1 [", round(exp$PC1[2] * 100,2),"%]"), 
         y = paste("PC2 [", round(exp$PC2[2] * 100,2),"%]")) +
    theme(legend.key.size = unit(1, "lines"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right", 
          legend.box = "vertical",
          legend.margin = margin()) +
    guides(fill = guide_legend(override.aes=list(shape=21, size = 1.7), order = 1),
           shape = guide_legend(order = 2, override.aes = list(size = 1.7))))

ggsave(paste0("./Figures/Final/PCA_DNA_SampleType.tiff"), pca.bi,
       width = 15, height = 13, unit = "cm")
ggsave(paste0("./Figures/Final/PCA_DNA_SampleType.png"), pca.bi,
       width = 15, height = 13, unit = "cm")

adonis(pb.mat ~ sample.type.year + Season, data = ord.df, method = "eu")

# Get species contribution to each axes
# 10 species contributing maximum to the positive and negative of each axes
# keep the ASVs with the highest contribution to PC1 and PC2
min.pc1 <- head(pca.species[order(pca.species$PC1),], 10)$Species
max.pc1 <- head(pca.species[order(pca.species$PC1, decreasing = T),], 10)$Species
min.pc2 <- head(pca.species[order(pca.species$PC2),], 10)$Species
max.pc2 <- head(pca.species[order(pca.species$PC2, decreasing = T),], 10)$Species

contrib.pc1 <- c(min.pc1, max.pc1)
contrib.pc2 <- c(min.pc2, max.pc2)

species.names <- as.data.frame(tax_mat(dna)[row.names(tax_mat(dna)) %in% c(contrib.pc1, contrib.pc2),])
species.names$ASV <- row.names(species.names)

pca.species <-
  rbind(pca.species[pca.species$Species %in% contrib.pc1,] %>% mutate(Contrib = "PC1"),
        pca.species[pca.species$Species %in% contrib.pc2,] %>% mutate(Contrib = "PC2"))

contrib.sp <- merge(species.names, pca.species, by.x = "ASV", by.y = "Species")
setDT(contrib.sp)

contrib.sp[!is.na(genus), labels := paste0(order,", ", family,"\n" ,genus)]
contrib.sp[is.na(genus), labels := paste0(order,"\n", family)]

contrib.sp[, c("sum.pc1", "sum.pc2") := list(sum(PC1), sum(PC2)), by = .(labels)]

(pc1.bar <- ggplot(contrib.sp[order(contrib.sp$PC1) & Contrib == "PC1",], 
                   aes(x = reorder(labels, sum.pc1), y = sum.pc1)) +
    coord_flip() +
    theme_pubr() +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(fill = "white", colour = "black") +
    labs(x = "", y = "Contribution to PC1") +
    theme(axis.text.y = element_text(size = 10, face = "italic")))

(pc2.bar <- ggplot(contrib.sp[order(contrib.sp$PC2) & Contrib == "PC2",], 
                   aes(x = reorder(labels, sum.pc2), y = sum.pc2)) +
    coord_flip() +
    theme_pubr() +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(fill = "white", colour = "black") +
    labs(x = "", y = "Contribution to PC2") +
    theme(axis.text.y = element_text(size = 10, face = "italic")))

ggsave("./Figures/Final/Contribution_species_PCA_DNA.png", p, width = 10, height = 6)

