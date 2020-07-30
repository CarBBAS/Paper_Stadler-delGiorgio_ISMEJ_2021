###---------------------------------------------###
#-   Script for all plots in the publication:   - #
#-   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   - #
###---------------------------------------------###

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
# correct one miscategorisation
met.df[met.df$DadaNames == "RO2111.60mD",]$sample.type.year <- "Deep"

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
# Therefore, we need to select an asymmetrical similarity distance coefficient
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

# PCoA was chosen for all ordinations for consistency

############
#---------------------------------------------------------------------------------------------------#
############
# Only DNA #
############
# subset only DNA samples
dna <- subset_samples(pb, DnaType == "DNA")

# extract ASV matrix
pb.mat <- t(otu_mat(dna))

# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)
# plot with custom function
dna.pcoa <- plot_bray(pb.bray.pcoa, .id = "DNA", colours = colvec, output = T)

# Test for significant difference between factors
ord.df <- dna.pcoa[["df"]]
adonis(pb.mat ~ sample.type.year + Season, data = ord.df, 
       sqrt.dist = T, method = "bray")
#Permutation: free
#Number of permutations: 999
#Terms added sequentially (first to last)
#                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  sample.type.year  16    50.986  3.1866  14.252 0.36475  0.001 ***
#  Season             2     4.726  2.3629  10.567 0.03381  0.001 ***
#  Residuals        376    84.073  0.2236         0.60145           
#Total            394   139.784                 1.00000

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
# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)

# plot with custom function
all.pcoa <- plot_bray(pb.bray.pcoa, .id = "All", colours = colvec, output = T)

# statistically test if factors are different
# extract ordination coordinates and factors
ord.df <- all.pcoa[["df"]]
# run PERMANOVA
adonis(pb.mat ~ sample.type.year + Season + DnaType, data = ord.df, 
                      sqrt.dist = T, method = "bray", parallel = cl)
#Permutation: free
#Number of permutations: 999
#Terms added sequentially (first to last)
#                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  sample.type.year  16    64.493  4.0308  17.201 0.30217  0.001 ***
#  Season             2     7.440  3.7199  15.874 0.03486  0.001 ***
#  DnaType            1     6.288  6.2883  26.834 0.02946  0.001 ***
#  Residuals        577   135.213  0.2343         0.63351           
#Total            596   213.434                 1.00000

# plot, arrange legends separately to allow different aesthetics in different columns
(p <- ggarrange(ggarrange(dna.pcoa$plot, 
          all.pcoa$plot,
          widths = c(0.5, 0.5),
          ncol = 2, 
          common.legend = T,
          align = "h",
          labels = "auto",
          font.label = list(size = 20)),
  ggarrange(all.pcoa$legends[[1]],
                        all.pcoa$legends[[2]], ncol = 2),
          widths = c(0.7,0.18))
)

# save
ggsave(paste0("./Figures/Final/PCoAs_SampleType.tiff"), p,
       width = 30, height = 12, unit = "cm")
ggsave(paste0("./Figures/Final/PCoAs_SampleType.png"), p,
       width = 30, height = 12, unit = "cm")


# Save 3D plot for all PCoA, third axis separates seasons
plot.df <- all.pcoa[["df"]]
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