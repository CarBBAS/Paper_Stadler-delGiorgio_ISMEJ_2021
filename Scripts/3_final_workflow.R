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
library(mvabund) # multivariate glm

#-----------#
# FUNCTIONS #
#-----------#
source("./Functions/custom_fun.R")
source("./Functions/SET_framework.R")

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
asv.tab <- select_newest("./Output", "201520162017_fin_css_otu99_table_")
asv.tab <- as.matrix(read.csv(
  paste0("./Output/", asv.tab),
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))

# row orders need to match between tax.tab and asv.tab
asv.tab <- asv.tab[order(row.names(asv.tab)),]

# Taxonomy table
tax.tab <- select_newest("./Output", "201520162017_tax_otu99_table_")
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
  select_newest(path = "./Output", file.pattern = "201520162017_meta_otu99_data_")
met.df <-
  read.csv(
    paste0("./Output/", met.df),
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  )


# merge some sample types
met.df$sample.type.year <- factor(met.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheicwater", 
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver", "Downriver","RO3", "RO2", "RO1","Deep",
                                                                      "Marine"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Soilwater", 
                                             "Groundwater","Stream", "Tributary",
                                             "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
                                             "Upriver","Downriver",
                                             "Reservoirs","Reservoirs", "Reservoirs","Reservoirs",
                                             "Estuary"))

met.df$Season <- factor(met.df$Season, levels = c("spring","summer","autumn"),
                        labels = c("Spring","Summer","Autumn"))

met.df$DnaType <- factor(met.df$DnaType, levels = c("DNA","cDNA"),
                        labels = c("DNA","RNA"))

# Construct phyloseq object
pb <- phyloseq(otu_table(asv.tab, taxa_are_rows = F),
               sample_data(met.df),
               tax_table(tax.tab))

#####################################################################################################
# For rarefied table

# full tidy data set
#rel.df <- select_newest("./Objects", "201520162017_css")
#rel.df <- readRDS(
#  paste0("./Objects/", rel.df))

## RAREFIED DATASET ##
# read in rarefied datasets
min_lib <- c(15000, 25000, 50000)
perm.rar <- select_newest("./Output",
                          "perm.rar_lib", by = min_lib)

perm.rar <- c(perm.rar, select_newest("./Output", "201520162017_CSS_asvtab"))

rar <- read.csv(paste0("./Output/", perm.rar[2]), sep = ";", stringsAsFactors = F)

data.table::setDT(rar)
# make mean column to count
rar <- setDF(dcast(rar, Sample ~ ASV, value.var = "iter.mean"))
rownames(rar) <- rar$Sample
rar$Sample <- NULL
rar < t(as.matrix(rar))
rar <- rar[order(row.names(rar)),]

tax.tab <- tax.tab[row.names(tax.tab) %in% colnames(rar),]
tax.tab <- tax.tab[order(match(row.names(tax.tab), colnames(rar))),]
met.df <- met.df[row.names(met.df) %in% row.names(rar),]
met.df <- met.df[order(match(row.names(met.df), row.names(rar))),]

pb <- phyloseq(otu_table(t(rar), taxa_are_rows = T),
               sample_data(met.df),
               tax_table(tax.tab))
# remove taxa without observation
pb <- prune_taxa(!taxa_sums(pb) == 0, pb)

#####################################################################################################


# create colour vector for later plotting
# ensure consistent colours for all sample types

# export factors for colouring
sample.factors <- levels(met.df$sample.type.year)

# more colour blind friendly
colvec <- c("#FCFDBFFF", #"#FEC589FF", #Soil
            "#FEC589FF", #"#FDEBACFF", #Sediment
            "#F9795DFF", #"#F9795DFF", #Soilwater
            "#DE4968FF", #"#DE4968FF", #Groundwater,
            "skyblue", #Stream
  "#AD347CFF",# Tributary, 
  "palegreen", #Riverine Lakes, 
  "#7AD151FF", #Headwater Ponds,
  "#FDE725FF",# Lake, 
  "#1F9F88FF", # Upriver, 
  "#375A8CFF", #Downriver,
  "orchid", #"#471063FF", #Reservoir, 
  "#050416FF") #Estuary)

names(colvec) <- as.character(sample.factors)

# create colour vector for plotting
#colvec <- c("red4","chocolate3","orangered2","orange3",
#            "khaki","cadetblue","darksalmon",
#            "darkolivegreen","darkolivegreen3","gold",
#            "royalblue","mediumorchid4", "violet","palevioletred2","navy","skyblue",
#            "seagreen3")

# old palette
#colvec <- c("red4","chocolate3","tomato4",
#            "khaki","cadetblue","darksalmon",
#            "darkolivegreen","darkolivegreen3","gold",
#            "royalblue", "orchid","skyblue",
#           "navy")


############
# Analysis #
############
####################################################################################

#- Figure 1 -#

# Q: Is the microbial assemblage different between habitat types and seasons?

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

#---------------------------------------------------------------------------------------------------#

############
# Only DNA #
############
# subset only DNA samples
dna <- subset_samples(pb, DnaType == "DNA")

# extract ASV matrix
pb.mat <- otu_mat(dna)
#pb.mat <- log2(pb.mat + 1)
# PCoA with Bray-Curtis

meta <- data.frame(Sample = as.character(row.names(sample_df(dna))),
                    sample_df(dna) %>% dplyr::select(sample.type.year, Season, Year, DnaType), 
                    stringsAsFactors = F)

# melt to calculate mean variance relationship
melt.mat <- melt.data.table(
  setDT(as.data.frame(pb.mat), keep.rownames = "Sample"),
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

plot.df <- melt.mat[, .(mean = mean(reads, na.rm = T),
             variance = var(reads, na.rm = T)), by = .(OTU)]

ggplot(plot.df, aes(x = log1p(mean), y = log(variance))) +
  geom_point()

# ASVs with high means also have high variances

#ord.asv <- plot.df$ASV[order(plot.df$mean, decreasing = T)]
#melt.mat$ASV <- factor(melt.mat$ASV, levels = ord.asv)
#melt.mat <- melt.mat[meta, c("sample.type.year", "Season") := list(i.sample.type.year,
#                                                               i.Season), on = .(Sample)]
#ggplot(melt.mat, aes(x = ASV, y = log2(reads + 1), colour = sample.type.year)) +
#  geom_point()


# make mvabund object of community matrix
dna.sp <- mvabund(pb.mat)

mod <- manyglm(dna.sp ~ meta$sample.type.year * meta$Season, family = "negative.binomial")
# warning but is integer
saveRDS(mod, "./Objects/manyglm.dna.negbinom.rds")

# check residuals, it's not optimal, but compared to other families, there is less of a pattern
png(filename="./Figures/General/manyglm_dna_residuals_binom.png")
plot(mod)
dev.off()

# test for habitat type and season effect
anova <- anova(mod)

pb.mat <- decostand(pb.mat, "hellinger")
pb.mori <- vegdist(pb.mat, method = "horn")
is.euclid(pb.mori) # FALSE
pb.mori <- sqrt(pb.mori) # make Euclidean
is.euclid(pb.mori) # TRUE

pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.mori)
# plot with custom function
dna.pcoa <- plot_bray(pb.bray.pcoa, .id = "DNA", colours = colvec, output = T)

p <- dna.pcoa$plot + guides(alpha = "none")

# save
ggsave(paste0("./Figures/Final/PCoA_log_DNA_SampleType.tiff"), p,
       width = 12, height = 10, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_log_DNA_SampleType.png"),  p,
       width = 12, height = 10, unit = "cm")


# PERMANOVA is sensitive towards unbalanced sampling designs
# First option: Remove samples from groups that have many more samples
ord.df <- dna.pcoa[["df"]]
ord.df$groups <- paste(ord.df$sample.type.year, ord.df$Season, sep = "_")

setDT(ord.df)

numbers <- ord.df[, .(n = .N), by = .(groups)]

too.many <- numbers$groups[numbers$n >= 3]
# we will drop all habitat type ~ season combination that have less than three samples
set.seed(3)
random <- ord.df[groups %in% too.many, .(Sample = sample(Sample, size = 3, replace = F)), by = .(groups)]$Sample

# redo PCoA
pb.bray <- vegdist(pb.mat[rownames(pb.mat) %in% random,], method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean

# Test for significant difference between factors
adonis(pb.mat[rownames(pb.mat) %in% random,] ~ sample.type.year * Season, data = ord.df[Sample %in% random,], 
       sqrt.dist = T, method = "bray")

adonis(pb.mat ~ sample.type.year * Season, data = ord.df, 
       sqrt.dist = T, method = "horn")

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#sample.type.year        12   10.8058 0.90048  4.1826 0.39986  0.001 ***
#  Season                   2    1.3304 0.66518  3.0896 0.04923  0.001 ***
#  sample.type.year:Season 12    3.2618 0.27182  1.2625 0.12070  0.008 ** 
#  Residuals               54   11.6259 0.21529         0.43021           
#Total                   80   27.0238                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# calculate multivaraite dispersions
mod <- betadisper(pb.bray, group = ord.df[Sample %in% random,]$groups, bias.adjust = T)
mod <- betadisper(pb.mori, group = ord.df$groups, bias.adjust = T)

## Perform test
anova(mod) # significant..., homogeneity of variance not fullfilled

## Permutation test for F
pmod <- permutest(mod, permutations = 999, pairwise = T)

TukeyHSD(mod)
plot(mod, label.cex = 0.5)
boxplot(mod)

# Plot for supplements
# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)
dna.pcoa <- plot_bray(pb.bray.pcoa, .id = "DNA", colours = colvec, output = T)

p <- dna.pcoa$plot + guides(alpha = "none")
p


############
# Only RNA #
############
# subset only DNA samples
rna <- subset_samples(pb, DnaType == "cDNA")

# extract ASV matrix
pb.mat <- t(otu_mat(rna))
pb.mat <- log2(pb.mat + 1)
# PCoA with Bray-Curtis
#pb.mat <- decostand(pb.mat, "hellinger")

pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)
# plot with custom function
rna.pcoa <- plot_bray(pb.bray.pcoa, .id = "DNA", colours = colvec, output = T)

p <- rna.pcoa$plot + guides(alpha = "none")

# save
ggsave(paste0("./Figures/Final/PCoA_log_RNA_SampleType.tiff"), p,
       width = 12, height = 10, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_log_RNA_SampleType.png"),  p,
       width = 12, height = 10, unit = "cm")

# Test for significant difference between factors
ord.df <- rna.pcoa[["df"]]
ord.df$groups <- paste(ord.df$sample.type.year, ord.df$Season, sep = "_")
adonis(pb.mat ~ sample.type.year * Season, data = ord.df, 
       sqrt.dist = T, method = "bray")

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#sample.type.year  10    17.787 1.77874  8.2872 0.28759  0.001 ***
#  Season             2     3.496 1.74812  8.1446 0.05653  0.001 ***
# sample.type.year:Season  13     4.919 0.37841  1.8683 0.07954  0.001 ***
#  Residuals               176    35.647 0.20254         0.57635           
#Total                   201    61.850                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# calculate multivaraite dispersions
mod <- betadisper(pb.bray, group = ord.df$groups)

## Perfmom test
anova(mod) #highly significant

## Permutation test for F
pmod <- permutest(mod, permutations = 999, pairwise = T)

TukeyHSD(mod)
plot(mod, label.cex = 0.5)
boxplot(mod)

####################
#---------------------------------------------------------------------------------------------------#
####################
# Both DNA and RNA #
####################
# extract species table with species in columns
pb.mat <- t(otu_mat(pb))

# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE

pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray)
# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)

# plot with custom function
all.pcoa <- plot_bray(pb.bray.pcoa, .id = "All", colours = colvec, output = T)
pcoa.23 <-  plot_bray(pb.bray.pcoa, .id = "All", colours = colvec, output = T, axes = "2+3")

pcoa.plot <- all.pcoa$plot + theme(legend.position = "left")

p <- ggarrange(all.pcoa$plot, pcoa.23$plot, ncol = 2, common.legend = T, legend = "right",
               align = "hv", labels = "auto")

# save
ggsave(paste0("./Figures/Final/PCoA_all_SampleType.tiff"), p,
       width = 20, height = 11, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_all_SampleType.png"),  p,
       width = 20, height = 11, unit = "cm")

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

#---------------------------------------------------------------------------------------------------#

#- Figure 2 -#

# Q2: How different are the the DNA and RNA assemblages of the same sample?

#-------------------------------------------------------#
# Extract pair-wise dissimilarity among DNA-RNA samples #
#-------------------------------------------------------#
# UPDATE: Instead of using distance between DNA and RNA wihin PCoA space, we use dissimilarity
# This approach is favoured as PCoA reduces sample complexity into two dimensions and misses variation
# Pair-wise dissimilarity is a more direct approach to grasp how different DNA and RNA of the same sample are

dissim.dr <- dissim.dnarna(pb, save.name = "All", output = T)

#-------------------------------------------------------------#
# Calculate distance between DNA and RNA points in PCoA space #
#-------------------------------------------------------------#
# use custom function to correct a few wrong sample names and match DNA-RNA counterpart samples
# calculating distance between points in two-dimensional space for both Bray-Curtis and Jaccard
# Use three distances, as variances explained by second and third axes are almost identical
dist.dr <- dist.dnarna(all.pcoa[["df"]], save.name = "3D", dimensions = 3, output = T)
#dist.dr <- dist.dnarna(dnarna.bray[["df"]], save.name = "All_2D", dimensions = 2)

# dissimilarity and distance show different patterns....
# what is behind this difference?

(p <-ggarrange(ggarrange(dist.dr$plot.main  + theme(axis.text.x = element_blank()), 
            dist.dr$plot.side  + theme(axis.text.x = element_blank()),
            widths = c(3,1),
            ncol = 2, nrow = 1, 
            common.legend = T,
            legend = "none",
            align = "h",
            font.label = list(size = 10)),
  ggarrange(dissim.dr$plot.main, 
            dissim.dr$plot.side,
            widths = c(3,1),
            ncol = 2, nrow = 1, 
            common.legend = T,
            legend = "none",
            align = "h",
            font.label = list(size = 10)),
  nrow = 2, heights = c(0.4, 0.5), labels = c("auto"), common.legend = T, align = "hv",
  legend.grob = get_legend(dist.dr$plot.main), legend = "right"))
# labels = c("b","c")
# add x axis title to be in the middle of two panels
(p <- annotate_figure(p, bottom = text_grob("Habitat Type")))

# save
ggsave(paste0("./Figures/Final/DNARNA_distdissim.tiff"), p,
       width = 17, height = 10, unit = "cm")
ggsave(paste0("./Figures/Final/DNARNA_distdissim.png"),  p,
       width = 17, height = 10, unit = "cm")

#(fin <- ggarrange(pcoa.plot, p,
#          ncol = 2, labels = c("a","")))

# save
ggsave(paste0("./Figures/Final/DNARNA_PCoA_distdissim.tiff"), fin,
       width = 32, height = 12, unit = "cm")
ggsave(paste0("./Figures/Final/DNARNA_PCoA_distdissim.png"),  fin,
       width = 32, height = 12, unit = "cm")

#---------------------------------------------------------------------------------------------------#

#- Figure 3 -#

# Q2: What is causing the dissimilarity patterns between DNA and RNA?

# Create abundance heatmap by abundance group
# 1. Define abundance groups

# Conventional grouping does not work with this data set
# Otherwise, there are no abundant ASVs as we cover too different ecosystem types

# Try to come up with a new classification
# 1. Calclate the mean abundance of each ASV for a sample type (e.g. reservoir, lake, stream, soil etc)
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

# After consulting the means and deviations of the maximum and minimum thresholds across samples
# (see previous script) we settle with:
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

# create new ID column with sample type and DNA Type
sum[, ID := paste(sample.type.year, DnaType, sep = "_")]

tax.order <- llply(ls.asvs, function(x){
  # extract abundance group specific ASVs
  df <- sum[ASV %in% x,]
  com.mat <- dcast(df, ID ~ ASV, value.var = "mean.css")
  setDF(com.mat)
  rownames(com.mat) <- com.mat$ID; com.mat$ID <- NULL # remove ID col
  # make PCoA to get ASV order in heatmap
  pb.bray <- vegdist(com.mat, method = "bray")
  pb.bray <- sqrt(pb.bray) # make Euclidean
  bray.pcoa <- ape::pcoa(pb.bray)
  # extract scores of first three axes
  sitesDF <- data.frame(bray.pcoa$vectors[,1:3])
  speciesDF <- data.frame(wascores(sitesDF, w = com.mat), stringsAsFactors = F)
  speciesDF <- na.omit(speciesDF)
  taxa.order <- row.names(speciesDF)[order(RadialTheta(speciesDF))]
  return(taxa.order)
}, .parallel = T)

sum$ID <- factor(sum$ID, levels = c("Soil_DNA", "Soil_cDNA",
                                    "Sediment_DNA", "Sediment_cDNA",
                                    "Soilwater_DNA","Soilwater_cDNA",
                                    "Hyporheicwater_DNA", "Hyporheicwater_cDNA",
                                    "Stream_DNA", "Stream_cDNA",
                                    "Tributary_DNA", "Tributary_cDNA",
                                    "HeadwaterLakes_DNA", "HeadwaterLakes_cDNA",
                                    "Lake_DNA", "Lake_cDNA",
                                    "Upriver_DNA", "Upriver_cDNA",
                                    "RO3_DNA", "RO3_cDNA",
                                    "RO2_DNA", "RO2_cDNA",
                                    "RO1_DNA", "RO1_cDNA",
                                    "Deep_DNA", "Deep_cDNA",
                                    "Downriver_DNA", "Downriver_cDNA",
                                    "Marine_DNA", "Marine_cDNA"),
                 labels = c("Soil DNA", "Soil RNA",
                            "Sediment DNA", "Sediment RNA",
                            "Soilwater DNA","Soilwater RNA",
                            "Hyporheicwater DNA", "Hyporheicwater RNA",
                            "Stream DNA", "Stream RNA",
                            "Tributary DNA", "Tributary RNA",
                            "Headwater Lakes DNA", "Headwater Lakes RNA",
                            "Lake DNA", "Lake RNA",
                            "Upriver DNA", "Upriver RNA",
                            "RO3 DNA", "RO3 RNA",
                            "RO2 DNA", "RO2 RNA",
                            "RO1 DNA", "RO1 RNA",
                            "Hypolimnion DNA", "Hypolimnion RNA",
                            "Downriver DNA", "Downriver RNA",
                            "Estuary DNA", "Estuary RNA"))

x.labs<-c("D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R")

heat.plots <- llply(tax.order, function(x){
  #if(length(x) > 5000){
  #  sub.x <- sample(x, 5000)
  #  df <- sum[sum$ASV %in% sub.x,]
  #} else {
  #  
  #}
  df <- sum[sum$ASV %in% x,]
  df$ASV <- factor(df$ASV, levels = x)
  
  # get max value of all
  max.lim <- round(max(log2(sum$mean.css+1)),0)
  
  #low <- "#000033"; high <- "#66CCFF"
  p <- ggplot(df, aes(x = ID, y = ASV, fill = log2(mean.css + 1))) +
    theme_bw() +
    geom_raster() +
    labs(x = "Habitat Type", y = "ASVs") +
    scale_fill_viridis_c(option = "magma", name = "Mean\nReads",
                         limits = c(min(log2(df$mean.css + 1)),
                                    max.lim),
                         breaks = c(min(log2(df$mean.css + 1)),
                                    max.lim/ 2,
                                    max.lim),
                         labels = c(as.character(min(log2(df$mean.css + 1))),
                                    as.character(max.lim / 2),
                                    as.character(max.lim))) +
    scale_x_discrete(labels = x.labs) +
    coord_cartesian(ylim = c(1,
                             length(levels(df$ASV))), clip = "off") +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.title.x = element_blank(),
          legend.position = "right", plot.title = element_text(size = 12)) +
    guides(fill = guide_colourbar(barwidth = unit(10, "pt"),
                                  barheight = unit(50, "pt")))
    #scale_fill_gradient(low = low, high = high)
  return(p)
}, .parallel = T)

(univ.rare <- heat.plots$universal.rare + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Universally rare"))

(spec <- heat.plots$specialist + 
  theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank()
) + labs(title = "Specialist"))

(rare.shif <- heat.plots$rare.shifter +
    theme(plot.margin=unit(c(2,2,18,2), "mm")) +
  labs(title = "Rare shifter") +
    annotate("segment", x = 0.75, xend = 2.25, y = -30, yen = -30) +
    annotate("text", x = 1.6, y = - 40, label = "Soil", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 2.75, xend = 4.25, y = -30, yen = -30) +
    annotate("text", x = 3.6, y = - 40, label = "Sediment", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 4.75, xend = 6.25, y = -30, yen = -30) +
    annotate("text", x = 5.6, y = - 40, label = "Soilwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 6.75, xend = 8.25, y = -30, yen = -30) +
    annotate("text", x = 7.6, y = - 40, label = "Hyporheicwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 8.75, xend = 10.25, y = -30, yen = -30) +
    annotate("text", x = 9.6, y = - 40, label = "Stream", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 10.75, xend = 12.25, y = -30, yen = -30) +
    annotate("text", x = 11.6, y = - 40, label = "Tributary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 12.75, xend = 14.25, y = -30, yen = -30) +
    annotate("text", x = 13.6, y = - 40, label = "Headwater Lakes", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 14.75, xend = 16.25, y = -30, yen = -30) +
    annotate("text", x = 15.6, y = -40, label = "Lake", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 16.75, xend = 18.25, y = -30, yen = -30) +
    annotate("text", x = 17.6, y = -40, label = "Upriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 18.75, xend = 20.25, y = -30, yen = -30) +
    annotate("text", x = 19.6, y = -40, label = "RO3", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 20.75, xend = 22.25, y = -30, yen = -30) +
    annotate("text", x = 21.6, y = -40, label = "RO2", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 22.75, xend = 24.25, y = -30, yen = -30) +
    annotate("text", x = 23.6, y = -40, label = "RO1", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 24.75, xend = 26.25, y = -30, yen = -30) +
    annotate("text", x = 25.6, y = -40, label = "Hypolimnion", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 26.75, xend = 28.25, y = -30, yen = -30) +
    annotate("text", x = 27.6, y = -40, label = "Downriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 28.75, xend = 30.25, y = -30, yen = -30) +
    annotate("text", x = 29.6, y = -40, label = "Estuary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5))


(ab.shif <- heat.plots$abundant.shifter + 
  theme(axis.title.y = element_blank(),
        plot.margin=unit(c(2,2,18,2), "mm")) + 
  labs(title = "Abundant shifter") +
    annotate("segment", x = 0.75, xend = 2.25, y = -70, yen = -70) +
    annotate("text", x = 1.6, y = - 90, label = "Soil", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 2.75, xend = 4.25, y = -70, yen = -70) +
    annotate("text", x = 3.6, y = - 90, label = "Sediment", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 4.75, xend = 6.25, y = -70, yen = -70) +
    annotate("text", x = 5.6, y = - 90, label = "Soilwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 6.75, xend = 8.25, y = -70, yen = -70) +
    annotate("text", x = 7.6, y = - 90, label = "Hyporheicwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 8.75, xend = 10.25, y = -70, yen = -70) +
    annotate("text", x = 9.6, y = - 90, label = "Stream", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 10.75, xend = 12.25, y = -70, yen = -70) +
    annotate("text", x = 11.6, y = - 90, label = "Tributary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 12.75, xend = 14.25, y = -70, yen = -70) +
    annotate("text", x = 13.6, y = - 90, label = "Headwater Lakes", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 14.75, xend = 16.25, y = -70, yen = -70) +
    annotate("text", x = 15.6, y = - 90, label = "Lake", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 16.75, xend = 18.25, y = -70, yen = -70) +
    annotate("text", x = 17.6, y = - 90, label = "Upriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 18.75, xend = 20.25, y = -70, yen = -70) +
    annotate("text", x = 19.6, y = - 90, label = "RO3", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 20.75, xend = 22.25, y = -70, yen = -70) +
    annotate("text", x = 21.6, y = - 90, label = "RO2", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 22.75, xend = 24.25, y = -70, yen = -70) +
    annotate("text", x = 23.6, y = - 90, label = "RO1", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 24.75, xend = 26.25, y = -70, yen = -70) +
    annotate("text", x = 25.6, y = - 90, label = "Hypolimnion", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 26.75, xend = 28.25, y = -70, yen = -70) +
    annotate("text", x = 27.6, y = - 90, label = "Downriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 28.75, xend = 30.25, y = -70, yen = -70) +
    annotate("text", x = 29.6, y = - 90, label = "Estuary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5))

p <- ggarrange(univ.rare, spec, rare.shif, ab.shif,
          ncol = 2, nrow = 2, align = "v", common.legend = T, heights = c(0.45, 0.55),
          labels = "auto", legend = "right")
p

# what's the N ?
lapply(tax.order, length)
# universal.rare 22503
# specialist 91
# abundant.shifter 525
# rare.shifter 230

# save
ggsave(paste0("./Figures/Final/Abundance_profiles_heat.tiff"), p,
       width = 25, height = 15, unit = "cm")
ggsave(paste0("./Figures/Final/Abundance_profiles_heat.png"),  p,
       width = 25, height = 15, unit = "cm")


#---------------------------------------------------------------------------------------------------#

#- Figure 4 -#

# Q2: What is causing the dissimilarity patterns between DNA and RNA?

# fit species onto PCoA ordination
pb.species <-
  data.frame(
    ASV = as.character(colnames(pb.mat)),
    wascores(pb.bray.pcoa$vectors[, 1:3],
             w = pb.mat),
    stringsAsFactors = F
  )

pb.species <- na.omit(pb.species)
setDT(pb.species)
# add abundance group names into species dataframe
for(i in 1:length(grp.names)){
  pb.species[ASV %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
  pb.species[ASV %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
}


# main plot
(p <- ggplot() +
    theme_cust() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(data = all.pcoa$df,
               aes(x = Axis.1, y = Axis.2, 
                   shape = Season, alpha = DnaType),
               size = 2.5, fill = "grey40") +
    scale_shape_manual(values = c(21,23,25)) +
    scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
    geom_point(data = pb.species[ab.names == "universal.rare",], 
               aes(x = Axis.1, y = Axis.2, fill = ab.names), shape = 21) +
    geom_point(data = pb.species[ab.names != "universal.rare",], 
               aes(x = Axis.1, y = Axis.2, fill = ab.names), shape = 21) +
    #scale_fill_manual(values = colvec, name = "Habitat Type") +
    coord_fixed(1) + # ensure aspect ratio
    labs(x = paste("PC1 [", unique(all.pcoa$df$x),"%]"), 
         y = paste("PC2 [", unique(all.pcoa$df$y),"%]")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
           alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
           fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
)

ggsave("./Figures/General/PCoA_all_ab.groups_withunivrare_pc12.png", p,
       height = 10, width = 15, units = "cm")

# Correlate species to axes?

cor.mat <- cor(all.pcoa$df[,2:4], pb.mat, method = "pearson")
cor.mat <- setDT(as.data.frame(cor.mat), keep.rownames = "Axes")
cor.df <- melt.data.table(
  cor.mat,
  id.vars = "Axes",
  measure.vars = patterns("^ASV_"),
  variable.name = "ASV",
  value.name = "pear.cor"
)


cor.heat <- llply(tax.order, function(x){
  df <- cor.df[cor.df$ASV %in% x,]
  df$ASV <- factor(df$ASV, levels = x)
  
  # get max value of all
  max.lim <- round(max(cor.df$pear.cor),0)
  
  p <- ggplot(df, aes(x = Axes, y = ASV, fill = pear.cor)) +
    theme_bw() +
    geom_raster() +
    labs(x = "PCoA Axes", y = "ASVs") +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "red",
                         name = "Pearson's r") +
    scale_x_discrete(labels = c("PC1", "PC2", "PC3")) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.title.x = element_blank(),
          legend.position = "right", plot.title = element_text(size = 12)) +
    guides(fill = guide_colourbar(barwidth = unit(10, "pt"),
                                  barheight = unit(50, "pt")))
  
  return(p)
}, .parallel = T)

p <- ggarrange(cor.heat$universal.rare +labs(title="Universally rare"), 
               cor.heat$specialist +labs(title="Specialist"),
               cor.heat$rare.shifter +labs(title="Rare shifter"), 
               cor.heat$abundant.shifter+labs(title="Abundant shifter"),
               ncol = 2, nrow = 2, align = "v", common.legend = T,
               labels = "auto", legend = "right")
p

ggsave("./Figures/General/heatplot_ab.groups_PC_pearcor.png", p,
       height = 12, width = 12)
# PC1: negative = soil, positive = water
# PC2: negative = RNA, positive = DNA
# PC3: negative = summer/autumn, positive = spring

###########################################################################


# does not work on private computer
sim <- with(sample_df(pb), vegan::simper(t(otu_mat(pb)), DnaType, permutations = 100, parallel = cl))

#### try simper #######

cor.df <- dcast(cor.df, ASV ~ Axes, value.var = "pear.cor")
colnames(cor.df)[2:4] <- c("cor.pc1","cor.pc2","cor.pc3")

cor.df <- cor.df[pb.species, c("cordi.pc1","cordi.pc2","cordi.pc3",
                     "ab.groups") := list(i.Axis.1, i.Axis.2, i.Axis.3,
                                          i.ab.names), on = .(ASV)]


ggplot(cor.df, aes(x = cor.pc1, y = cor.pc2, fill = ab.groups)) +
  theme_bw() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
  geom_point(data = subset(cor.df, ab.groups == "universal.rare"), shape = 21, alpha = 0.5) +
  geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"), shape = 21) +
  geom_point(data = subset(cor.df, ab.groups == "rare.shifter"), shape = 21) +
  geom_point(data = subset(cor.df, ab.groups == "specialist"), shape = 21) +
  scale_fill_viridis_d() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Terrestrial",""
                   ,"","Water"),
  hjustvar = c(-0.1,0,1,1.2) ,
  vjustvar = c(-0.2,1,0,1.2))

#hjust: lower values right, larger values left
#vjust lower values up, higher values down

#Bottom Left (h0,v0)","Top Left (h0,v1)"
#,"Bottom Right h1,v0","Top Right h1,v1"

(p1 <- ggplot() +
  theme_bw() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
  geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
             aes(x = cordi.pc1, y = cor.pc1, colour = ab.groups), shape = 21, alpha = 0.5, colour = "grey80") +
  geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"),
             aes(x = cordi.pc1, y = cor.pc1, colour = ab.groups), shape = 21, alpha = 0.8) +
  geom_point(data = subset(cor.df, ab.groups == "rare.shifter"),
             aes(x = cordi.pc1, y = cor.pc1, colour = ab.groups), shape = 21, alpha = 0.8) +
  geom_point(data = subset(cor.df, ab.groups == "specialist"),
             aes(x = cordi.pc1, y = cor.pc1, colour = ab.groups), shape = 21, alpha = 0.8) +
    geom_smooth(data = cor.df, 
                aes(x = cordi.pc1, y = cor.pc1, group = ab.groups,
                     colour = ab.groups, fill = ab.groups), size = 1.5,
                method = 'gam', formula = y ~ s(x, bs = "cs")) +
  geom_text(data=annotations,aes(x=xpos,y=ypos,
                                 hjust=hjustvar,vjust=vjustvar,
                                 label=annotateText)) +
  scale_fill_viridis_d() +
  labs(x = "Species coordinates [PC1]", 
       y = expression(paste("Pearsons's ", italic("r")," of ASVs with PC1"))) +
  scale_colour_viridis_d() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
)


annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("RNA",""
                   ,"","DNA"),
  hjustvar = c(-0.1,0,1,1.2) ,
  vjustvar = c(-0.2,1,0,1.2))

(p2 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cordi.pc2, y = cor.pc2, colour = ab.groups), shape = 21, alpha = 0.5, colour = "grey80") +
    geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"),
               aes(x = cordi.pc2, y = cor.pc2, colour = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "rare.shifter"),
               aes(x = cordi.pc2, y = cor.pc2, colour = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "specialist"),
               aes(x = cordi.pc2, y = cor.pc2, colour = ab.groups), shape = 21, alpha = 0.8) +
    geom_smooth(data = cor.df, 
                aes(x = cordi.pc2, y = cor.pc2, group = ab.groups,
                    colour = ab.groups, fill = ab.groups), size = 1.5,
                method = 'gam', formula = y ~ s(x, bs = "cs")) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    scale_fill_viridis_d() +
  labs(x = "Species coordinates [PC2]", 
       y = expression(paste("Pearsons's ", italic("r")," of ASVs with PC2"))) +
  scale_colour_viridis_d() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  )


annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Summer-Autumn",""
                   ,"","Spring"),
  hjustvar = c(-0.05,0,1,1.2) ,
  vjustvar = c(-0.2,1,0,1.2))

(p3 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cordi.pc3, y = cor.pc3, colour = ab.groups), shape = 21, alpha = 0.5, colour = "grey80") +
    geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"),
               aes(x = cordi.pc3, y = cor.pc3, colour = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "rare.shifter"),
               aes(x = cordi.pc3, y = cor.pc3, colour = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "specialist"),
               aes(x = cordi.pc3, y = cor.pc3, colour = ab.groups), shape = 21, alpha = 0.8) +
    geom_smooth(data = cor.df, 
                aes(x = cordi.pc3, y = cor.pc3, group = ab.groups,
                    colour = ab.groups, fill = ab.groups), size = 1.5,
                method = 'gam', formula = y ~ s(x, bs = "cs")) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    scale_fill_viridis_d() +
    labs(x = "Species coordinates [PC3]", 
         y = expression(paste("Pearsons's ", italic("r")," of ASVs with PC3"))) +
    scale_colour_viridis_d() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)


annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Terrestrial - RNA","Terrestrial - DNA"
                   ,"Water - RNA","Water - DNA"),
  hjustvar = c(-0.05, -0.05, 
               1.05, 1.05) ,
  vjustvar = c(-0.2, 1.2, 
               -0.2, 1.2))

(pp1 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.2) +
    geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"), 
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "rare.shifter"),
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "specialist"), 
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    labs(x = expression(paste("Pearsons's ", italic("r")," of ASVs with PC1")), 
         y = expression(paste("Pearsons's ", italic("r")," of ASVs with PC2"))) +
    scale_fill_viridis_d() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Terrestrial - Summer","Terrestrial - Spring"
                   ,"Water - Summer","Water - Spring"),
  hjustvar = c(-0.05, -0.05, 
               1.05, 1.05) ,
  vjustvar = c(-0.2, 1.2, 
               -0.2, 1.2))

(pp2 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.2) +
    geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"), 
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "rare.shifter"),
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "specialist"), 
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.8) +
    labs(x = expression(paste("Pearsons's ", italic("r")," of ASVs with PC1")), 
         y = expression(paste("Pearsons's ", italic("r")," of ASVs with PC3"))) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    scale_fill_viridis_d() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

p.cordcor <- ggarrange(p1, p2, p3, ncol = 3, common.legend = T, align = "hv")

p.corcor <- ggarrange(pp1, pp2, ncol = 2, common.legend = T, align = "hv")

ggsave("./Figures/Final/abundance_grps_pear.cor_cord.png", p.cordcor,
       height = 12, width = 30, units = "cm")

ggsave("./Figures/Final/abundance_grps_pear.cors_pcs.png", p.corcor,
       height = 10, width = 20, units = "cm")


##########################################################################
#-----------------------#
# Correlate to richness #
#-----------------------#
# read in rarefied datasets

alpha.df <- read.csv("./Output/alpha_div_summary.csv", sep = ";", dec = ".", stringsAsFactors = F)

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

# get data from All PCoA
pdataframe <- dist.dr$raw.df

# get ID
setDT(alpha.df); setDT(pdataframe)
alpha.df[pdataframe, c("ID","DnaType") := 
           list(i.ID, i.DnaType), on = .(Sample)]

# split DNA and RNA
dna.alpha <- alpha.df[alpha.df$DnaType == "DNA",]
rna.alpha <- alpha.df[alpha.df$DnaType == "RNA",]

###
# DNA
###

# Execute regressions first on DNA diversity
df <- melt(dna.alpha, id.vars = c("ID","Data"),
                  measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
                  variable.name = "Index",
                  value.name = "Diversity")

# Only focus on CSS, non rarefied data and two diversity indices
df <- df[Data == "css" & (Index == "Shannon" | Index == "Pielou"),]

# calculate mean of duplicates
df[, .(Diversity = mean(Diversity, na.rm = T)), by = .(ID, Index)]

# merge with distance
reg.df <- df[dist.dr$indiv.df, c("distance", "sample.type.year") := 
                       list(i.distance, i.sample.type.year), on = .(ID)]

# remove any NAs, e.g. samples from 2015
reg.df <- reg.df[!is.na(distance),]
setorderv(reg.df, c("Index","Diversity")) # rearrange dataframe

# decide which model is best (e.g. linear or polynomial)
z <- reg.df[Index == "Shannon"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
      distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 0.67
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.73
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.70
anova(lm0,lm1) # preferred model is lm1, higher R2 and lowered RSS
# polynomials do not add much, and avoid overfitting and choose a parsimonious model
anova(lm1,lm2) 
anova(lm2,lm3) 
rm(z)


z <- reg.df[Index == "Pielou"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
      distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 0.70
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.85
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.81
anova(lm0,lm1) 
anova(lm1,lm2) # preferred model is lm2, higher R2 and lowered RSS, and significance < 0.01
anova(lm2,lm3) 
rm(z)

lin.ls <- dlply(reg.df, .(Index), function(z){
  setDT(z)
  means <- z[, .(Diversity = mean(Diversity, na.rm = T),
                      distance = mean(distance, na.rm = T)), by = .(sample.type.year)]
  if(unique(z$Index) == "Shannon"){
    lin <- lm(means$Diversity ~ means$distance)
  } else if(unique(z$Index) == "Pielou") {
    lin <-  lm(means$Diversity ~ poly(means$distance, 2))
  }
  
  # check linear assumptions
  #plot(lin) # normality not good
  # large sample sizes, normality does not affect results too much (central limit theorem)
  # homoscedasticity and independence important
  #summary(lin)
  #confint(lin, level = 0.95)
  
  # get data for plotting
  x <- data.frame(x = sort(means$distance))
  pred <- predict(lin, newdata = x, se = T)
  ci <- pred$se.fit[order(means$distance)] * qt(0.95 / 2 + 0.5, pred$df)
  y <- pred$fit[order(means$distance)]
  ymin <- y - ci
  ymax <- y + ci
  
  plot.df <- data.frame(x = sort(means$distance),
                        y = y,
                        ymin = ymin,
                        ymax = ymax,
                        se = pred$se.fit[order(means$distance)])
  
  # extract only colours that are in data frame
  colvec <- colvec[names(colvec) %in% as.character(levels(means$sample.type.year))]
  
  p <- ggplot() +
    theme_cust("pubr") +
    geom_point(data = z, aes(x = distance, y = Diversity), 
               colour = "gray40", alpha = 0.3, size = 2) +
    geom_line(data = plot.df, aes(x = x, y = y), inherit.aes = F, size = 1) +
    geom_line(data = plot.df, aes(x = x, y = ymax), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_line(data = plot.df, aes(x = x, y = ymin), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_point(data = means,
               aes(x = distance, y = Diversity, fill = sample.type.year), shape = 21, size = 3) +
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
       fitted = plot.df,
       plot = p)
  
})

(p <- ggarrange(lin.ls[[1]]$plot, lin.ls[[2]]$plot,
                ncol = 2, common.legend = T, legend = "right"))

(p <- annotate_figure(p, 
                      bottom = text_grob("Distance in ordination space",
                                         just = "centre")))
ggsave("./Figures/Final/Richness_distance_nonlin_reg.png", p,
       width = 22, height = 11, unit = "cm")


####
## RNA
####

# Execute regressions now with RNA diversity
df <- melt(rna.alpha, id.vars = c("ID","Data"),
           measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
           variable.name = "Index",
           value.name = "Diversity")

# Only focus on CSS, non rarefied data and two diversity indices
df <- df[Data == "css" & (Index == "Shannon" | Index == "Pielou"),]

# calculate mean of duplicates
df[, .(Diversity = mean(Diversity, na.rm = T)), by = .(ID, Index)]

# merge with distance
reg.df <- df[dist.dr$indiv.df, c("distance", "sample.type.year") := 
               list(i.distance, i.sample.type.year), on = .(ID)]

# remove any NAs, e.g. samples from 2015
reg.df <- reg.df[!is.na(distance),]
setorderv(reg.df, c("Index","Diversity")) # rearrange dataframe

# decide which model is best (e.g. linear or polynomial)
z <- reg.df[Index == "Shannon"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
           distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 0.018
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.33
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.34
anova(lm0,lm1)
anova(lm1,lm2) # preferred model is lm2, higher R2 and lowered RSS
# linear model itself is very bad
anova(lm2,lm3) 
rm(z)


z <- reg.df[Index == "Pielou"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
           distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 0.05
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.53
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.61
anova(lm0,lm1) 
anova(lm1,lm2)
anova(lm2,lm3) # preferred model is lm2, higher R2 and lowered RSS
# poly 3 has higher R2 but is not significant
rm(z)

lin.ls <- dlply(reg.df, .(Index), function(z){
  setDT(z)
  means <- z[, .(Diversity = mean(Diversity, na.rm = T),
                 distance = mean(distance, na.rm = T)), by = .(sample.type.year)]
  # both are poly 2
  lin <-  lm(means$Diversity ~ poly(means$distance, 2))

  # check linear assumptions
  #plot(lin) # normality not good
  # large sample sizes, normality does not affect results too much (central limit theorem)
  # homoscedasticity and independence important
  #summary(lin)
  #confint(lin, level = 0.95)
  
  # get data for plotting
  x <- data.frame(x = sort(means$distance))
  pred <- predict(lin, newdata = x, se = T)
  ci <- pred$se.fit[order(means$distance)] * qt(0.95 / 2 + 0.5, pred$df)
  y <- pred$fit[order(means$distance)]
  ymin <- y - ci
  ymax <- y + ci
  
  plot.df <- data.frame(x = sort(means$distance),
                        y = y,
                        ymin = ymin,
                        ymax = ymax,
                        se = pred$se.fit[order(means$distance)])
  
  # extract only colours that are in data frame
  colvec <- colvec[names(colvec) %in% as.character(levels(means$sample.type.year))]
  
  p <- ggplot() +
    theme_cust("pubr") +
    geom_point(data = z, aes(x = distance, y = Diversity), 
               colour = "gray40", alpha = 0.3, size = 2) +
    geom_line(data = plot.df, aes(x = x, y = y), inherit.aes = F, size = 1) +
    geom_line(data = plot.df, aes(x = x, y = ymax), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_line(data = plot.df, aes(x = x, y = ymin), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_point(data = means,
               aes(x = distance, y = Diversity, fill = sample.type.year), shape = 21, size = 3) +
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
       fitted = plot.df,
       plot = p)
  
})

(p <- ggarrange(lin.ls[[1]]$plot, lin.ls[[2]]$plot,
                ncol = 2, common.legend = T, legend = "right"))

(p <- annotate_figure(p, 
                      bottom = text_grob("Distance in ordination space",
                                         just = "centre")))
ggsave("./Figures/Final/Richness.RNA_distance_nonlin_reg.png", p,
       width = 22, height = 11, unit = "cm")
