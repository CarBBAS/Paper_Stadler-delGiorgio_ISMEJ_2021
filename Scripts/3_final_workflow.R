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

# Set seed for session and reproducibility of permutations
set.seed(3)

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
#asv.tab <- asv.tab[, 1:1000]

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
#tax.tab <- tax.tab[,1:1000]

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

#pb <- prune_taxa(!taxa_sums(pb) == 0, pb)
pb <- prune_samples(!sample_sums(pb) == 0, pb)

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
pb.mat <- log2(pb.mat + 1)
# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray) # 372 registers
# plot with custom function
dna.pcoa <- plot_bray(pb.bray.pcoa, .id = "DNA", colours = colvec, output = T)

p <- dna.pcoa$plot + guides(alpha = "none")

(zoom.terr <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
    theme_cust() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(size = 2.5, alpha = 0.2, fill = "gray20") +
    new_scale_fill() +
    geom_point(data = dna.pcoa$df[dna.pcoa$df$sample.type.year == "Soil" | 
                                    dna.pcoa$df$sample.type.year == "Soilwater" |
                                    dna.pcoa$df$sample.type.year == "Sediment",], 
               aes(x = Axis.1, y = Axis.2, fill = Season, shape = sample.type.year), size = 3) +
    scale_fill_viridis_d(option = "cividis", name = "Season") +
    scale_shape_manual(values = c(21,23, 25), "Habitat type") +
    coord_cartesian(ylim = c(-0.15, 0.1), xlim = c(-0.42, -0.22), expand = F) + 
    labs(x = paste("PC1"), 
         y = paste("PC2")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
)


(zoom.est <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
  theme_cust() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
  geom_point(size = 2.5, alpha = 0.2, fill = "gray20") +
  new_scale_fill() +
  geom_point(data = dna.pcoa$df[dna.pcoa$df$sample.type.year == "Estuary",], 
             aes(x = Axis.1, y = Axis.2, fill = abs(distance.from.mouth), shape = Season), size = 3) +
  scale_fill_viridis_c(name = "Distance from \nmouth (km)", direction = -1) +
  scale_shape_manual(values = c(21,23,25)) +
  scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
  coord_cartesian(ylim = c(-0.15, 0.25), xlim = c(-0.05, 0.3), expand = F) + # ensure aspect ratio
  labs(x = paste("PC1"), 
       y = paste("PC2")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
  #       alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
  #       fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
)
(zoom.river <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
  theme_cust() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
  geom_point(size = 2.5, alpha = 0.2, shape = 21, fill = "gray20") +
  #scale_fill_manual(values = colvec, name = "Habitat Type") +
  new_scale_fill() +
  geom_point(data = dna.pcoa$df[dna.pcoa$df$sample.type.year == "Downriver" | dna.pcoa$df$sample.type.year == "Upriver",], 
             aes(x = Axis.1, y = Axis.2, fill = abs(distance.from.mouth), shape = sample.type.year), size = 3) +
  scale_fill_continuous(type = "viridis", name = "Distance from \nmouth (km)", direction = -1) +
  #scale_fill_viridis_b(name = "Distance from mouth", direction = -1) +
  scale_shape_manual(values = c(21,23), name = "Habitat Type") +
  coord_cartesian(xlim = c(0, 0.32), ylim = c(-0.18,0.3), expand = F) + # ensure aspect ratio
  labs(x = paste("PC1"), 
       y = paste("PC2")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
  #       alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
  #       fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
)

# annotate collage boxes within main plot
p <- p + 
  annotate(geom = "rect", xmin = -0.42, ymin = - 0.15, 
           xmax = -0.21, ymax = 0.1, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "rect", xmin = 0, xmax = 0.32,
           ymin = - 0.13, ymax = 0.3, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "rect", xmin = -0.04, xmax = 0.28,
           ymin = - 0.12, ymax = 0.22, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "text", x = -0.405, y = 0.085, label = "b", alpha = 0.8, colour = "grey50") +
  annotate(geom = "text", x = 0.015, y = 0.285, label = "c", alpha = 0.8, colour = "grey50") +
  annotate(geom = "text", x = -0.025, y = 0.2, label = "d", alpha = 0.8, colour = "grey50")

collage <- ggarrange(p, 
          ggarrange(zoom.terr, zoom.river, zoom.est, ncol = 3, labels = c("b","c","d"), align = "hv"),
          nrow = 2, labels = c("a"), heights = c(0.6, 0.4))

# save
ggsave(paste0("./Figures/Final/PCoA_log_DNA_SampleType.tiff"), p,
       width = 12, height = 10, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_log_DNA_SampleType.png"),  p,
       width = 12, height = 10, unit = "cm")

ggsave(paste0("./Figures/Final/PCoA_log_DNA_collage.tiff"), collage,
       width = 25, height = 15, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_log_DNA_collage.png"),  collage,
       width = 25, height = 15, unit = "cm")


# Run PERMANOVA
# sensitive towards unbalanced sampling designs = bias.adjust
ord.df <- dna.pcoa[["df"]]
ord.df$groups <- paste(ord.df$sample.type.year, ord.df$Season, sep = "_")

setDT(ord.df)
# Test for significant difference between factors
perm.mod <- adonis(pb.mat ~ sample.type.year + Season, 
                   permutations = 9999, data = ord.df, sqrt.dist = T, method = "bray",
                   parallel = cl)
# habitat = F[12] = 17.09, R^2 = 0.35107, p = 1e-04
# season = F[2] = 11.08, R^2 = 0.03793, p = 1e-04
# iterations = 9999

# calculate multivariate dispersions
mod <- betadisper(pb.bray, 
                  group = as.character(ord.df$sample.type.year),
                  bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
mod2 <- betadisper(pb.bray, 
                  group = as.character(ord.df$Season),
                  bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
mod3 <- betadisper(pb.bray, 
                   group = ord.df$groups,
                   bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed

## Perform test
perm1 <- anova(mod, permutations = 9999, pairwise = T, parallel = cl) 
perm2 <- anova(mod2, permutations = 9999, pairwise = T, parallel = cl) 
perm3 <- anova(mod3, permutations = 9999, pairwise = T, parallel = cl) 
# runs into errors with permutest, anova.cca is same function but without errors...

# habitat alone = F[12] = 32.17, p = 1e-04
# season alone = F[2] = 42.79, p = 1e-04
# combined = F[27] = 17.77, p = < 2.2e-16

############
# Only RNA #
############
# subset only RNA samples
rna <- subset_samples(pb, DnaType == "RNA")

# extract ASV matrix
pb.mat <- otu_mat(rna)
pb.mat <- log2(pb.mat + 1)
# PCoA with Bray-Curtis
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

####################
#---------------------------------------------------------------------------------------------------#
####################
# Both DNA and RNA #
####################
# extract species table with species in columns
pb.mat <- otu_mat(pb)
# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray)
# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray) # 572 registers

#coor.dist <- as.matrix(dist(pb.bray.pcoa$vectors)) # same as original pb.bray matrix

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
# Run PERMANOVA
# sensitive towards unbalanced sampling designs = bias.adjust
ord.df <- all.pcoa[["df"]]
ord.df$groups <- paste(ord.df$sample.type.year, ord.df$Season, ord.df$DnaType, sep = "_")

setDT(ord.df)
# Test for significant difference between factors
perm.mod <- adonis(pb.mat ~ sample.type.year + Season + DnaType, 
                   permutations = 9999, data = ord.df, sqrt.dist = T, method = "bray",
                   parallel = cl)
# habitat = F[12] = 20.65, R^2 = 0.28824, p = 1e-04
# season = F[2] = 14.862, R^2 = 0.03458, p = 1e-04
# nucleic acid type = F[1] = 26.072, R^2 = 0.03033, p = 1e-04
# iterations = 9999

# calculate multivariate dispersions
mod <- betadisper(pb.bray, 
                  group = ord.df$sample.type.year,
                  bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
mod2 <- betadisper(pb.bray, 
                   group = ord.df$Season,
                   bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
mod3 <- betadisper(pb.bray, 
                   group = ord.df$DnaType,
                   bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
mod4 <- betadisper(pb.bray, 
                   group = ord.df$groups,
                   bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed

## Perform test
perm1 <- anova(mod, permutations = 9999, pairwise = T, parallel = cl) 
perm2 <- anova(mod2, permutations = 9999, pairwise = T, parallel = cl) 
perm3 <- anova(mod3, permutations = 9999, pairwise = T, parallel = cl) 
perm4 <- anova(mod4, permutations = 9999, pairwise = T, parallel = cl) 
# runs into errors with permutest and parallel, anova.cca is same function but without errors...

# habitat alone = F[12] = 31.577, p = 1e-04
# season alone = F[2] = 44.262, p = 1e-04
# nucleic acid type alone = F[1] = 1.6476, p = 0.1998
# combined = F[48] = 10.96, p = < 2.2e-16


####################
#- 3D interactive -#
####################
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
# convert distance matrix into long format
dissim.df <- melt.dist(pb.bray)

# Get meta data to rename DNA and RNA data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(DR.names, DnaType, Year, Season, sample.type.year), 
                   stringsAsFactors = F)

# add meta data for both sample x and sample y
dissim.df <- merge(dissim.df, meta, by.x =  "Sample.x", by.y = "Sample")
dissim.df <- merge(dissim.df, meta %>% select(DnaType, Sample, Year, DR.names), by.x =  "Sample.y", by.y = "Sample")

# omit all samples of 2015 (no RNA was taken, and sample name strategy changed -> creates duplicates)
dissim.df <- dissim.df[dissim.df$Year.x != 2015,]
dissim.df <- dissim.df[dissim.df$Year.y != 2015,]

# keep all rows where Sample.x and Sample.y are the same
dissim.df <- dissim.df[dissim.df$DR.names.x == dissim.df$DR.names.y,]
dissim.df <- dissim.df[!(dissim.df$DnaType.x == dissim.df$DnaType.y),] # omit all distances between same DnaType

dissim.dr <- dissim.df %>% select(ID = DR.names.x, Year = Year.x, Season, sample.type.year, dist)
setDT(dissim.dr)

# add new column to split plot into main and side panel
dissim.dr[, panels := "main"]
dissim.dr[sample.type.year == "Tributary" |
          sample.type.year == "Lake" |
          sample.type.year == "HeadwaterLakes" |
          sample.type.year == "Sediment", panels := "side"]

# calculate confidence interval and means of sample type and season combinations
sum.dissim <- dissim.dr[, .(mean =  mean(dist, na.rm = T),
                       conf.int = conf.int(dist),
                       stdev = sd(dist, na.rm = T)), by = .(sample.type.year, Season, panels)]

#-------------------------------------------------------------#
# Calculate distance between DNA and RNA points in PCoA space #
#-------------------------------------------------------------#
# use custom function to correct a few wrong sample names and match DNA-RNA counterpart samples
# calculating distance between points in two-dimensional space for both Bray-Curtis and Jaccard
# Use three distances, as variances explained by second and third axes are almost identical

pb.scores <- data.frame(Sample = as.character(row.names(pb.bray.pcoa$vectors)),
                        pb.bray.pcoa$vectors, stringsAsFactors = F)  # get first three axes
# merge with a selection of meta data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, 
                                                   DnaType, distance.from.mouth, DR.names), 
                   stringsAsFactors = F)
data <- merge(pb.scores, meta, by = "Sample")
data$Sample <- as.character(data$Sample)

setDT(data)
# melt datatable
temp <- melt(data, id.vars = c("DR.names","DnaType"), measure.vars = patterns("^Axis."),
             variable.name = "Axis", value.name = "Coordinates")
temp <- dcast(temp, DR.names + Axis ~ DnaType, value.var = c("Coordinates"))
# remove NAs
temp <- na.omit(temp)

# Calculate distance of temp axes
temp[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
temp <- temp[, .(sum.dist = sum(pnt.dist)), by = .(DR.names)] # sum temp axes
temp <- temp[, dist := sqrt(sum.dist)]

# Calculate distance for a selection of axes
temp.1d <- melt(data, id.vars = c("DR.names","DnaType"), measure.vars = patterns("^Axis."),
             variable.name = "Axis", value.name = "Coordinates")

temp.1d <- dcast(temp.1d[Axis == "Axis.3" | Axis == "Axis.2" | Axis == "Axis.1",], 
                 DR.names + Axis ~ DnaType, value.var = "Coordinates")
temp.1d <- temp.1d[, .(dist = abs(DNA - RNA)), by = .(DR.names, Axis)]

#temp.2d[, distance.12 := sqrt(abs((DNA_Axis.1 - RNA_Axis.1))^2 + abs((DNA_Axis.2 - RNA_Axis.2))^2)]
#temp.2d[, distance.13 := sqrt(abs((DNA_Axis.1 - RNA_Axis.1))^2 + abs((DNA_Axis.3 - RNA_Axis.3))^2)]

# combine back with categories
dist.dr <- temp[data, c("sample.type.year",
                        "Year", "Season") := list(i.sample.type.year,
                                                  i.Year, i.Season), on = .(DR.names)]
dist.1d <- temp.1d[data, c("sample.type.year",
                           "Year", "Season") := list(i.sample.type.year,
                                                     i.Year, i.Season), on = .(DR.names)]

# calculate confidence interval and means of sample type and season combinations
sum.dist <- dist.dr[, .(mean =  mean(dist, na.rm = T),
                       conf.int = conf.int(dist),
                       stdev = sd(dist, na.rm = T)),
                   by = .(sample.type.year, Season)]

sum.dist1d <- dist.1d[, .(mean =  mean(dist, na.rm = T),
                       conf.int = conf.int(dist),
                       stdev = sd(dist, na.rm = T)),
                   by = .(sample.type.year, Season)]

# See if distance and dissimilarity are correlated
# All axes
all.ax <- merge(dissim.dr, dist.dr, by.x = "ID", by.y = "DR.names")

(p <- ggplot(all.ax, aes(x = dist.x, y = dist.y, fill = sample.type.year.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = colvec) +
  labs(x = "Pair-wise Bray-Curtis dissimilarity", y = "Pair-wise distance in ordination across all axes"))

ggsave("./Figures/General/allaxesdist_dissim_cor.png", p)

# Only axis of interest
sec.ax <- merge(dissim.dr, dist.1d[Axis == "Axis.2",], by.x = "ID", by.y = "DR.names")
sec.ax <- merge(sec.ax, dist.dr, by.x = "ID", by.y = "DR.names")
sec.ax[, rel.dist := dist - dist.y]

(p <- ggplot(sec.ax, aes(x = dist.x, y = dist.y, fill = sample.type.year.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = colvec) +
  labs(x = "Pair-wise Bray-Curtis dissimilarity", y = "Pair-wise distance in ordination (PC2)")
)

ggsave("./Figures/General/allaxesdist_dissim_cor.png", p)
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

ggsave(paste0("./Figures/Final/DNARNA_dist.23.tiff"), p,
       width = 17, height = 10, unit = "cm")
ggsave(paste0("./Figures/Final/DNARNA_dist.23.png"),  p,
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
# 1. Calculate the mean abundance of each ASV for a sample type (e.g. reservoir, lake, stream, soil etc)
# 2. Define abundance groups based on rank abundance curves:
#-- * For each sample type, we create a rank abundance curve
#-- * Log-transform and take the derivative of the curve
#-- * Use second derivative of log-curve to define where the curve starts to bend.

# calculate the mean abundance of each ASV for all sample types

rel.df <- select_newest("./Objects", "201520162017_css_otu99_")
rel.df <- readRDS(
  paste0("./Objects/", rel.df))

setDT(rel.df)
# calculate means by sample type
means <- rel.df[, .(mean.css = mean(css.reads, na.rm = T),
                    sd.css = sd(css.reads, na.rm = T)), by = .(sample.type.year, DnaType, OTU)]
# order the abundances to make ranks
means <- means[mean.css != 0,] # remove all 0 observations
means <- means[order(mean.css, decreasing = T)]
means[, rank.abun := 1:.N, by = .(DnaType, sample.type.year)]
means[, log.mean := log1p(mean.css), by = .(DnaType, sample.type.year)]

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
# Get derivatives for an example and plot for supplementary material
x <- means[sample.type.year == "Soilwater" & DnaType == "DNA"]

spl <- smooth.spline(x$rank.abun, x$log.mean, spar = 0.7)
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
ggsave("./Figures/Final/abundance_class_ex_otu99.png", p,
       width = 22, height = 9, unit = "cm")

rm(x, p, raw, logged, pred, sec, first, spl)

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#####################
# Apply to the data #
#####################

classif.thres <- ddply(means, .(DnaType, sample.type.year), function(x){
  spl <- smooth.spline(x$rank.abun, x$log.mean, spar = 0.7)
  #pred <- predict(spl)
  #first <- predict(spl, deriv = 1) # first derivative
  sec <- predict(spl, deriv = 2) # second derivative
  setDT(x)
  x[rank.abun <= x$rank.abun[localMaxima(sec$y)[1]], ab.group := "Abundant", by = .(OTU)]
  x[rank.abun > x$rank.abun[localMaxima(sec$y)[1]] &
      rank.abun <= x$rank.abun[localMaxima(sec$y)[2]], ab.group := "Medium", by = .(OTU)]
  x[rank.abun > x$rank.abun[localMaxima(sec$y)[2]], ab.group := "Rare", by = .(OTU)]
  
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
# (see previous script) we settle with:
# Abundant >= 47 css reads
# Medium < 47 & >= 5 css reads
# Rare < 5 css reads

means[mean.css >= 47, ab.group := 1] # Abundant
means[mean.css < 47 & mean.css >= 5, ab.group := 2] # Medium
means[mean.css < 5 & mean.css > 0, ab.group := 3] # Rare

# code to resuscitate absent rows
temp <- dcast(means, DnaType + OTU ~ sample.type.year, value.var = c("ab.group")) # simple wide first
temp[is.na(temp)] <- 0 # much faster to overwrite NAs this way
ag.class <- melt(temp, id.vars = c("DnaType", "OTU"),
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
  univ.abun <- unique(as.character(x[ab.group == 1, abundant.counts := .N, by = .(OTU)][abundant.counts == length(levels(factor(sample.type.year))),]$OTU))
  # 2. Universally medium
  univ.med <- unique(as.character(x[ab.group == 2, abundant.counts := .N, by = .(OTU)][abundant.counts == length(levels(factor(sample.type.year))),]$OTU))
  # 3. Universally rare (all present)
  #univ.rare <- as.character(x[ab.group == 3, rare.counts := .N, by = .(OTU)][rare.counts == length(levels(factor(sample.type.year))),]$OTU)
  # none
  
  # 3. Universally rare (not all present)
  univ.rare <- unique(as.character(x[ab.group == 3 | ab.group == 0, rare.counts := .N, by = .(OTU)][rare.counts == length(levels(factor(sample.type.year))),]$OTU))
  
  # 4. Cosmopolitan
  x[, domain.count := .N, by = .(OTU, ecosys.domain)]
  abun.dom <- unique(as.character(x[ab.group == 1, abundant.counts.by.dom := .N, by = .(ecosys.domain, OTU)][abundant.counts.by.dom == domain.count,]$OTU))
  rare <- as.character(unique(x[ab.group == 3 | ab.group == 0,]$OTU))
  cosmo <- as.character(unique(x[OTU %in% abun.dom & !(OTU %in% rare),]$OTU))
  
  # 5. Abundant in all of one ecosystem domain, but never abundant in other domains == Specialist
  abun.dom <- x[ab.group == 1 | ab.group == 2,
                abun.counts.by.dom := .N, by = .(ecosys.domain, OTU)][abun.counts.by.dom == domain.count, abun.domain := ecosys.domain, by = .(OTU)]
  spec.asvs <- as.character(unique(abun.dom[!is.na(abun.domain),]$OTU))
  abun.dom[OTU %in% spec.asvs & is.na(abun.domain), abun.domain := "No" ]
  abun.dom <- as.character(unique(abun.dom[OTU %in% spec.asvs & 
                                             abun.domain != ecosys.domain &
                                             (ab.group == 1 | ab.group == 2),
                                           abun.counts := .N, by = .(OTU, ecosys.domain)][OTU %in% spec.asvs & abun.domain != ecosys.domain, .(sum.abun.counts = sum(abun.counts, na.rm = T)), by = .(OTU)][sum.abun.counts == 0,]$OTU))
  
  # 6. Shifters
  rest <- c(univ.abun, univ.med, univ.rare, cosmo, abun.dom)
  remain <- x[OTU %in% unique(x$OTU)[!(unique(x$OTU) %in% rest)]]
  
  remain[,
         abundant.counts := nrow(.SD[ab.group == 2 | ab.group == 3]), by = .(OTU)]
  remain[, rare.counts := nrow(.SD[ab.group == 1]), by = .(OTU)]
  remain[, absent.counts := nrow(.SD[ab.group == 0]), by = .(OTU)]
  
  # present in all samples
  # Mostly abundant, sometimes rare
  present.abundant.shifter <- as.character(unique(remain[abundant.counts >= rare.counts & absent.counts == 0,]$OTU)) # 2534
  # Mostly rare, sometimes abundant
  present.rare.shifter <- as.character(unique(remain[absent.counts == 0 &
                                               rare.counts > abundant.counts,]$OTU))
  
  # Rare 
  abundant.shifter <- as.character(unique(remain[!(OTU %in% c(present.abundant.shifter, present.rare.shifter)) &
                                                             abundant.counts >= (rare.counts + absent.counts),]$OTU)) # 1175
  # Mostly rare, sometimes abundant
  rare.shifter <- as.character(unique(remain[!(OTU %in% c(present.abundant.shifter, present.rare.shifter)) &
                                                       (rare.counts + absent.counts) > abundant.counts,]$OTU))
  
  #sanity check
  #length(univ.abun) + length(univ.med) + length(univ.rare) + length(cosmo) + length(abun.dom) + length(abundant.shifter) + length(rare.shifter) +  length(present.abundant.shifter) + length(present.rare.shifter) == length(unique(x$OTU))
  
  list(universal.abundant = univ.abun,
       universal.medium = univ.med,
       universal.rare = univ.rare,
       present.as = present.abundant.shifter,
       present.rs = present.rare.shifter,
       cosmopolitan = cosmo,
       specialist = abun.dom,
       abundant.shifter = abundant.shifter,
       rare.shifter = rare.shifter
       )
}, .parallel = T)

# Only extract abundance classification based on DNA
dna.ab.group <- abun.list[["DNA"]]

# only keep those abundance groups that have ASVs assigned
ls.asvs <- dna.ab.group[which(sapply(dna.ab.group, length) != 0L)]
grp.names <- names(ls.asvs)

# extract raw DNA - RNA relationship
rel.df <- rel.df[Year != 2015,]

# create new ID column with sample type and DNA Type
rel.df[, ID := paste(sample.type.year, DnaType, sep = "_")]

# what's the N ?
lapply(ls.asvs, length)
# universal.rare 17459
# present abundant shifter 36
# present rare shifter 7
# specialist 144
# abundant.shifter 1175
# rare.shifter 1364

#---------------------------------------------------------------------------------------------------#

#- Figure 4 -#

# Q2: Where in the rank abundance curve does community reshuffling occur?

# fit species onto PCoA ordination
pb.species <-
  data.frame(
    OTU = as.character(colnames(pb.mat)),
    wascores(pb.bray.pcoa$vectors[, 1:3],
             w = pb.mat),
    stringsAsFactors = F
  )

pb.species <- na.omit(pb.species)
setDT(pb.species)
# add abundance group names into species dataframe
for(i in 1:length(grp.names)){
  pb.species[OTU %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
  pb.species[OTU %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
}

# Correlate species to axes
cor.mat <- cor(all.pcoa$df[,2:4], pb.mat, method = "pearson")
cor.mat <- setDT(as.data.frame(cor.mat), keep.rownames = "Axes")
cor.df <- melt.data.table(
  cor.mat,
  id.vars = "Axes",
  measure.vars = patterns("OTU_"),
  variable.name = "OTU",
  value.name = "pear.cor"
)

# does not work on private computer
sim <- with(sample_df(pb), vegan::simper(t(otu_mat(pb)), DnaType, permutations = 1, parallel = cl))

sim <- big.simper(otu_mat(pb), sample_df(pb)$DnaType, permutations = 100, parallel = 2)

cor.df <- dcast(cor.df, OTU ~ Axes, value.var = "pear.cor")
colnames(cor.df)[2:4] <- c("cor.pc1","cor.pc2","cor.pc3")

cor.df <- cor.df[pb.species, c("cordi.pc1","cordi.pc2","cordi.pc3",
                     "ab.groups") := list(i.Axis.1, i.Axis.2, i.Axis.3,
                                          i.ab.names), on = .(OTU)]


ggplot(cor.df, aes(x = cor.pc1, y = cor.pc2, fill = ab.groups)) +
  theme_bw() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
  geom_point(data = subset(cor.df, ab.groups == "universal.rare"), shape = 21, alpha = 0.5) +
  geom_point(data = subset(cor.df, ab.groups != "universal.rare"), shape = 21) +
  #geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"), shape = 21) +
  #geom_point(data = subset(cor.df, ab.groups == "rare.shifter"), shape = 21) +
  #geom_point(data = subset(cor.df, ab.groups == "specialist"), shape = 21) +
  scale_fill_viridis_d() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

####################################################################
# compare slopes

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("","Terrestrial"
                   ,"","Water"),
  hjustvar = c(0,-0.08,1,1.2) ,
  vjustvar = c(0,1.2,0,1.2))

#hjust: lower values right, larger values left
#vjust lower values up, higher values down

#Bottom Left (h0,v0)","Top Left (h0,v1)"
#,"Bottom Right h1,v0","Top Right h1,v1"

(p1 <- ggplot() +
  theme_bw() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
  geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
             aes(x = cordi.pc1, y = abs(cor.pc1), colour = ab.groups),
             shape = 21, alpha = 0.5, colour = "grey80") +
  geom_point(data = subset(cor.df, ab.groups != "universal.rare"),
             aes(x = cordi.pc1, y = abs(cor.pc1), colour = ab.groups), shape = 21, alpha = 0.7) +
  #geom_smooth(data = cor.df[ab.groups != "present.rs",], 
  #              aes(x = cordi.pc1, y = abs(cor.pc1), group = ab.groups,
  #                   colour = ab.groups), size = 1.5,
  #              method = 'gam', formula = y ~ s(x, bs = "cs"), se = F) +
    geom_smooth(data = cor.df[cordi.pc1 > 0,], 
                aes(x = cordi.pc1, y = abs(cor.pc1), group = ab.groups,
                    colour = ab.groups), size = 1.5,
                method = 'lm', formula = y ~ x, se = F) +
    geom_smooth(data = cor.df[cordi.pc1 < 0,], 
                aes(x = cordi.pc1, y = abs(cor.pc1), group = ab.groups,
                    colour = ab.groups), size = 1.5,
                method = 'lm', formula = y ~ x, se = F) +
  geom_text(data=annotations,aes(x=xpos,y=ypos,
                                 hjust=hjustvar,vjust=vjustvar,
                                 label=annotateText)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_viridis_d(option = "viridis") +
    scale_colour_viridis_d(option = "viridis") +
  labs(x = "Species coordinates [PC1]", 
       y = expression(paste(" | Pearsons's ", italic("r")," | of OTUs with PC1"))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
)


annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("","RNA"
                   ,"","DNA"),
  hjustvar = c(0,-0.08,1,1.2) ,
  vjustvar = c(0,1.2,0,1.2))

(p2 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cordi.pc2, y = abs(cor.pc2), colour = ab.groups),
               shape = 21, alpha = 0.5, colour = "grey80") +
    geom_point(data = subset(cor.df, ab.groups != "universal.rare"),
               aes(x = cordi.pc2, y = abs(cor.pc2), colour = ab.groups), shape = 21, alpha = 0.7) +
    #geom_smooth(data = cor.df[ab.groups != "present.rs",], 
    #              aes(x = cordi.pc2, y = abs(cor.pc2), group = ab.groups,
    #                   colour = ab.groups), size = 1.5,
    #              method = 'gam', formula = y ~ s(x, bs = "cs"), se = F) +
    geom_smooth(data = cor.df[cordi.pc2 > 0,], 
                aes(x = cordi.pc2, y = abs(cor.pc2), group = ab.groups,
                    colour = ab.groups), size = 1.5,
                method = 'lm', formula = y ~ x, se = F) +
    geom_smooth(data = cor.df[cordi.pc2 < 0,], 
                aes(x = cordi.pc2, y = abs(cor.pc2), group = ab.groups,
                    colour = ab.groups), size = 1.5,
                method = 'lm', formula = y ~ x, se = F) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_viridis_d(option = "viridis") +
    scale_colour_viridis_d(option = "viridis") +
    labs(x = "Species coordinates [PC2]", 
         y = expression(paste(" | Pearsons's ", italic("r")," | of OTUs with PC2"))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("","Summer-Autumn"
                   ,"","Spring"),
  hjustvar = c(0,-0.08,1,1.2) ,
  vjustvar = c(0,1.2,0,1.2))

(p3 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cordi.pc3, y = abs(cor.pc3), colour = ab.groups),
               shape = 21, alpha = 0.5, colour = "grey80") +
    geom_point(data = subset(cor.df, ab.groups != "universal.rare"),
               aes(x = cordi.pc3, y = abs(cor.pc3), colour = ab.groups), shape = 21, alpha = 0.7) +
    #geom_smooth(data = cor.df[ab.groups != "present.rs",], 
    #              aes(x = cordi.pc3, y = abs(cor.pc3), group = ab.groups,
    #                   colour = ab.groups), size = 1.5,
    #              method = 'gam', formula = y ~ s(x, bs = "cs"), se = F) +
    geom_smooth(data = cor.df[cordi.pc3 > 0,], 
                aes(x = cordi.pc3, y = abs(cor.pc3), group = ab.groups,
                    colour = ab.groups), size = 1.5,
                method = 'lm', formula = y ~ x, se = F) +
    geom_smooth(data = cor.df[cordi.pc3 < 0,], 
                aes(x = cordi.pc3, y = abs(cor.pc3), group = ab.groups,
                    colour = ab.groups), size = 1.5,
                method = 'lm', formula = y ~ x, se = F) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_viridis_d(option = "viridis") +
    scale_colour_viridis_d(option = "viridis") +
    labs(x = "Species coordinates [PC3]", 
         y = expression(paste(" | Pearsons's ", italic("r")," | of OTUs with PC3"))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

# combine
ggarrange(p1, p2, p3, ncol = 3, common.legend = T)


#############################################################################

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
ggsave

##
set.seed(3)
t <- sample(1:nrow(pb.mat),10)
submat <- pb.mat[t,]

submat <- submat[,colSums(submat) > 0]

meta <- data.frame(Sample = as.character(row.names(sample_df(dna))),
                   sample_df(dna) %>% dplyr::select(sample.type.year, Season, Year, DnaType), 
                   stringsAsFactors = F)
meta <- meta[t,]

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
dna.sp <- mvabund(submat)

mod <- manyglm(dna.sp ~ meta$sample.type.year * meta$Season, family = "negative.binomial")
# warning but is integer
saveRDS(mod, "./Objects/manyglm.dna.negbinom.log.rds")

# check residuals, it's not optimal, but compared to other families, there is less of a pattern
png(filename="./Figures/General/manyglm_dna_residuals_binom_log.png")
plot(mod)
dev.off()

# test for habitat type and season effect
anova.mod <- anova(mod)
saveRDS(anova.mod, "./Objects/manyglm.dna.negbinom.anova.rds")
print("DONE")
pb.mat <- decostand(pb.mat, "hellinger")
pb.mori <- vegdist(pb.mat, method = "horn")
is.euclid(pb.mori) # FALSE
pb.mori <- sqrt(pb.mori) # make Euclidean
is.euclid(pb.mori) # TRUE
anova <- anova(mod)

## Run bayesian ordination
# test control options, for quick building. Not final
mcmc.control. <- list(n.burnin = 10, 
                      n.iteration = 400, 
                      n.thin = 30, 
                      seed = 3)

fit.lvmbinom <- boral(y = pb.mat, 
                      family = "negative.binomial", 
                      num.lv = 2, 
                      mcmc.control = mcmc.control.,
                      row.eff = "fixed")

####################################################################

####################
# Only terrestrial #
####################
ter <- subset_samples(pb, sample.type.year == "Soil" |
                        sample.type.year == "Sediment" |
                        sample.type.year == "Soilwater" |
                        sample.type.year == "Groundwater")

# remove ASVs that do not appear in this dataset
ter <- prune_taxa(taxa_sums(ter) != 0, ter)

pb.mat <- otu_mat(ter)
pb.mat <- log2(pb.mat + 1)
# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # TRUE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)
# plot with custom function
dna.pcoa <- plot_bray(pb.bray.pcoa, .id = "DNA", colours = colvec, output = T)

ggsave("./Figures/General/terr_PCoA.png", dna.pcoa$plot)

dissim.dr <- dissim.dnarna(ter, save.name = "All", output = T)

dist.dr <- dist.dnarna(dna.pcoa[["df"]], save.name = "terr", output = T)
#dist.dr <- dist.dnarna(dnarna.bray[["df"]], save.name = "All_2D", dimensions = 2)

test <- merge(dissim.dr$original.df, dist.dr$indiv.df, by.x = "ID", by.y = "DR.names")
#dist.dr$df[Axis == "Axis.2",]
ggplot(test, aes(x = dist, y = distance.1D, fill = sample.type.year.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = colvec) +
  labs(x = "Pair-wise Bray-Curtis dissimilarity", y = "Pair-wise distance in ordination")

p <- ggplot(test, aes(x = dist, y = distance.1D, fill = Season.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  labs(x = "Pair-wise Bray-Curtis dissimilarity", y = "Pair-wise distance in ordination")

ggsave("./Figures/General/terr_distdissimcor.png", p)
cor.test(test$dist, test$distance.1D)
