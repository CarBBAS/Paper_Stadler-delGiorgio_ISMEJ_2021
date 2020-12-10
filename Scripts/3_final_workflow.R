###---------------------------------------------###
#-   Script for all plots in the publication:   - #
#-   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   - #
###---------------------------------------------###
set.seed(3)
# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("phyloseq", "plyr", "tidyverse", "data.table", # wrangling
              "ggpubr", "ggnewscale", "cowplot", "plotly", # plotting,
              "kableExtra",# "xlsx", # making tables for manuscript (LaTeX) and export as excel (for ISME)
              "doMC", # parallel computing
              "vegan", "ape", "ade4") # statistics

### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Note: For package "xlsx" Java has to be installed
# For Linux users install dependencies: default-jre, default-jdk


### Functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")
#source("./Functions/SET_framework.R")


### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# prepare for parallel processing (with 'parallel') for 'vegan' functions
numCores <- detectCores()
cl <- makeCluster(numCores, type = "FORK") # using forking

# Set seed for session and reproducibility of permutations
# (NMDS were done, but not part of main workflow anymore)
set.seed(3)


# 2. Read and prepare data ---------------------------------------------------------------

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

names(colvec) <- as.character(sample.factors) # assign names for later easy selection

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

# set theme for plotting
theme_set(theme_bw())



# 3. Analysis ----------------------------------------------------------------------------------

## Figure 1: DNA -----------------------------------------------------------------------------------
# Q: Is the microbial assemblage different between habitat types and seasons?

#--------------------------#
#- Multivariate analysis -#
#--------------------------#
# We are dealing with large environmental gradients, thus we expect a high proportion of zeros
# Therefore, we need to select an asymmetrical similarity distance coefficient
# If we're dealing with samples from fairly homogeneous environmental conditions (short envir. gradients)
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

# subset only DNA samples
dna <- subset_samples(pb, DnaType == "DNA")

# remove ASVs that do not appear in this dataset
dna <- prune_taxa(taxa_sums(dna) != 0, dna)
dna <- prune_samples(sample_sums(dna) != 0, dna) # remove samples with no reads

# extract ASV matrix
pb.mat <- otu_mat(dna) # convert phyloseq obj to matrix
pb.mat <- log2(pb.mat + 1) # log transform, add 1 = for zeros
# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE
# we need euclidean distance to compute the distance metrics later on

ncol(pb.mat) # OTUs
nrow(pb.mat) # Samples

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray) # 372 registers
# plot with custom function (= made to avoid repetitive code)
# custom function is in ./Functions/custom_fun.R
dna.pcoa <- plot_pcoa(pb.bray.pcoa, physeq = dna, colours = colvec, output = T)

# main PCoA
p <- dna.pcoa$plot + guides(alpha = "none")

# zoom into terrestrial part
# different colouring by Season, habitat type as shape
zoom.terr <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
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


# zoom into estuaries
# different colouring by distance from mouth, shape is Season
zoom.est <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
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

# zoom into rivers (upriver vs. downriver)
# different colouring by distance from mouth, shape is Season
zoom.river <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
  theme_cust() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
  geom_point(size = 2.5, alpha = 0.2, shape = 21, fill = "gray20") +
  #scale_fill_manual(values = colvec, name = "Habitat Type") +
  new_scale_fill() +
  geom_point(data = dna.pcoa$df[dna.pcoa$df$sample.type.year == "Downriver" | 
                                  dna.pcoa$df$sample.type.year == "Upriver",], 
             aes(x = Axis.1, y = Axis.2, 
                 fill = abs(distance.from.mouth), shape = sample.type.year), size = 3) +
  stat_ellipse(data = subset(dna.pcoa$df, (sample.type.year == "Upriver" | 
                                             sample.type.year == "Downriver")  &
                               Season == "Spring"),
                 aes(x = Axis.1, y = Axis.2, group = Year), 
               linetype = "dashed", colour = "grey50") +
  #stat_ellipse(data = subset(dna.pcoa$df, (sample.type.year == "Upriver" | sample.type.year == "Downriver")),
  #               aes(x = Axis.1, y = Axis.2, group = paste(Year,Season), linetype = as.character(Year))) +
  annotate(geom = "text", x = c(0.28,0.28), y = c(0.255,0.12), label = c("2015", "2016"), 
           size = 3, colour = "grey50") +
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


# annotate collage boxes within main plot
p <- p + 
  annotate(geom = "rect", xmin = -0.42, ymin = - 0.15, 
           xmax = -0.21, ymax = 0.1, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "rect", xmin = 0, xmax = 0.32,
           ymin = - 0.13, ymax = 0.3, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "rect", xmin = -0.04, xmax = 0.28,
           ymin = - 0.12, ymax = 0.22, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "text", x = c(-0.405,0.015,-0.025), 
           y = c(0.085,0.285,0.2), label = c("b","c","d"), alpha = 0.8, colour = "grey50")

# combine all plots into one
(collage <- ggarrange(p, 
          ggarrange(zoom.terr, zoom.river, zoom.est, ncol = 3, labels = c("b","c","d"), align = "hv"),
          nrow = 2, labels = c("a"), heights = c(0.6, 0.4))
)

# save
ggsave(paste0("./Figures/Final/PCoA_log_DNA_SampleType.tiff"), p,
       width = 12, height = 10, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_log_DNA_SampleType.png"),  p,
       width = 12, height = 10, unit = "cm")

ggsave(paste0("./Figures/Final/PCoA_log_DNA_collage.tiff"), collage,
       width = 25, height = 15, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_log_DNA_collage.png"),  collage,
       width = 25, height = 15, unit = "cm")


# PERMANOVA -------------------------------------------------------------------------------------
# sensitive towards unbalanced sampling designs = bias.adjust
# extract PCoA axes and meta data
ord.df <- dna.pcoa[["df"]]
# create groups that will be tested in PERMANOVA
ord.df$groups <- paste(ord.df$sample.type.year, ord.df$Season, sep = "_")

setDT(ord.df)
# Test for significant difference between factors
perm.mod <- adonis(pb.mat ~ sample.type.year + Season, 
                   permutations = 9999, data = ord.df, sqrt.dist = T, method = "bray",
                   parallel = cl)
#-- results
# habitat = F[12] = 17.07, R^2 = 0.35083, p = 1e-04
# season = F[2] = 11.09, R^2 = 0.03798, p = 1e-04
# iterations = 9999

# calculate multivariate dispersions
# only habitat type
mod <- betadisper(pb.bray, 
                  group = as.character(ord.df$sample.type.year),
                  bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
# only Seasons
mod2 <- betadisper(pb.bray, 
                  group = as.character(ord.df$Season),
                  bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
# both factors
mod3 <- betadisper(pb.bray, 
                   group = ord.df$groups,
                   bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed

## Perform test
perm1 <- anova(mod, permutations = 9999, pairwise = T, parallel = cl) 
perm2 <- anova(mod2, permutations = 9999, pairwise = T, parallel = cl)
perm3 <- anova(mod3, permutations = 9999, pairwise = T, parallel = cl) 
# runs into errors with permutest, anova.cca is same function but runs without errors...

#-- results
# habitat alone = F[12] = 32.32, p = < 2.2e-16
# season alone = F[2] = 42.77, p = < 2.2e-16
# combined = F[27] = 17.81, p = < 2.2e-16

# individual comparisons (exploring differences in dispersion among groups)
TukeyHSD(mod2)

## Figure 2: DNA-RNA -----------------------------------------------------------------------------------
# Q2: Are the DNA and RNA assemblages different?

# remove OTUs that only appear in RNA
pb <- prune_taxa(taxa_names(pb) %in% taxa_names(dna), pb)

# extract species table with species in columns
pb.mat <- otu_mat(pb)
# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray)
# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)


ncol(pb.mat) # 20182 OTUs
nrow(pb.mat) # 572 registers

ncol(pb.bray.pcoa$vectors) # 571 axes

# Sanity check
# back calculate BC dissimilarity matrix from PCoA
#coor.dist <- as.matrix(dist(pb.bray.pcoa$vectors)) # same as original pb.bray matrix

# plot with custom function
all.pcoa <- plot_pcoa(pb.bray.pcoa, physeq = pb, axes = c(1:3), colours = colvec, output = T)
# plot second and third axes
pcoa.23 <- plot_pcoa(pb.bray.pcoa, physeq = pb, plot.axes = c(3,2), colours = colvec, output = T)

pcoa.plot <- all.pcoa$plot + theme(legend.position = "left")

(p <- ggarrange(all.pcoa$plot, pcoa.23$plot, ncol = 2, common.legend = T, legend = "right",
               align = "hv", labels = "auto")
)

# save
ggsave(paste0("./Figures/Final/PCoA_all_SampleType.tiff"), p,
       width = 20, height = 11, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_all_SampleType.png"),  p,
       width = 20, height = 11, unit = "cm")

# check screeplot
ggplot(pb.bray.pcoa$values[pb.bray.pcoa$values$Eigenvalues > 1,],
       aes(x = as.numeric(row.names(pb.bray.pcoa$values[pb.bray.pcoa$values$Eigenvalues >= 1,])), y = Eigenvalues)) +
         geom_col()

# how many axes to have 75% of variance captured?
nrow(pb.bray.pcoa$values[pb.bray.pcoa$values$Cumul_eig <= 0.75,])
# 186 axes

# explore other axes
#for(i in 2:ncol(pb.bray.pcoa$vectors)){
#  temp <- plot_pcoa(pb.bray.pcoa, physeq = pb, plot.axes = c(1,i), colours = colvec, output = T)
#  ggsave(paste0("./Figures/PCoA/PCoA_DR_ax1_",i,".png"), temp$plot,
#         width = 12, height = 10, unit = "cm")
#}


# PERMANOVA -------------------------------------------------------------------------------------
# sensitive towards unbalanced sampling designs = bias.adjust
ord.df <- all.pcoa[["df"]]
ord.df$groups <- paste(ord.df$sample.type.year, ord.df$Season, ord.df$DnaType, sep = "_")

setDT(ord.df)
# Test for significant difference between factors
perm.mod <- adonis(pb.mat ~ sample.type.year + Season + DnaType, 
                   permutations = 9999, data = ord.df, sqrt.dist = T, method = "bray",
                   parallel = cl)
#-- results
# habitat = F[12] = 20.637, R^2 = 0.28819, p = 1e-04
# season = F[2] = 14.845, R^2 = 0.03455, p = 1e-04
# nucleic acid type = F[1] = 25.988, R^2 = 0.03024, p = 1e-04
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

#-- results
# habitat alone = F[12] = 31.611, p =  < 2.2e-16
# season alone = F[2] = 44.548, p = < 2.2e-16
# nucleic acid type alone = F[1] = 1.6727, p = 0.1964
# combined = F[48] = 10.921, p = < 2.2e-16

tuk <-TukeyHSD(mod4)
t <- data.frame(Groups = rownames(tuk$group),unlist(tuk$group))
t <- t %>% separate(Groups, into = c("habitat","season","nucacid",
                                "habitat_y","season_y","nucacid_y"), sep = "_|-")
row.names(t) <- NULL
# all pair-wise comparisons that are significant
t <- t[t$p.adj < 0.05,]

## 3D plot ---------------------------------------------------------------------------------------
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

## Figure 3: Distance -----------------------------------------------------------------------------------
# Q: How different are the the DNA and RNA assemblages of the same sample?

# Calculate incidence based dissimilarity
# PCoA with Sorensen, incidence based equivalent of Bray Curtis
pb.soren <- vegdist(pb.mat, method = "bray", binary = T)
is.euclid(pb.soren) # FALSE
pb.soren <- sqrt(pb.soren) # make Euclidean
is.euclid(pb.soren)
# make PCoA
pb.soren.pcoa <- ape::pcoa(pb.soren) 

soren.12 <- plot_pcoa(pb.soren.pcoa, physeq = pb, axes = c(1:2), colour = colvec, output = T)
# plot second and third axes
soren.32 <- plot_pcoa(pb.soren.pcoa, physeq = pb, plot.axes = c(3,2), colour = colvec, output = T)
# less strong, but axis 3 separates nucleic acid types

# how many axes to have 75% of variance captured?
nrow(pb.soren.pcoa$values[pb.soren.pcoa$values$Cumul_eig <= 0.75,])
# 236 axes

# plot Sorensen plots for supplementary
(p <- ggarrange(soren.12$plot, soren.32$plot, ncol = 2, common.legend = T, legend = "right",
               align = "hv", labels = "auto"))

# save
ggsave(paste0("./Figures/Final/PCoA_all_Sorensen.tiff"), p,
       width = 20, height = 11, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_all_Sorensen.png"),  p,
       width = 20, height = 11, unit = "cm")

#-------------------------------------------------------#
# Extract pair-wise dissimilarity among DNA-RNA samples #
#-------------------------------------------------------#
# convert distance matrix into long format
dissim.df <- rbind(melt.dist(pb.bray) %>% mutate(Metric = "Bray"),
                   melt.dist(pb.soren) %>% mutate(Metric = "Sorensen"))

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

dissim.dr <- dissim.df %>% select(Metric, ID = DR.names.x, Year = Year.x, Season, sample.type.year, dist)
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
                       stdev = sd(dist, na.rm = T)), by = .(Metric,sample.type.year, Season, panels)]

cast.dis <- dcast(dissim.dr, ID + sample.type.year + Season ~ Metric, value.var = "dist")


(pw.sorbray <- ggplot(cast.dis, aes(x =  Bray, y = Sorensen)) +
  theme_cust() +
  geom_point(aes(fill = sample.type.year, shape = Season), size = 2, alpha = 0.9) +
  scale_fill_manual(values =  colvec[names(colvec) %in% as.character(levels(cast.dis$sample.type.year))],
                    name = "Habitat Type") +
  scale_shape_manual(values = c(21, 23, 25)) +
  labs(x = expression(paste("DNA-RNA ", italic(D)["BC"])), #paste("Pairwise ", italic(D)["BC"]^"DNA:RNA")
       y = expression(paste("DNA-RNA ", italic(D)["S"]))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate(geom = "text", x = 0.53, y = 0.57, label = "1:1", angle = 45, size = 3) +
  theme(legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(0.4, "cm")) +
  guides(fill = guide_legend(override.aes=list(shape=21, size = 2), order = 1),
         shape = guide_legend(override.aes=list(size = 2)))
)

ggplot(cast.dis, aes(y = Sorensen, x =  Bray)) +
  theme_bw() +
  geom_point(aes(fill = Season), size = 3, shape = 21) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = expression(paste(italic(D)["BC"])), y = expression(paste(italic(D)["S"]))) +
  annotate(geom = "text", x = 0.55, y = 0.57, label = "1:1", angle = 45) +
  theme(legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes=list(shape=21), order = 1))

#-------------------------------------------------------------#
# Calculate distance between DNA and RNA points in PCoA space #
#-------------------------------------------------------------#
# use custom function to correct a few wrong sample names and match DNA-RNA counterpart samples
# calculating distance between points in two-dimensional space for both Bray-Curtis and Sorensen

pb.scores <- rbind(data.frame(Sample = as.character(row.names(pb.bray.pcoa$vectors)), Metric = "Bray",
                        pb.bray.pcoa$vectors, stringsAsFactors = F),
                   data.frame(Sample = as.character(row.names(pb.soren.pcoa$vectors)), Metric = "Sorensen",
                              pb.soren.pcoa$vectors, stringsAsFactors = F))# get first three axes
# merge with a selection of meta data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, 
                                                   DnaType, distance.from.mouth, DR.names), 
                   stringsAsFactors = F)
data <- merge(pb.scores, meta, by = "Sample")
data$Sample <- as.character(data$Sample)

setDT(data)
setcolorder(data, c("Sample","Metric",
                    "sample.type.year","Season","Year", "DnaType","distance.from.mouth","DR.names",
                    colnames(data)[!(colnames(data) %in% c("Sample","Metric",
                                                           "sample.type.year","Season","Year", "DnaType","distance.from.mouth","DR.names"))]))
# melt datatable
temp <- melt.data.table(data, id.vars = c("DR.names","DnaType", "Metric"), measure.vars = patterns("^Axis."),
             variable.name = "Axis", value.name = "Coordinates")
temp <- dcast(temp, DR.names + Axis + Metric ~ DnaType, value.var = c("Coordinates"))
# remove NAs
temp <- na.omit(temp)

# Calculate distance of all axes
temp[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
temp <- temp[, .(sum.dist = sum(pnt.dist)), by = .(Metric, DR.names)] # sum temp axes
temp <- temp[, dist := sqrt(sum.dist)]

# Calculate distance among axes capturing 75% of variation
# For Bray-Curtis
bray.df <- data.frame(Sample = as.character(row.names(pb.bray.pcoa$vectors)), Metric = "Bray",
                      pb.bray.pcoa$vectors, stringsAsFactors = F)
bray.df <- bray.df[,1:which(colnames(bray.df) == paste("Axis", 
                                                  nrow(pb.bray.pcoa$values[pb.bray.pcoa$values$Cumul_eig <= 0.75,]), sep = "."))]

# For Sorensen
soren.df <- data.frame(Sample = as.character(row.names(pb.soren.pcoa$vectors)), Metric = "Sorensen",
                       pb.soren.pcoa$vectors, stringsAsFactors = F)
soren.df <- soren.df[,1:which(colnames(soren.df) == paste("Axis", 
                                                        nrow(pb.soren.pcoa$values[pb.soren.pcoa$values$Cumul_eig <= 0.75,]), sep = "."))]

# merge with a selection of meta data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, 
                                                   DnaType, distance.from.mouth, DR.names), 
                   stringsAsFactors = F)
bray.df <- setDT(merge(meta, bray.df, by = "Sample"))
soren.df <- setDT(merge(meta, soren.df, by = "Sample"))

# melt and combine bray-curtis and sorensen results into one data frame
temp.75 <- rbind(melt.data.table(bray.df, id.vars = c("DR.names","DnaType", "Metric"), measure.vars = patterns("^Axis."),
                      variable.name = "Axis", value.name = "Coordinates"),
      melt.data.table(soren.df, id.vars = c("DR.names","DnaType", "Metric"), measure.vars = patterns("^Axis."),
                      variable.name = "Axis", value.name = "Coordinates"))

temp.75 <- dcast(temp.75, DR.names + Axis + Metric ~ DnaType, value.var = c("Coordinates"))
# remove NAs
temp.75 <- na.omit(temp.75)
temp.75[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
temp.75 <- temp.75[, .(sum.dist = sum(pnt.dist)), by = .(Metric, DR.names)] # sum temp axes
temp.75 <- temp.75[, dist := sqrt(sum.dist)] # take sqrt

# combine back with categories
# all axes
dist.dr <- temp[data, c("sample.type.year",
                        "Year", "Season") := list(i.sample.type.year,
                                                  i.Year, i.Season), on = .(DR.names)]
# 75% variance
dist.75 <- temp.75[data, c("sample.type.year",
                           "Year", "Season") := list(i.sample.type.year,
                                                     i.Year, i.Season), on = .(DR.names)]

# Make Bray to Sorensen ratio (only 75%)
diff.df <- dcast(dist.75, DR.names + sample.type.year + Season ~ Metric, value.var = "dist")
diff.df <- diff.df[, bcs.ratio := Sorensen / Bray]
sum.diff <- diff.df[, .(mean =  mean(bcs.ratio, na.rm = T),
                       conf.int = conf.int(bcs.ratio),
                       stdev = sd(bcs.ratio, na.rm = T),
                       n = .N),
                   by = .(sample.type.year, Season)]

# calculate confidence interval and means of sample type and season combinations
sum.dist <- dist.dr[, .(mean =  mean(dist, na.rm = T),
                       conf.int = conf.int(dist),
                       stdev = sd(dist, na.rm = T),
                       n = .N),
                   by = .(Metric, sample.type.year, Season)]

sum.dist75 <- dist.75[, .(mean =  mean(dist, na.rm = T),
                          conf.int = conf.int(dist),
                          stdev = sd(dist, na.rm = T),
                          n = .N),
                      by = .(Metric, sample.type.year, Season)]

# add new column to split plot into main and side panel
sum.ls <- lapply(list(sum.dist, sum.dist75, sum.diff), function(x){
  x[, panels := "main"]
  x[sample.type.year == "Tributary" |
               sample.type.year == "Lake" |
               sample.type.year == "Riverine \nLakes" |
               sample.type.year == "Sediment", panels := "side"]
  x[, sample.type.year := factor(sample.type.year, levels = c("Soil","Sediment",
                                                              "Soilwater",
                                                              "Stream", "Tributary",
                                                              "Riverine \nLakes", "Headwater \nPonds", "Lake",
                                                              "Upriver",
                                                              "Reservoirs","Downriver",
                                                              "Estuary"))]
})



# Proof that distance is extracting part of the individual pairwise dissimilarity
# If distance is calculated over all PCoA dimension, we back-calculate the pairwise dissimilarity
all.ax <- merge(dissim.dr[dissim.dr$Metric == "Bray",], 
                dist.dr[dist.dr$Metric == "Bray",], 
                by.x = "ID", by.y = "DR.names")

#all.ax <- merge(dissim.dr[dissim.dr$Metric == "Sorensen",], 
#                dist.dr[dist.dr$Metric == "Sorensen",], 
#                by.x = "ID", by.y = "DR.names") # works for both metrics

(p <- ggplot(all.ax, aes(x = dist.x, y = dist.y, fill = sample.type.year.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  labs(x = expression(paste("Pairwise ", italic("D")["BC"])),
       y = expression(paste(italic("m")["BC"]^"n = 100%")))
)

ggsave("./Figures/General/allaxesdist_dissim_cor.png", p)

# Only axis of interest
#sec.ax <- merge(dissim.dr, dist.1d[Axis == "Axis.2",], by.x = "ID", by.y = "DR.names")
#sec.ax <- merge(sec.ax, dist.dr, by.x = "ID", by.y = "DR.names")
#sec.ax[, rel.dist := dist - dist.y]

#(p <- ggplot(sec.ax, aes(x = dist.x, y = dist.y, fill = sample.type.year.x)) +
#  theme_bw() +
#  geom_point(shape = 21) +
#  scale_fill_manual(values = colvec) +
#  labs(x = "Pair-wise Bray-Curtis dissimilarity", y = "Pair-wise distance in ordination (PC2)")
#)

#ggsave("./Figures/General/allaxesdist_dissim_cor.png", p)
# dissimilarity and distance show different patterns....

# Scatter plots
plot.df <- sum.ls[[3]][!is.nan(sum.ls[[3]]$mean),]
d <- scatter_panels(plot.df, labs = c("Habitat Type", 
                                      expression(paste(italic("m")["S"], " : ", italic("m")["BC"]))))

(dm.plot <- ggarrange(d$main + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
                       
                        theme(axis.title.y = element_text(margin = margin(r = 10)),
                              axis.text.x = element_text(size = 7)),
                    d$side + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
                      theme(
                        axis.text.x = element_text(size = 7)),
                    widths = c(3,1),
                    ncol = 2, nrow = 1, 
                    #common.legend = T,
                    align = "h",
                    font.label = list(size = 10),
                    legend = "none")
)

# only mBC
bray.df <- sum.ls[[2]][!is.nan(sum.ls[[2]]$mean),]
bray.df <- bray.df[Metric == "Bray",]
lim.point <- max(bray.df$mean + bray.df$stdev, na.rm = T)
p <- scatter_panels(bray.df, labs = c("Habitat Type", 
                                      expression(paste(italic("m")["BC"]))))

# only mS
soren.df <- sum.ls[[2]][!is.nan(sum.ls[[2]]$mean),]
soren.df <- soren.df[Metric == "Sorensen",]
lim.point <- max(soren.df$mean + soren.df$stdev, na.rm = T)
s <- scatter_panels(soren.df, labs = c("Habitat Type", 
                                      expression(paste(italic("m")["S"]))))

left <- ggarrange(p$main + 
            geom_text(aes(y = lim.point + 0.03, colour = Season, label = n), position = position_dodge(width=0.7),
                      size = 2, show.legend = F) +
            lims(y = c(min(bray.df$mean - bray.df$stdev, na.rm = T), lim.point + 0.03)) +
            theme(axis.text.x = element_blank(), legend.position = "none"),
            s$main + 
              lims(y = c(min(soren.df$mean - soren.df$stdev, na.rm = T), lim.point + 0.03)) +
              theme(axis.text.x = element_blank(), legend.position = "none"),
          d$main + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
            theme(axis.title.y = element_text(margin = margin(r = 10)),
                  axis.text.x = element_text(size = 7), legend.position = "none"),
          nrow = 3, align = "v", heights = c(0.34,0.33,0.45), labels = "auto")

right <- ggarrange(p$side + geom_text(aes(y = lim.point + 0.03, colour = Season, label = n), position = position_dodge(width=0.7),
                                     size = 2, show.legend = F) +
                    lims(y = c(min(bray.df$mean - bray.df$stdev, na.rm = T), lim.point + 0.03)) +
                    theme(axis.text.x = element_blank()),
                   s$side +
                     lims(y = c(min(soren.df$mean - soren.df$stdev, na.rm = T), lim.point + 0.03)) +
                     theme(axis.text.x = element_blank()),
                  d$side + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
                    theme(
                      axis.text.x = element_text(size = 7)),
                  nrow = 3, align = "v", heights = c(0.34,0.33,0.45), common.legend = T, legend = "right")

(scat.plots <- ggarrange(left, right, ncol = 2, widths = c(2,1)))


(distdissm.col <- ggarrange(pw.sorbray, scat.plots, ncol = 2, #labels = c("a",""),
                            widths = c(1.5,2)))

ggsave("./Figures/Final/distBC_distratio_scat.png", scat.plots,
       width =15, height = 10, units = "cm") #width =15, height = 8.2,
ggsave("./Figures/Final/dissim_braysor.png", pw.sorbray,
       width =12, height = 8.2, units = "cm")

#ggsave("./Figures/Final/dissimBCS_distBC_distratio.png", distdissm.col,
#       width = 23, height = 8.2, units = "cm")
#ggsave("./Figures/Final/distBC_distratio_scatter.tiff", scat.plots,
#       width = 20, height = 9, units = "cm")

#ggsave("./Figures/Final/distBC_distratio_scatter.png", scat.plots,
#       width = 15, height = 9, units = "cm")
#ggsave("./Figures/Final/distBC_distratio_scatter.tiff", scat.plots,
#       width = 15, height = 9, units = "cm")
#expression(paste(Delta, " ", italic("m")["S"], " - ", italic("m")["BC"])))


# How do we interpret the mS:mBC ratio? --------------------------------------------------------
# extract sample pairs above 1 and below 1 of the mS:mBC ratio
big <- as.character(diff.df[bcs.ratio > 1,]$DR.names)
small <- as.character(diff.df[bcs.ratio < 1,]$DR.names)

# melt community matrix for tidy data set
commat <- melt.data.table(
  setDT(as.data.frame(otu_mat(pb)), keep.rownames = "Sample"),
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

# add meta variables
commat[setDT(sample_df(pb), keep.rownames = "Sample"), 
       c("DR.names", "DnaType") := list(i.DR.names, i.DnaType), on = .(Sample)]

# keep only the ones that have both DNA and RNA
#commat <- commat[duplicated(DR.names),]

# cast so that we can calculate difference between DNA and RNA of each sample
temp <- dcast(commat, DR.names + OTU ~ DnaType, value.var = c("reads"))
temp[, diff := DNA - RNA]
temp <- temp[!is.na(diff),]

# calculate summary variables
abs.diff <- temp[, .(mean.diff = mean(abs(diff), na.rm = T),
                     max.diff = max(abs(diff), na.rm = T),
                     var.diff = var(abs(diff), na.rm = T)), by = .(DR.names)]

# calculate how many taxa are shared by sample
n.shared <- temp[DNA > 0 & RNA > 0, ]
n.shared <- n.shared[, .(n = .N), by = .(DR.names)]

# calculate the abundance difference mean of the OTUs that have been classified as abundant
abun.diff <- temp[DNA >= 47,]
abun.diff <- abun.diff[, .(mean.abun.diff = mean(abs(diff))), by = .(DR.names)]

# combine all to one data set
all.met <- abs.diff[n.shared, n.shared := i.n, on = .(DR.names)]
all.met <- all.met[diff.df, ratio := i.bcs.ratio, on = .(DR.names)]
all.met <- all.met[abun.diff, mean.abun.diff := i.mean.abun.diff, on = .(DR.names)]

# give >1 and <1 ms:mBC categories
all.met[DR.names %in% big, ratio.cat := ">1"]
all.met[DR.names %in% small, ratio.cat := "<1"]
#all.met <- all.met[!is.na(ratio.cat),]

# quick plots
ggplot(all.met, aes(x = ratio.cat, y = mean.abun.diff)) +
  geom_boxplot()
ggplot(all.met, aes(x = ratio.cat, y = log(mean.diff))) +
  geom_boxplot()
ggplot(all.met, aes(x = ratio.cat, y = n.shared)) +
  geom_boxplot()

# Kruskal-Wallis -------------------------------------------------------------------------------------
# Statistically test whether there is a difference between these categories
# 1. Mean abundance ----------------------------------------------------------------------------------
# Is there a difference in abundance difference between the categories identified by the mS:mBC ratio?

# Check ANOVA assumptions
# Do we have outliers?
(outliers <- all.met %>%
    group_by(ratio.cat) %>%
    rstatix::identify_outliers(mean.diff) %>%
    select(DR.names, ratio.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]
#aov.df <- all.met
# Check normality
model <- lm(mean.diff ~ ratio.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis instead of one-way ANOVA
(stat.kw <- kruskal.test(mean.diff ~ ratio.cat, data = aov.df)) # significantly different (*** level, p < 0.0001)

# Difference in mean abundance difference among >1 and <1 mS:mBC is statistically significant

(m.ab <- ggplot(aov.df, aes(x = ratio.cat, y = mean.abun.diff)) +
  geom_boxplot(width = 0.3) +
  labs(y = expression(paste("Mean | ", Delta, " OTU CSS reads |")), x = expression(paste(m[S]:m[BC]))) +
  annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(750,730,730), yend = c(750,750,750)) +
  annotate(geom = "text", x = 1.5, y = 800, label = paste("Kruskal-Wallis, p ", abbrev.p(stat.kw$p.value)[1])))

# 2. Mean abundance of dominant taxa ------------------------------------------------------------------
# Is there a difference in abundance of dominant OTUs between the categories identified by the mS:mBC ratio?
# Do we have outliers?
(outliers <- all.met %>%
    group_by(ratio.cat) %>%
    rstatix::identify_outliers(mean.abun.diff) %>%
    select(DR.names, ratio.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

# Check normality
model <- lm(mean.abun.diff ~ ratio.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(mean.abun.diff ~ ratio.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(dom.diff <- ggplot(aov.df, aes(x = ratio.cat, y = mean.abun.diff)) +
    geom_boxplot(width = 0.3) +
    labs(y = expression(atop(paste("Mean | ", Delta, " OTU CSS reads |"),"of abundant taxa")), x = expression(paste(m[S]:m[BC]))) +
    annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(600,580,580), yend = c(600,600,600)) +
    annotate(geom = "text", x = 1.5, y = 640, label = paste("Kruskal-Wallis, p ", abbrev.p(stat.kw$p.value)[1])))


# 3. Number of shared OTUs ------------------------------------------------------------------
# Is there a difference in number of shared taxa between the categories identified by the mS:mBC ratio?
# Do we have outliers?
(outliers <- all.met %>%
    group_by(ratio.cat) %>%
    rstatix::identify_outliers(n.shared) %>%
    select(DR.names, ratio.cat, is.outlier, is.extreme))

# None are extreme, do not exclude
aov.df <- all.met

# Check normality
model <- lm(n.shared ~ ratio.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(n.shared ~ ratio.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(n.sh <- ggplot(aov.df, aes(x = ratio.cat, y = n.shared)) +
    geom_boxplot(width = 0.3) +
    labs(y = "Mean number of shared OTUs", x = expression(paste(m[S]:m[BC]))) +
    annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(770,750,750), yend = c(770,770,770)) +
    annotate(geom = "text", x = 1.5, y = 800, label = paste("Kruskal-Wallis, p ", abbrev.p(stat.kw$p.value)[1])))

(krusk.plots <- ggarrange(m.ab, dom.diff, n.sh, ncol = 3, align = "hv", labels = "auto"))

ggsave("./Figures/Final/ratio_krusk_tests.png", krusk.plots, width = 22, height = 10, units = "cm")



## Figure 4: Abundance groups-------------------------------------------------------------------
# Q: Who in the rank abundance curve is driving the dissimilarity patterns between DNA and RNA?

# Define abundance groups-----------------------------------------------------------------------

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

# smooth and get derivative
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


# Apply abundance groups -----------------------------------------------------------------------

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
#ag.class[, ecosys.domain := "Freshwater"]
#ag.class[sample.type.year == "Soil" | sample.type.year == "Soilwater" | sample.type.year == "Hyporheicwater"
#         | sample.type.year == "Wellwater" | sample.type.year == "Sediment",
#         ecosys.domain := "Soily"]
#ag.class[sample.type.year == "Marine", ecosys.domain := NA]

# Custom function to identify ASVs for each group

abun.list <- dlply(ag.class, .(DnaType), function(x){
  setDT(x)
  # 1. Universally abundant
  univ.abun <- unique(as.character(x[ab.group == 1, abun.counts := .N, by = .(OTU)][abun.counts == length(levels(factor(sample.type.year))),]$OTU))
  # 2. Universally moderate/medium
  univ.med <- unique(as.character(x[ab.group == 2, med.counts := .N, by = .(OTU)][med.counts == length(levels(factor(sample.type.year))),]$OTU))
  # 3. Universally rare (all present)
  univ.rare <- as.character(x[ab.group == 3, rare.counts := .N, by = .(OTU)][rare.counts == length(levels(factor(sample.type.year))),]$OTU)
  #remove cols
  x[, c("abun.counts", "med.counts", "rare.counts") := list(NULL, NULL, NULL)]
  
  # 4. Abundant (not all present)
  abundant <- unique(as.character(x[ab.group == 1 | ab.group == 0, abun.counts := .N, by = .(OTU)][abun.counts == length(levels(factor(sample.type.year))),]$OTU))
  # 5. Moderate (not all present)
  moderate <- unique(as.character(x[ab.group == 2 | ab.group == 0, med.counts := .N, by = .(OTU)][med.counts == length(levels(factor(sample.type.year))),]$OTU))
  # 6. Rare (not all present)
  rare <- unique(as.character(x[ab.group == 3 | ab.group == 0, rare.counts := .N, by = .(OTU)][rare.counts == length(levels(factor(sample.type.year))),]$OTU))
  
  # 6. Habitat specialist
  # Abundant in only one habitat, never abundant in other habitats
  hab.spec <- as.character(unique(x[ab.group == 1 | ab.group == 2,
                abun.counts.by.hab := nrow(.SD), by = .(OTU)][abun.counts.by.hab == 1,]$OTU))
  
  # remove OTUs classified as specialists from abundant, moderate and rare
  abundant <- abundant[!(abundant %in% hab.spec)]
  moderate <- moderate[!(moderate %in% hab.spec)]
  rare <- rare[!(rare %in% hab.spec)] # there should be none, but just to be safe
  
  #remove cols
  x[, c("abun.counts", "med.counts", "rare.counts", "abun.counts.by.hab") := list(NULL, NULL, NULL, NULL)]
  
  # 7. Shifters
  rest <- c(univ.abun, univ.med, univ.rare, abundant, moderate, rare, hab.spec)
  remain <- x[OTU %in% unique(x$OTU)[!(unique(x$OTU) %in% rest)]]
  
  remain[,
         abun.counts := nrow(.SD[ab.group == 3]), by = .(OTU)]
  remain[,med.counts := nrow(.SD[ab.group == 2 ]), by = .(OTU)]
  remain[, rare.counts := nrow(.SD[ab.group == 1]), by = .(OTU)]
  remain[, absent.counts := nrow(.SD[ab.group == 0]), by = .(OTU)]
  
  # present in all samples = cosmopolitan, shifts between medium, abundant, rare
  cosmopolitan <- as.character(unique(remain[absent.counts == 0,]$OTU))
  
  # Shifts between medium and abundant
  upper.shifter <- as.character(unique(remain[!(OTU %in% cosmopolitan) & 
                                                   med.counts > 0 & rare.counts == 0 & abun.counts > 0,]$OTU))
  # Shifts between rare and medium
  lower.shifter <- as.character(unique(remain[!(OTU %in% cosmopolitan) & 
                                                 med.counts > 0 & rare.counts > 0 & abun.counts == 0,]$OTU))
  
  # Shifter, shifts between rare, medium and abundant or between rare an abundant
  shifter <- as.character(unique(remain[!(OTU %in% c(cosmopolitan, upper.shifter, lower.shifter)) & 
                                          rare.counts > 0 & abun.counts > 0,]$OTU))
  
  #sanity check
  if(length(univ.abun) + length(univ.med) + length(univ.rare) + length(abundant) + length(moderate) + length(rare) +
     length(hab.spec) +
     length(cosmopolitan) + length(shifter) +  length(upper.shifter) + length(lower.shifter) != length(unique(x$OTU))){
    stop("Sum of abundance groups does not equal input")
  }
  
  list(universal.abundant = univ.abun,
       universal.medium = univ.med,
       universal.rare = univ.rare,
       abundant = abundant,
       moderate = moderate,
       rare = rare,
       specialist = hab.spec,
       cosmopolitan = cosmopolitan,
       shifter = shifter,
       upper.shifter = upper.shifter,
       lower.shifter = lower.shifter
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
unlist(lapply(dna.ab.group, length))
# abundant 1
# moderate 4
# rare 17459
# specialist 1561
# cosmopolitan 46
# shifter 341
# upper.shifter 755
# lower.shifter 17

# create table

tab1 <- data.frame(abundance.groups = c("Universal abundant", "Universal moderate", "Universal rare",
                                        "Abundant", "Moderate", "Rare",
                                        "Specialist",
                                        "Cosmopolitan",
                                        "Shifter",
                                        "Upper shifter",
                                        "Lower shifter"),
           criteria = c("Abundant* in all habitats, never absent",
                        "Moderate* in all habitats, never absent",
                        "Rare* in all habitats, never absent",
                        "Only abundant observations, can be absent",
                        "Only moderate observations, can be absent",
                        "Only rare observations, can be absent",
                        "Only abundant in one habitat type, never abundant in any other habitats",
                        "Shifts between abundant, moderate and rare, but present in all habitats",
                        "Shifts between abundant, moderate and rare",
                        "Shifts in the upper fraction (abundant - moderate) of the rank abundance curve",
                        "Shifts in the lower fraction (moderate - rare) of the rank abundance curve"),
           otus = unlist(lapply(dna.ab.group, length)))
row.names(tab1) <- NULL

colnames(tab1) <- c("Abundance groups", "Criteria", "Categorized OTUs")

# Copy paste to separete LaTeX document
kable(tab1, format = "latex", booktabs = T, caption = "Abundance groups") %>%
  kable_styling(position = "center", latex_options = "scale_down" ) %>%
  row_spec(0) %>%
  footnote(symbol = linebreak(c("Abundant: 47 >= CSS reads, Moderate: 5 >= CSS reads < 47, Rare: 0 > CSS reads < 5")))

# ISME needs excel tables
tab1[nrow(tab1)+1,] <- c("Abundant: 47 >= CSS reads, Moderate: 5 >= CSS reads < 47, Rare: 0 > CSS reads < 5", NA, NA)

# Write the first table into a new workbook
write.xlsx(tab1, file = "./Manuscript/Tables/ISME_tables.xlsx",
           sheetName = "Tab1", append = FALSE, row.names = F, showNA = F)
# Add a second data set in a new worksheet
#write.xlsx(tab2, file = "./Manuscript/Tables/ISME_tables.xlsx", 
#           sheetName="Tab2", append=TRUE)

# Where in the rank abundance curve does community reshuffling occur? ----------------------------
# calculate change in abundance between DNA and RNA between each OTU

cast.abun <- dcast(rel.df, DR.names + OTU ~ DnaType, value.var = "css.reads")
cast.abun <- cast.abun[!is.na(DNA) & !is.na(cDNA),]
cast.abun[, diff.abun := DNA - cDNA]

dist.cast <- dcast(dist.75, DR.names + sample.type.year + Season ~ Metric, value.var = "dist")
#dist.cast <- dcast(dist.3d, DR.names + sample.type.year + Season ~ Metric, value.var = "dist.3d")
#dist.cast <- dcast(dist.1d[(Axis == "Axis.2" & Metric == "Bray") | (Axis == "Axis.3" & Metric == "Sorensen")
#                             ,], DR.names + sample.type.year + Season ~ Metric, value.var = "dist")
dist.cast <- dist.cast[!is.na(dist.cast$Bray),]
dist.cast[, bcs.ratio := Sorensen / Bray]

cast.abun[dist.cast, c("Bray", "Sorensen", "bcs.ratio", "sample.type.year","Season") :=
            list(i.Bray, i.Sorensen, i.bcs.ratio, i.sample.type.year, i.Season), on = .(DR.names)]

#cast.abun[dist.1d[Metric == "Bray" & Axis == "Axis.2",], c("dist.pc2") := list(i.dist), on = .(DR.names)]

cast.abun <- cast.abun[!(diff.abun == 0),] # omit no abundance difference as they do not contribute to differences
cast.abun[cDNA > 0 & DNA == 0,]
# add abundance group names into species dataframe
for(i in 1:length(grp.names)){
  cast.abun[OTU %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
  cast.abun[OTU %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
}

# omit OTUs only in RNA
cast.abun <- cast.abun[!is.na(ab.names),]

#cast.abun[ab.names == "abundant" | ab.names == "moderate", ab.names := "ab.mod"]
cast.abun <- cast.abun[ab.names != "abundant",]

cast.abun$ab.groups <- factor(cast.abun$ab.names, levels = c('moderate', 'rare',
                                                             'specialist',
                                                             'cosmopolitan',
                                                             'shifter',
                                                             'upper.shifter', 'lower.shifter'),
                           labels = c('Moderate', 'Rare',
                                      'Specialist', 'Cosmopolitan', 'Shifter',
                                      'Upper shifter', 'Lower shifter'))

#library(thematic)
#okabe.ito.sel <- okabe_ito(n = 8)
okabe.ito.sel <- c("#000000",
                   "#0072B2", "#009E73", "#CC79A7", "#56B4E9","#F0E442", "#E69F00")


within.p <- ggplot(cast.abun, aes(x = Bray, y = log1p(abs(diff.abun)))) + #dist.pc2
  theme_pubr() +
  geom_point(colour = "grey80", alpha = 0.8, size = 0.5) +
  labs(x = expression(paste(italic("m")["BC"])),
       y = expression(paste("| ", Delta, " OTU CSS reads |"))) +
  annotate(geom = "rect", xmin = 0.19, ymin = 0, 
           xmax = 1, ymax = 8, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  theme(axis.text = element_text(size = 3), axis.title = element_text(size = 5))

(main <- ggplot(cast.abun,  aes(x = Bray, y = log1p(abs(diff.abun)))) + # dist.pc2
    theme_pubr() +
    geom_point(colour = "grey80", alpha = 0.2) +
    #geom_smooth(aes(group = ab.groups, colour = ab.groups, fill = ab.groups),
    #            method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.4) +
    stat_summary_bin(fun=mean, geom="line", aes(colour = ab.groups), bins = 10, size = 1) +
    scale_colour_manual(values = okabe.ito.sel, name = "Spatial\nabundance\ngroups") +
    scale_fill_manual(values = okabe.ito.sel, name = "Spatial\nabundance\ngroups") +
    labs(x = expression(paste(italic("m")["BC"])), #^"PC2"
         y = expression(paste("| ",Delta, " OTU CSS reads |"))) +
    coord_cartesian(ylim = c(0,10)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    #annotate(geom = "text", x = 0.5/2, y = 9, label = "DNA > RNA") +
    annotation_custom(ggplotGrob(within.p), xmin = 0.2, xmax = 0.6, ymin = 6.7, ymax = 10) +
    theme(legend.position = "right"))

ggsave("./Figures/Final/abgroups_mBC_75.png", main,
       width = 14, height = 9 , units = "cm")

# OTU optima in PCoA ----------------------------------------------------------------------------

# fit species onto PCoA ordination
pb.species <-
  data.frame(
    OTU = as.character(colnames(pb.mat)),
    wascores(pb.bray.pcoa$vectors[, 1:3],
             w = pb.mat),
    stringsAsFactors = F
  ) %>% mutate(Metric = "Bray")

#strength.species <- eigengrad(pb.bray.pcoa$vectors[, 1:3],
#          w = pb.mat)

pb.species <- na.omit(pb.species)
setDT(pb.species)
# add abundance group names into species dataframe
for(i in 1:length(grp.names)){
  pb.species[OTU %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
  pb.species[OTU %in% ls.asvs[[grp.names[i]]], ab.names :=  grp.names[i]]
}

sp.melt <- melt(pb.species, id.vars = c("OTU","ab.names"), measure.vars = patterns("^Axis."),
     variable.name = "Axis", value.name = "Scores")

sp.melt <- sp.melt[!is.na(sp.melt$ab.names),]

sp.melt$ab.groups <- factor(sp.melt$ab.names, levels = c('abundant.shifter', 'present.as', 
                                                             'rare.shifter', 'present.rs', 
                                                             'specialist', 'universal.rare'),
                              labels = c('Abundant shifter', 'Universal \nabundant shifter',
                                         'Rare shifter', 'Universal \nrare shifter',
                                         'Specialist', 'Universal rare'))

sp.melt$ab.groups <- factor(sp.melt$ab.names, levels = c('abundant', 'moderate', 'rare',
                                                             'specialist',
                                                             'cosmopolitan',
                                                             'shifter',
                                                             'upper.shifter', 'lower.shifter'),
                              labels = c('Abundant' , 'Moderate', 'Rare',
                                         'Specialist', 'Cosmopolitan', 'Shifter',
                                         'Upper shifter', 'Lower shifter'))

sp.melt$Axis <- factor(sp.melt$Axis, levels = c("Axis.1", "Axis.2", "Axis.3"),
                       labels = c("PC1", "PC2", "PC3"))

labels <- data.frame(Axis = rep(c("PC1", "PC2", "PC3"), times = 2),
                     x = rep(4.5, times = 6), y = c(rep(0.32, times = 3), rep(-0.43, times = 3)),
                     labels = c("Aquatic", "DNA", "Spring",
                                "Terrestrial", "RNA", "Summer/Autumn"))

# combine with PCoA
(vio <- ggplot(sp.melt, aes(x = ab.groups, y = Scores)) +
  theme_pubr() +
  geom_hline(yintercept = 0, colour = "grey20", size = 0.3, linetype = "dotted") +
  facet_grid(.~Axis, scales = "free") +
  geom_boxplot(colour = "grey50", outlier.alpha = 0, width = 0.7) +
  geom_violin(aes(fill = ab.groups), scale = "width", alpha = 0.6) +
  geom_text(data = labels, aes(x = x, y = y, label = labels), colour = "grey30", size = 3.5) +
  scale_fill_manual(values = okabe.ito.sel, name = "Abundance groups") +
  labs(x = "Spatial abundance groups", y = "Species scores") +
  theme(
   axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey30"),
   strip.text = element_text(colour = "white")
  ) +
    guides(fill = "none"))

ggsave("./Figures/Final/abgroups_dist_violin.png", vio,
       width = 15, height = 10, units = "cm")

(prow <- plot_grid(main + theme(legend.position = "bottom",
                                legend.title = element_text(size = 9),
                                legend.text = element_text(size = 8),
                                legend.key.size = unit(0.4, "cm")) +
                     guides(colour = guide_legend(nrow = 4)), vio, labels = "auto", 
          align = "hv", axis = "bt",
          rel_widths = c(1,2)))

ggsave("./Figures/Final/abgroups_rollmean_violin.png", prow,
       width = 25, height = 11, units = "cm")

#-----------------------------------------------------------------------------------------------

## Figure 5. -----------------------------------------------------------------------------------

# Q5: ?

# read in rarefied datasets
min_lib <- c(15000, 25000, 50000)

perm.rar <- select_newest("./Output",
                          "perm.rar_lib_otu99_", by = min_lib)

perm.rar <- c(perm.rar, select_newest("./Output", "201520162017_fin_css_otu99_table_")) # "201520162017_CSS_asvtab"

alpha <- llply(as.list(perm.rar), function(x){
  if(grepl(x, pattern = "fin_css") == T){
    rar <- read.csv(paste0("./Output/", x), sep = ";", dec = ".", stringsAsFactors = F)
    rar <- round(rar, 0)
    data.table::setDT(rar, keep.rownames = "Sample")
    rar <- setDF(rar)
  } else {
    rar <- read.csv(paste0("./Output/", x), sep = ";", stringsAsFactors = F)
    data.table::setDT(rar)
    # make mean column to count
    rar[, iter.mean := round(iter.mean, 0)]
    rar <- setDF(dcast(rar, Sample ~ OTU, value.var = "iter.mean"))
  }
  
  
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

# combine with DR.names and meta data
sumdf <- readRDS("./Objects/summary.meta.with.oldnames.rds")

sumdf <- sumdf[,c("Sample","DR.names", "DnaType","Season"), with = T]
sumdf <- setDT(sumdf %>% distinct()) ; setDT(alpha.df)
alpha.df[sumdf, c("DR.names","Season","DnaType") := list(i.DR.names,
                                                         i.Season,
                                                         i.DnaType), on = .(Sample)]

write.table(alpha.df, paste0("./Output/alpha_div_otu99_summary_", Sys.Date(),".csv"), sep = ",", dec = ".", row.names = F)

#-----------------------------------------------------------------------------------------------
alpha.df <- select_newest("./Output/", "alpha_div_otu99_summary_")
alpha.df <- read.csv(paste0("./Output/", alpha.df), sep = ",", stringsAsFactors = F)

# make to data table
setDT(alpha.df)

# split DNA and RNA
dna.alpha <- alpha.df[alpha.df$DnaType == "DNA",]
rna.alpha <- alpha.df[alpha.df$DnaType == "cDNA",]

temp <- alpha.df[DR.names %in% unique(alpha.df[DnaType == "cDNA",]$DR.names),][Data == "css",]
t <- temp[DR.names %in% DR.names[duplicated(DR.names)],]
cast.alpha <- dcast(temp,
                     DR.names + Season ~ DnaType, value.var = c("Shannon", "Simpson", "Pielou","Chao1"))

ggplot(cast.alpha, aes(x = Pielou_DNA, y = Pielou_cDNA, colour = Season)) +
  geom_point()


###
# DNA
###

# Execute regressions first on DNA diversity
df <- melt(dna.alpha, id.vars = c("DR.names","Data"),
                  measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
                  variable.name = "Index",
                  value.name = "Diversity")

# Only focus on CSS, non rarefied data and two diversity indices
df <- df[Data == "css" & (Index == "Shannon" | Index == "Pielou"),]

# merge with distance
reg.df <- df[dist.75[Metric == "Sorensen",], c("distance", "sample.type.year", "Season") := 
                       list(i.dist, i.sample.type.year, i.Season), on = .(DR.names)]

# remove any NAs, e.g. samples from 2015
reg.df <- reg.df[!is.na(distance),]
setorderv(reg.df, c("Index","Diversity")) # rearrange dataframe

# decide which model is best (e.g. linear or polynomial)
z <- reg.df[Index == "Shannon"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
      distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 -0.03645 p > 0.05
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.47 p > 0.05
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.18 p > 0.05
anova(lm0,lm1) # preferred model is lm1, higher R2 and lowered RSS
# polynomials do not add much, and avoid overfitting and choose a parsimonious model
anova(lm1,lm2) # preferred model lm2
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
    lin <- lm(means$Diversity ~ poly(means$distance, 2))
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
df <- melt(rna.alpha, id.vars = c("DR.names","Data"),
           measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
           variable.name = "Index",
           value.name = "Diversity")

# Only focus on CSS, non rarefied data and two diversity indices
df <- df[Data == "css" & (Index == "Shannon" | Index == "Pielou"),]

# calculate mean of duplicates
df[, .(Diversity = mean(Diversity, na.rm = T)), by = .(DR.names, Index)]

# merge with distance
reg.df <- df[dist.75, c("distance", "sample.type.year") := 
               list(i.dist, i.sample.type.year), on = .(DR.names)]

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
  lin <-  lm(means$Diversity ~ means$distance, 2))

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
