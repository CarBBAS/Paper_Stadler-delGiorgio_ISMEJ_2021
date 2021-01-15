#-- Script for the publication:
#-- Title
#-- Responsible author: Masumi Stadler

# This script is the third of a series of scripts that were used to analyse the data
# used in the publication.

# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("phyloseq", "plyr", "tidyverse", "data.table", # wrangling
              "ggpubr", "ggnewscale", "cowplot", "plotly", "ggpmisc", "gridExtra", # plotting,
              "kableExtra", "xlsx", # making tables for manuscript (LaTeX) and export as excel (for ISME)
              "doMC", # parallel computing
              "vegan", "ape", "ade4", "rstatix") # statistics

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
otu.tab <- select_newest("./Output", "201520162017_fin_css_otu99_table_")
otu.tab <- as.matrix(read.csv(
  paste0("./Output/", otu.tab),
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))
#otu.tab <- otu.tab[, 1:1000]

# row orders need to match between tax.tab and otu.tab
otu.tab <- otu.tab[order(row.names(otu.tab)),]

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
# orders need to match between tax.tab and otu.tab
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
pb <- phyloseq(otu_table(otu.tab, taxa_are_rows = F),
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
  "gray40") #Estuary)

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
  pb.mat <- decostand(pb.mat, "hellinger")
  
  # Try NMDS, no convergence, do not take patterns seriously
  # Just to see what main patterns are in only two axes
  #nmds <- ordinate(dna, "NMDS", "bray")
  #plot_ordination(dna, nmds, type="samples", color="sample.type.year", shape="Season") +
  #  scale_colour_manual(values = colvec)
  
  # PCoA with Bray-Curtis
  pb.bray <- vegdist(pb.mat, method = "bray")
  is.euclid(pb.bray) # FALSE
  pb.bray <- sqrt(pb.bray) # make Euclidean
  is.euclid(pb.bray) # TRUE
# we need euclidean distance for the PERMANOVA later on

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
    coord_cartesian(ylim = c(-0.12, 0.08), xlim = c(-0.42, -0.218), expand = F) + 
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
  coord_cartesian(ylim = c(-0.15, 0.25), xlim = c(-0.12, 0.29), expand = F) + # ensure aspect ratio
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
  annotate(geom = "text", x = c(0.285,0.285), y = c(0.255,0.12), label = c("2015", "2016"), 
           size = 3, colour = "grey50") +
  scale_fill_continuous(type = "viridis", name = "Distance from \nmouth (km)", direction = -1) +
  #scale_fill_viridis_b(name = "Distance from mouth", direction = -1) +
  scale_shape_manual(values = c(21,23), name = "Habitat Type") +
  coord_cartesian(xlim = c(0.03, 0.32), ylim = c(-0.2,0.3), expand = F) + # ensure aspect ratio
  labs(x = paste("PC1"), 
       y = paste("PC2")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
  #       alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
  #       fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))


# annotate collage boxes within main plot
p <- p + 
  annotate(geom = "rect", xmin = -0.42, ymin = - 0.12, 
           xmax = -0.21, ymax = 0.08, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "rect", xmin = 0.03, xmax = 0.32,
           ymin = - 0.2, ymax = 0.3, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "rect", xmin = -0.1, xmax = 0.28,
           ymin = - 0.12, ymax = 0.245, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
  annotate(geom = "text", x = c(-0.405,0.045,-0.085), 
           y = c(0.065,0.285,0.23), label = c("b","c","d"), alpha = 0.8, colour = "grey50")

# combine all plots into one
(collage <- ggarrange(p, 
          ggarrange(zoom.terr, zoom.river, zoom.est, ncol = 3, labels = c("b","c","d"), align = "hv"),
          nrow = 2, labels = c("a"), heights = c(0.6, 0.4))
)

# save
ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_SampleType.tiff"), p,
       width = 12, height = 10, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_SampleType.png"),  p,
       width = 12, height = 10, unit = "cm")

ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_collage.tiff"), collage,
       width = 25, height = 15, unit = "cm")
ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_collage.png"),  collage,
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
# habitat = F[12] = 18.01, R^2 = 0.36321, p = 1e-04
# season = F[2] = 11.014, R^2 = 0.03701, p = 1e-04
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
# habitat alone = F[12] = 36.76, p = < 2.2e-16
# season alone = F[2] = 40.62, p = < 2.2e-16
# combined = F[27] = 18.62, p = < 2.2e-16

# look at dispersion
mod # habitats
mod2 # seasons

# individual comparisons (exploring differences in dispersion among groups)
tuk <-TukeyHSD(mod)
t <- data.frame(Groups = rownames(tuk$group),unlist(tuk$group))
t <- t %>% separate(Groups, into = c("habitat","season","nucacid",
                                     "habitat_y","season_y","nucacid_y"), sep = "_|-")
row.names(t) <- NULL
# all pair-wise comparisons that are significant
t <- t[t$p.adj < 0.05,]

## Figure 2: DNA-RNA -----------------------------------------------------------------------------------
# Q2: Are the DNA and RNA assemblages different?

# remove OTUs that only appear in RNA
pb <- prune_taxa(taxa_names(pb) %in% taxa_names(dna), pb)

# removing of OTUs that only appear in DNA did not have a difference on results

# extract species table with species in columns
pb.mat <- otu_mat(pb)
#pb.mat <- decostand(pb.mat, "hellinger")
# creates strong horse-shoe

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

# extract legends
# add grey box to colour legend
leg1 <- get_legend(all.pcoa$plot + theme(legend.key = element_rect(fill = "gray50")) +
                                           guides(colour = guide_legend(order = 3, size = 2.5,
                                                                 override.aes = list(fill = "grey20",
                                                                                     shape = 21)),
                                                  fill = F, shape = F))
leg2 <- get_legend(all.pcoa$plot + guides(colour = F) +
                     theme(legend.margin = margin(3,5.5,3,5.5)))
# change layout to remove gap between legends
legs <- gtable_rbind(leg2, leg1)
legs$layout[4,c(1,3)] <- c(9,9)
                           
pcoa.plot <- all.pcoa$plot + theme(legend.position = "none")
pcoa.plot2 <- pcoa.23$plot + theme(legend.position = "none")

(p <- ggarrange(pcoa.plot, pcoa.plot2, ncol = 2,
               align = "hv", labels = "auto", legend.grob = legs,
               legend = "right"))

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

tuk <-TukeyHSD(mod)
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

# extract legends
# add grey box to colour legend
# same as bray curtis

# plot Sorensen plots for supplementary
(p <- ggarrange(soren.12$plot + theme(legend.position = "none"),
                soren.32$plot + theme(legend.position = "none"),
                ncol = 2, common.legend = T, legend.grob = legs,
                legend = "right",
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

# Make data frame to calculate residuals to 1:1 line
dist.resid <- dcast(dist.75, DR.names + sample.type.year + Season ~ Metric, value.var = "dist")

# add residuals to one to one line to data frame
dist.resid$resid <- resid(lm(dist.resid$Sorensen-dist.resid$Bray ~ 0))

# Make Bray to Sorensen ratio (only 75%)
#diff.df <- dcast(dist.75, DR.names + sample.type.year + Season ~ Metric, value.var = "dist")
#diff.df <- diff.df[, bcs.ratio := Sorensen / Bray]
#sum.diff <- diff.df[, .(mean =  mean(bcs.ratio, na.rm = T),
#                       conf.int = conf.int(bcs.ratio),
#                       stdev = sd(bcs.ratio, na.rm = T),
#                       n = .N),
#                   by = .(sample.type.year, Season)]

# calculate mean of residuals by habitat and season
sum.diff <- dist.resid[, .(mean =  mean(resid, na.rm = T),
                       conf.int = conf.int(resid),
                       stdev = sd(resid, na.rm = T),
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

dist.resid[, n := .N, by = .(sample.type.year, Season)]
sum.mets <- dist.resid[, .(Bray = mean(Bray),
                          Bray.sd = sd(Bray),
                          Bray.se = sd(Bray) / sqrt(.N),
                          Sorensen = mean(Sorensen),
                          Soren.sd = sd(Sorensen),
                          Soren.se = sd(Sorensen) / sqrt(.N)),
                      by = .(sample.type.year, Season)]

# add new column to split plot into main and side panel
sum.ls <- lapply(list(sum.dist, sum.dist75, sum.diff, sum.mets), function(x){
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

# 1:1 plot with trajectories
plot.df <- sum.ls[[4]]

(one.to.one <- ggplot(plot.df, aes(x = Bray, y = Sorensen)) +
    theme_cust(base_theme = "pubr", border = T) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = Sorensen - Soren.se, ymax = Sorensen + Soren.se,
                      group = paste(sample.type.year, Season, sep ="_")), colour = "gray40", alpha = 0.5) +
    geom_errorbarh(aes(xmin = Bray - Bray.se, xmax = Bray + Bray.se,
                      group = paste(sample.type.year, Season, sep ="_")), colour = "gray40", alpha = 0.5) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
    scale_fill_manual(values = colvec, name = "Habitat Type") +
    scale_shape_manual(values = c(21, 23, 25)) +
  labs(x = expression(paste(italic(m)["BC"])), y = expression(paste(italic(m)["S"]))) +
  annotate(geom = "text", x = 0.33, y = 0.37, label = "1:1", angle = 45) +
    guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
           fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
  )
#guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
#       fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2))))

# Scatter plots
plot.df <- sum.ls[[3]][!is.nan(sum.ls[[3]]$mean),]
# Rename level for plotting
levels(plot.df$sample.type.year)[levels(plot.df$sample.type.year) == "Riverine \nLakes"] <- "Riv. Lakes"
names(colvec)[7] <- "Riv. Lakes"

# add limit points
plot.df$stdev[is.na(plot.df$stdev)] <- 0
plot.df$lim.points <- plot.df$mean + plot.df$stdev

d <- scatter_panels(plot.df, labs = c("Habitat Type", 
                                      "Residuals"))
#expression(paste(italic("m")["S"], " : ", italic("m")["BC"])
           
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
#bray.df <- sum.ls[[2]][!is.nan(sum.ls[[2]]$mean),]
#bray.df <- bray.df[Metric == "Bray",]
#lim.point <- max(bray.df$mean + bray.df$stdev, na.rm = T)
#p <- scatter_panels(bray.df, labs = c("Habitat Type", 
#                                      expression(paste(italic("m")["BC"]))))

# only mS
#soren.df <- sum.ls[[2]][!is.nan(sum.ls[[2]]$mean),]
#soren.df <- soren.df[Metric == "Sorensen",]
#lim.point <- max(soren.df$mean + soren.df$stdev, na.rm = T)
#s <- scatter_panels(soren.df, labs = c("Habitat Type", 
#                                      expression(paste(italic("m")["S"]))))

# make a ratio legend
resid.leg <-
  get_legend(one.to.one + 
  geom_text(aes(label = Bray <= 0.5, size = Bray <= 0.5), alpha = 0) +
  scale_size_manual(values = c(1,3),
                    name =  "Residuals",
                    labels = c("Mass effect","Selection")) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10)) +
    guides(size = guide_legend(override.aes = list(label = c("> 0", "< 0"), size = 3, alpha = 1))))


#left <- ggarrange(p$main +
#                    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lim.point, ymax = Inf),
#                              fill = "gray90", colour = "gray90") +
#            geom_text(aes(y = lim.point + 0.03, colour = Season, label = n), 
# position = position_dodge(width=0.7),
#                      size = 2.5, show.legend = F) +
#            lims(y = c(min(bray.df$mean - bray.df$stdev, na.rm = T), lim.point + 0.03)) +
#            theme(axis.text.x = element_blank(), legend.position = "none"),
#            s$main + 
#              lims(y = c(min(soren.df$mean - soren.df$stdev, na.rm = T), lim.point + 0.03)) +
#              theme(axis.text.x = element_blank(), legend.position = "none"),
#          d$main + geom_point(aes(shape = mean >= 1), alpha = 0) +
#            geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
#            theme(axis.title.y = element_text(margin = margin(r = 10)),
#                  axis.text.x = element_text(size = 7), legend.position = "none"),
#          nrow = 3, align = "v", heights = c(0.34,0.33,0.45), labels = "auto")

#right <- ggarrange(p$side + 
#                     geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lim.point, ymax = Inf),
#                               fill = "gray90", colour = "gray90") +
#                     geom_text(aes(y = lim.point + 0.03, colour = Season, label = n), position = position_dodge(width=0.7),
#                               size = 2.5, show.legend = F) +
#                     lims(y = c(min(bray.df$mean - bray.df$stdev, na.rm = T), lim.point + 0.03)) +
#                     theme(axis.text.x = element_blank()),
#                   s$side +
#                     lims(y = c(min(soren.df$mean - soren.df$stdev, na.rm = T), lim.point + 0.03)) +
#                     theme(axis.text.x = element_blank()),
#                   d$side + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
#                     theme(
#                       axis.text.x = element_text(size = 7)),
#                   nrow = 3, align = "v", heights = c(0.34,0.33,0.45), common.legend = T,
#                   legend.grob = ratio.leg,
#                   legend = "right")

#(scat.plots <- ggarrange(left, right, ncol = 2, widths = c(2,1)))

(conti <- ggarrange(one.to.one,
                   d$main +
                     #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lim.point + 0.01,
                    #               ymax = lim.point + 0.03),
                    #           fill = "gray90", colour = "gray90") +
                     geom_text(aes(y = lim.points + 0.03, group = Season, label = n), 
                               position = position_dodge(width=0.7),
                               size = 2.5, show.legend = F, colour = "gray40") +
                     lims(y = c(min(plot.df$mean - plot.df$stdev, na.rm = T),
                                max(plot.df$mean + plot.df$stdev, na.rm = T) + 0.03)) +
                     theme(axis.title = element_text(size = 14),
                           axis.text = element_text(size = 10),
                           legend.text = element_text(size = 9),
                           legend.title = element_text(size = 10)),
                  d$side + 
                    #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lim.point + 0.01,
                    #              ymax = lim.point + 0.03),
                    #          fill = "gray90", colour = "gray90") +
                    geom_text(aes(y = lim.points + 0.03, group = Season, label = n), 
                              position = position_dodge(width=0.7),
                              size = 2.5, show.legend = F, colour = "gray40") +
                    lims(y = c(min(plot.df$mean - plot.df$stdev, na.rm = T),
                               max(plot.df$mean + plot.df$stdev, na.rm = T) + 0.03)) +
                    theme(axis.title = element_text(size = 14),
                            axis.text = element_text(size = 10),
                          legend.text = element_text(size = 9),
                          legend.title = element_text(size = 10)),
                  ncol = 3, align = "h", widths = c(1,1.5,0.5), labels = c("a","b", NULL),
                  common.legend = T,
                  legend.grob = resid.leg, legend = "right"))

ggsave("./Figures/Final/distBCS_residcontinuum.png", conti,
       width =23, height = 9, units = "cm")


(distdissm.col <- ggarrange(pw.sorbray, scat.plots, ncol = 2, #labels = c("a",""),
                            widths = c(1.5,2)))

ggsave("./Figures/Final/distBC_distresid_scat.png", scat.plots,
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


# How do we interpret the residuals? --------------------------------------------------------
# extract sample pairs above 0 and below 0 of the residuals
over <- as.character(dist.resid[resid > 0,]$DR.names)
below <- as.character(dist.resid[resid < 0,]$DR.names)

#sanity check
length(over) + length(below) == length(unique(dist.resid$DR.names)) # ok

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
# doesn't make a difference

# cast so that we can calculate difference between DNA and RNA of each sample
temp <- dcast(commat, DR.names + OTU ~ DnaType, value.var = c("reads"))
temp[, diff := DNA - RNA]
temp <- temp[!is.na(diff),]

# calculate summary variables
abs.diff <- temp[, .(mean.diff = mean(abs(diff), na.rm = T),
                     max.diff = max(abs(diff), na.rm = T),
                     var.diff = var(abs(diff), na.rm = T)), by = .(DR.names)]

# calculate abundance difference of shared taxa
shar.diff <- temp[DNA > 0 & RNA > 0,]
shar.diff <- shar.diff[, diff := abs(DNA - RNA)][, .(shar.diff = mean(diff)), by = .(DR.names)]

# calculate how many taxa are shared by sample
n.shared <- temp[DNA > 0 & RNA > 0, ]
n.shared <- n.shared[, .(n = .N), by = .(DR.names)]

# calculate richness difference
dna.rich <- temp[DNA > 0,][, .(n = .N), by = .(DR.names)]
rna.rich <- temp[RNA > 0,][, .(n = .N), by = .(DR.names)]
rich.diff <- dna.rich[rna.rich, n.rna := i.n, on = .(DR.names)][, rich.diff := n - n.rna]

# extract distances
#dists <- dcast(dist.75, DR.names ~ Metric, value.var = "dist")

# calculate replacement
# binary transformation
pa.temp <- setDF(temp); setDT(pa.temp)
pa.temp <- pa.temp[, c("DNA",
                    "RNA") := list(ifelse(DNA > 0, 1, 0),
                                   ifelse(RNA > 0, 1, 0))]
replac <- pa.temp[DNA != RNA,][, .(n = .N), by = .(DR.names)]

# calculate abundance replacement (abundance difference of unshared taxa)
ab.replac <- temp[which(pa.temp$DNA != pa.temp$RNA),][, replac.diff := abs(DNA - RNA)][, .(replac.diff = mean(replac.diff)), by = .(DR.names)]

# combine all to one data set
all.met <- abs.diff[n.shared, n.shared := i.n, on = .(DR.names)]
all.met <- all.met[rich.diff, rich.diff := abs(i.rich.diff), on = .(DR.names)]
all.met <- all.met[replac, replac := i.n, on = .(DR.names)]
all.met <- all.met[ab.replac, replac.diff := i.replac.diff, on = .(DR.names)]
all.met <- all.met[shar.diff, on = .(DR.names)]
#all.met <- all.met[dists, c("mBC", "mS") := list(i.Bray, i.Sorensen), on = .(DR.names)]
all.met <- all.met[dist.resid, c("mBC","mS", "resid") := list(i.Bray,
                                                              i.Sorensen,
                                                              i.resid), on = .(DR.names)]

# give >0 and <0 resid categories
all.met[DR.names %in% over, resid.cat := ">0"]
all.met[DR.names %in% below, resid.cat := "<0"]
#all.met <- all.met[!is.na(ratio.cat),]

# add categories
all.met[setDT(sample_df(pb) %>%
                select(DR.names, sample.type.year, Season) %>%
                group_by(DR.names) %>% distinct()), c("sample.type.year","Season") := list(i.sample.type.year,
                                                                                           i.Season),
        on = .(DR.names)]

# explore metrics
all.met[, .(rich.diff = mean(rich.diff, na.rm = T)), by = .(sample.type.year, Season)]
# all positive, meaning that DNA bigger than richness > RNA richness
all.met[, .(shar.diff = mean(shar.diff, na.rm = T)), by = .(sample.type.year, Season)]
# all positive, meaning that DNA > RNA abundance


# quick plots
ggplot(all.met, aes(x = resid.cat, y = n.shared)) +
  geom_boxplot()
ggplot(all.met, aes(x = resid.cat, y = log(mean.diff))) +
  geom_boxplot()
ggplot(all.met, aes(x = resid.cat, y = log(shar.diff))) +
  geom_boxplot()
ggplot(all.met, aes(x = resid.cat, y = rich.diff)) +
  geom_boxplot()
ggplot(all.met, aes(x = resid.cat, y = replac)) +
  geom_boxplot()
ggplot(all.met, aes(x = resid.cat, y = log(replac.diff))) +
  geom_boxplot()

# Kruskal-Wallis -------------------------------------------------------------------------------------
# Statistically test whether there is a difference between these categories
# 1. Mean abundance ----------------------------------------------------------------------------------
# Is there a difference in abundance difference between the categories identified by the mS:mBC ratio?

# Check ANOVA assumptions
# Do we have outliers?
(outliers <- all.met %>%
    group_by(resid.cat) %>%
    rstatix::identify_outliers(mean.diff) %>%
    select(DR.names, resid.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]
#aov.df <- all.met
# Check normality
model <- lm(mean.diff ~ resid.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis instead of one-way ANOVA
(stat.kw <- kruskal.test(mean.diff ~ resid.cat, data = aov.df)) # significantly different (*** level, p < 0.0001)

# Difference in mean abundance difference among >1 and <1 mS:mBC is statistically significant

(m.ab <- ggplot(aov.df, aes(x = resid.cat, y = mean.diff)) +
    theme_pubr() +
    theme(axis.title.x = element_blank()) +
  geom_boxplot(width = 0.3) +
  labs(y = expression(paste(Delta," DNA-RNA abundance")), 
       x = "Residuals") +
  #annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(3.5,3.4,3.4), yend = c(3.5,3.5,3.5)) +
  annotate(geom = "text", x = 1.5, y = 3.6,
           label = paste(abbrev.p(stat.kw$p.value)[2])))
#expression(atop(paste(bar(x)," DNA-RNA abundance"),"difference of all OTUs"))

# 2. Mean abundance of shared taxa ------------------------------------------------------------------
# Is there a difference in abundance of shared OTUs?
# Do we have outliers?
(outliers <- all.met %>%
    group_by(resid.cat) %>%
    rstatix::identify_outliers(shar.diff) %>%
    select(DR.names, resid.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

# Check normality
model <- lm(shar.diff ~ resid.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(shar.diff ~ resid.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(sh.diff <- ggplot(aov.df, aes(x = resid.cat, y = shar.diff)) +
    theme_pubr() +
    theme(axis.title.x = element_blank()) +
    geom_boxplot(width = 0.3) +
    labs(y = expression(atop(paste(Delta," DNA-RNA abundance"),"of shared OTUs")),
         x = "Residuals") +
    #annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(430,420,420), yend = c(430,430,430)) +
    annotate(geom = "text", x = 1.5, y = 450,
             label = paste(abbrev.p(stat.kw$p.value)[2])))

# 3. Mean abundance of unshared taxa ------------------------------------------------------------------
# Is there a difference in abundance of unshared OTUs?
# Do we have outliers?
(outliers <- all.met %>%
   group_by(resid.cat) %>%
   rstatix::identify_outliers(replac.diff) %>%
   select(DR.names, resid.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

# Check normality
model <- lm(replac.diff ~ resid.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(replac.diff ~ resid.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(unsh.diff <- ggplot(aov.df, aes(x = resid.cat, y = replac.diff)) +
    theme_pubr() +
    theme(axis.title.x = element_blank()) +
    geom_boxplot(width = 0.3) +
    labs(y = expression(atop(paste(Delta," DNA-RNA abundance"),"of unshared OTUs")),
         x = "Residuals") +
    #annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(40,39,39), yend = c(40,40,40)) +
    annotate(geom = "text", x = 1.5, y = 42,
             label = paste0("italic(p) =", 
                            abbrev.p(stat.kw$p.value)[1]), parse = TRUE)) # double == for parse



# 4. Number of shared OTUs ------------------------------------------------------------------
# Is there a difference in number of shared taxa between the categories identified by the mS:mBC ratio?
# Do we have outliers?
(outliers <- all.met %>%
    group_by(resid.cat) %>%
    rstatix::identify_outliers(n.shared) %>%
    select(DR.names, resid.cat, is.outlier, is.extreme))

# None are extreme, do not exclude
aov.df <- all.met

# Check normality
model <- lm(n.shared ~ resid.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(n.shared ~ resid.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(n.sh <- ggplot(aov.df, aes(x = resid.cat, y = n.shared)) +
    theme_pubr() +
    theme(axis.title.x = element_blank()) +
    geom_boxplot(width = 0.3) +
    labs(y = "Number of shared OTUs", x = "Residuals") +
    annotate(geom = "text", x = 1.5, y = 800,
             label = paste(abbrev.p(stat.kw$p.value)[2])))


# 5. Number of replaced OTUs ------------------------------------------------------------------
# Is there a difference in number of replaced taxa between the categories identified by the residuals?
# Do we have outliers?
(outliers <- all.met %>%
   group_by(resid.cat) %>%
   rstatix::identify_outliers(replac) %>%
   select(DR.names, resid.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

# Check normality
model <- lm(replac ~ resid.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(replac ~ resid.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(n.r <- ggplot(aov.df, aes(x = resid.cat, y = replac)) +
    theme_pubr() +
    theme(axis.title.x = element_blank()) +
    geom_boxplot(width = 0.3) +
    labs(y = "Number of replaced OTUs", x = "Residuals") +
    #annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(2000,1950,1950), 
    #         yend = c(2000,2000,2000)) +
    annotate(geom = "text", x = 1.5, y = 2100, 
             label = paste(abbrev.p(stat.kw$p.value)[2])))

# 6. Richness difference ------------------------------------------------------------------
# Is there a difference in richness between the categories identified by the residuals?
# Do we have outliers?
(outliers <- all.met %>%
   group_by(resid.cat) %>%
   rstatix::identify_outliers(rich.diff) %>%
   select(DR.names, resid.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

# Check normality
model <- lm(rich.diff ~ resid.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(rich.diff ~ resid.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(n.rich <- ggplot(aov.df, aes(x = resid.cat, y = rich.diff)) +
    theme_pubr() +
    theme(axis.title.x = element_blank()) +
    geom_boxplot(width = 0.3) +
    labs(y = expression(paste(Delta, " DNA-RNA Richness")), x = "Residuals") +
    #annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(1350, 1300, 1300),
    #         yend = c(1350,1350,1350)) +
    annotate(geom = "text", x = 1.5, y = 1400,
             label = paste(abbrev.p(stat.kw$p.value)[2])))


krusk.plots <- ggarrange(m.ab, sh.diff,
                         unsh.diff,
                          n.sh,
                          n.r,
                          n.rich, 
                          nrow = 2, ncol = 3, align = "hv", labels = "auto")
(krusk.plots <- annotate_figure(krusk.plots, bottom = "Residuals"))

ggsave("./Figures/Final/resid_krusk_tests.png", krusk.plots, width = 25, height = 15, units = "cm")

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

cast.abun[dist.resid, c("Bray", "Sorensen", "resid", "sample.type.year","Season") :=
            list(i.Bray, i.Sorensen, i.resid, i.sample.type.year, i.Season), on = .(DR.names)]

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
  theme(axis.text = element_text(size = 5), axis.title = element_text(size = 8))

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
    #annotation_custom(ggplotGrob(within.p), xmin = 0.2, xmax = 0.6, ymin = 6.7, ymax = 10) +
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

# calculate relative abundance 
# colour blind friendly, derived from https://medialab.github.io/iwanthue/
col_vector<-c("#350070","#dcd873","#06a644","#b9ce40","#003499","#6b8aff","#f5d249","#ec83f6","#7beb8b",
              "#c947b1","#01b072","#df3587","#006f26","#ff83d0","#215a00","#99a1ff","#668200",
              "#e77e28","#019cf8","#b5221d","#bd0a35","#b2e294","#840066",
              "#314800","#ffb0ed","#954700","#3d0e52","#ff9b61","#59003b","#ff6e83","#aa7dbf","#620009")

# Calculate relative abundance of each OTU in spatial abundance group
tax.tab <- setDT(as.data.frame(tax.tab), keep.rownames = "OTU")
cast.abun[tax.tab, c("phylum", "class") := list(i.phylum, i.class), on = .(OTU)]

cast.abun[, c("sum.read.dna",
              "sum.read.rna") := list(sum(round(DNA,0)),
                                      sum(round(cDNA,0))), by = .(ab.groups, sample.type.year)]
abund.phyla <- cast.abun[, .(rel.dna = sum(round(DNA,0)) / sum.read.dna,
                             rel.rna = sum(round(cDNA,0)) / sum.read.rna),
                         by = .(phylum, ab.groups, sample.type.year)] %>% 
  distinct()
abund.phyla[, rel.rna := rel.rna * -1] # make negatvive for visualization

#melt
abund.melt <- melt(abund.phyla, 
                   id.vars = c("phylum", "ab.groups", "sample.type.year"),
                   measure.vars = c("rel.dna","rel.rna"),
                   variable.name = "DnaType", value.name = "rel.abund")

ggplot(abund.melt, aes(x = sample.type.year, y = rel.abund, fill = phylum)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_hline(yintercept = 0)+
  facet_wrap(~ab.groups) +
  coord_flip() +
  scale_fill_manual(values = rev(col_vector)) +
  scale_y_continuous(breaks = pretty(abund.melt$rel.abund),
                     labels = abs(pretty(abund.melt$rel.abund)))
# left is RNA, right is DNA


#---------------------#
#------- Done! -------#
#---------------------#
sessionInfo()

#R version 4.0.3 (2020-10-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.5 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

#locale:
#  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C               LC_TIME=en_CA.UTF-8       
#[4] LC_COLLATE=en_CA.UTF-8     LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
#[7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggpmisc_0.3.7     rstatix_0.6.0     ade4_1.7-16       ape_5.4-1        
#[5] vegan_2.5-6       lattice_0.20-41   permute_0.9-5     doMC_1.3.7       
#[9] iterators_1.0.13  foreach_1.5.1     kableExtra_1.3.1  plotly_4.9.2.1   
#[13] cowplot_1.1.0     ggnewscale_0.4.3  ggpubr_0.4.0      data.table_1.13.2
#[17] forcats_0.5.0     stringr_1.4.0     dplyr_1.0.2       purrr_0.3.4      
#[21] readr_1.4.0       tidyr_1.1.2       tibble_3.0.4      ggplot2_3.3.2    
#[25] tidyverse_1.3.0   plyr_1.8.6        phyloseq_1.32.0  

#loaded via a namespace (and not attached):
#  [1] colorspace_2.0-0     ggsignif_0.6.0       ellipsis_0.3.1       rio_0.5.16          
#[5] rprojroot_2.0.2      XVector_0.28.0       fs_1.5.0             rstudioapi_0.13     
#[9] farver_2.0.3         fansi_0.4.1          lubridate_1.7.9.2    xml2_1.3.2          
#[13] codetools_0.2-16     splines_4.0.3        knitr_1.30           pkgload_1.1.0       
#[17] polynom_1.4-0        jsonlite_1.7.1       broom_0.7.2          cluster_2.1.0       
#[21] dbplyr_2.0.0         compiler_4.0.3       httr_1.4.2           backports_1.2.0     
#[25] assertthat_0.2.1     Matrix_1.2-18        lazyeval_0.2.2       limma_3.44.3        
#[29] cli_2.2.0            htmltools_0.5.0      prettyunits_1.1.1    tools_4.0.3         
#[33] igraph_1.2.6         gtable_0.3.0         glue_1.4.2           reshape2_1.4.4      
#[37] Rcpp_1.0.5           carData_3.0-4        Biobase_2.48.0       cellranger_1.1.0    
#[41] vctrs_0.3.5          Biostrings_2.56.0    multtest_2.44.0      nlme_3.1-149        
#[45] xfun_0.19            testthat_3.0.0       openxlsx_4.2.3       rvest_0.3.6         
#[49] lifecycle_0.2.0      gtools_3.8.2         zlibbioc_1.34.0      MASS_7.3-53         
#[53] scales_1.1.1         hms_0.5.3            biomformat_1.16.0    metagenomeSeq_1.30.0
#[57] rhdf5_2.32.4         RColorBrewer_1.1-2   yaml_2.2.1           curl_4.3            
#[61] gridExtra_2.3        stringi_1.5.3        desc_1.2.0           S4Vectors_0.26.1    
#[65] caTools_1.18.0       BiocGenerics_0.34.0  zip_2.1.1            shape_1.4.5         
#[69] bitops_1.0-6         matrixStats_0.57.0   rlang_0.4.9          pkgconfig_2.0.3     
#[73] Wrench_1.6.0         evaluate_0.14        Rhdf5lib_1.10.1      htmlwidgets_1.5.2   
#[77] labeling_0.4.2       tidyselect_1.1.0     magrittr_2.0.1       R6_2.5.0            
#[81] gplots_3.1.0         IRanges_2.22.2       generics_0.1.0       DBI_1.1.0           
#[85] pillar_1.4.7         haven_2.3.1          foreign_0.8-79       withr_2.3.0         
#[89] mgcv_1.8-33          survival_3.2-7       abind_1.4-5          modelr_0.1.8        
#[93] crayon_1.3.4         car_3.0-10           KernSmooth_2.23-17   utf8_1.1.4          
#[97] rmarkdown_2.5        progress_1.2.2       locfit_1.5-9.4       grid_4.0.3          
#[101] readxl_1.3.1         reprex_0.3.0         digest_0.6.27        webshot_0.5.2       
#[105] glmnet_4.0-2         stats4_4.0.3         munsell_0.5.0        viridisLite_0.3.0   
