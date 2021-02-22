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


# Phantom taxa ---------------------------------------------------------------------------------------
otu.tab <- setDT(as.data.frame(otu_table(pb), strings.As.Factors = F), keep.rownames = "Sample")
# melt OTU table into long format
otu.tab <- melt.data.table(
  otu.tab,
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

# Join OTU and meta table
sumdf <- left_join(
  otu.tab,
  met.df[met.df$DadaNames %in% otu.tab$Sample,] %>%
    dplyr::select(
      DadaNames,
      DR.names,
      Year,
      Season,
      DnaType,
      sample.type,
      sample.type.year,
    ),
  by = c("Sample" = "DadaNames")
)

# set back to data.table, order data.table by catchment.area
setDT(sumdf)

# check for duplicates
any(duplicated(sumdf[DnaType == "DNA" & OTU == "OTU_100",]$DR.names) == T)
any(duplicated(sumdf[DnaType == "RNA" & OTU == "OTU_100",]$DR.names) == T)

# CSS reads creates 0.5 reads overwrite all with 1 (integer, makes them 0)
range(sumdf[reads > 0,]$reads, na.rm = T)
sumdf[reads == 0.5, reads := 1]

# cast into wide format
castdf <- dcast(sumdf, DR.names + OTU ~ DnaType, value.var = "reads")

# fill NAs with 0, those are taxa not found in RNA or DNA
castdf[is.na(DNA), DNA := 0]
castdf[is.na(RNA), RNA := 0]

# how many samples have an OTU with RNA > 0 and DNA == 0? (phantom taxa)
castdf[DNA > 0 | RNA > 0, n.all := .N, .(DR.names)] # number of any observations above 0 by sample
temp <- castdf[RNA > 0 & DNA == 0,]
temp <- temp[, .(n = .N, n.all = unique(n.all)), .(DR.names)][, prop := n * 100 / n.all]

setDT(met.df)
View(met.df[DR.names %in% temp[prop != 100 & prop > 50,]$DR.names,] %>% select(DR.names,DnaType,LibrarySize))

# overwrite all DNA observations == 0, where RNA > 0 with 1
castdf[RNA > 0 & DNA == 0, DNA := 1]
castdf[RNA > 0 & DNA == 0,] # check if overwrite was successful

# format back to long format
temp <- melt.data.table(castdf,
                id.vars = c("DR.names","OTU"),
                measure.vars = c("DNA","RNA"),
                variable.name = "DnaType",
                value.name = "reads")

# original and casted data frame have not the same length
# possible that samples dropped when there was not both DNA and RNA
nrow(temp) == nrow(sumdf)

# do row merge/overwrite
sumdf[temp, cor.reads := i.reads, on = .(DR.names, OTU, DnaType)]

# cast into wide format
fin.df <- dcast(sumdf, Sample ~ OTU, value.var = "cor.reads")
# make matrix and give row.names
fin.df <- as.matrix(setDF(fin.df, rownames = fin.df$Sample)[,-1])

# same dimensions as original df?
dim(fin.df)
dim(otu_table(pb)) # yes

write.table(fin.df, paste0("./Output/201520162017_fin_css_otu99_phantomcor_",Sys.Date(),".csv"),
            sep = ";", dec = ".", row.names = T)
rm(temp, castdf, sumdf, otu.tab)

# read in
fin.df <- select_newest("./Output", "201520162017_fin_css_otu99_phantomcor_")
fin.df <- as.matrix(read.csv(
  paste0("./Output/", fin.df),
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))

# Update phyloseq object
pb <- phyloseq(otu_table(fin.df, taxa_are_rows = F),
               sample_data(met.df),
               tax_table(tax.tab))



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
dna.pcoa <- plot_pcoa(pb.bray.pcoa, physeq = dna,  colours = colvec, output = T)

# main PCoA
p <- dna.pcoa$plot + guides(colour = "none")

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

rm(collage, p, zoom.river, zoom.est, zoom.terr)

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
# habitat = F[12] = 18.02, R^2 = 0.363330, p = 1e-04
# season = F[2] = 10.95, R^2 = 0.03679, p = 1e-04
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
# habitat alone = F[12] = 37.97, p = < 2.2e-16
# season alone = F[2] = 41.39, p = < 2.2e-16
# combined = F[27] = 18.99, p = < 2.2e-16

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
# habitat = F[12] = 20.719, R^2 = 0.28898, p = 1e-04
# season = F[2] = 14.923, R^2 = 0.03469, p = 1e-04
# nucleic acid type = F[1] = 25.891, R^2 = 0.03009, p = 1e-04
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
# habitat alone = F[12] = 31.253, p =  < 2.2e-16
# season alone = F[2] = 44.633, p = < 2.2e-16
# nucleic acid type alone = F[1] = 1.594, p = 0.2073
# combined = F[48] = 10.889, p = < 2.2e-16

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
# 199 axes

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

cast.dis <- dcast(dissim.dr, ID + sample.type.year + Season + panels ~ Metric, value.var = "dist")

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
# calculating distance between points in n-dimensional space for both Bray-Curtis and Sorensen

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

# Calculate delta of the two metrics
diff.df <- dcast(dist.75, DR.names + sample.type.year + Season ~ Metric, value.var = "dist")
diff.df <- diff.df[, delta := Bray - Sorensen]

diff.df <- melt(diff.df, id.vars = c("DR.names","sample.type.year","Season"),
                measure.vars = c("delta","Sorensen"),
                variable.name = "Metric", value.name = "dist")
sum.delta <- diff.df[, .(mean =  mean(dist, na.rm = T),
                        conf.int = conf.int(dist),
                        stdev = sd(dist, na.rm = T),
                        n = .N),
                    by = .(sample.type.year, Season, Metric)]


# calculate confidence interval and means of sample type and season combinations
# All axes
sum.dist <- dist.dr[, .(mean =  mean(dist, na.rm = T),
                       conf.int = conf.int(dist),
                       stdev = sd(dist, na.rm = T),
                       n = .N),
                   by = .(Metric, sample.type.year, Season)]

# only 75%
sum.dist75 <- dist.75[, .(mean =  mean(dist, na.rm = T),
                          conf.int = conf.int(dist),
                          stdev = sd(dist, na.rm = T),
                          n = .N),
                      by = .(Metric, sample.type.year, Season)]

# add new column to split plot into main and side panel
sum.ls <- lapply(list(sum.dist, sum.dist75, sum.delta), function(x){
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


##########################
# Lollipop plots of distance
temp <- sum.ls[[3]][!is.nan(sum.ls[[3]]$mean),]
# Rename level for plotting
levels(temp$sample.type.year)[levels(temp$sample.type.year) == "Riverine \nLakes"] <- "Riv. Lakes"
names(colvec)[7] <- "Riv. Lakes"

plot.df <- dcast(temp, sample.type.year + Season + panels ~ Metric, value.var = "mean")

temp$Metric <- factor(temp$Metric, levels = c('delta','Sorensen'), labels = c("Abundance",
                                                                                  "Incidence"))

# Make a fake plot with all the legend units, habitat type, season and metric
legs <- get_legend(ggplot(temp, aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_point(aes(y=mean, colour = sample.type.year),
             size = 3) +
  scale_colour_manual(values = colvec, name = "Habitat type") +
  
  geom_point(aes(y=mean, fill = Metric, shape = Season), size = 3, colour = 'black') +
  scale_fill_manual(values = c("gray40","white"), name = "Metric") +
    scale_shape_manual(values = c(21, 23, 25)) +
  guides(shape = guide_legend(order = 3, override.aes=list(size = 2)),
         colour = guide_legend(order = 2,
                               override.aes=list(size = 2, shape = 21,
                                                 fill = colvec[names(colvec) %in% temp$sample.type.year], 
                                                 colour = "black")),
         fill = guide_legend(order = 1, override.aes=list(shape = 21, size = 2))))


m <- ggplot(plot.df[panels == "main",], aes(x = sample.type.year)) +
  theme_cust(base_theme = "pubr") +
  facet_grid(.~Season, scales = "free_x", space = "free") +
  geom_linerange(aes(ymin = Sorensen, ymax = Sorensen + delta,
                     colour = sample.type.year, group = Season),
                 colour = "gray30", linetype = "solid", position = position_dodge(0.7)) +
  geom_jitter(aes(y = Sorensen, colour = sample.type.year,
                  shape = Season), colour = "black",
              fill = "white", position = position_dodge(0.7), size = 3) +
  scale_fill_manual(values = colvec, name = "Habitat type") +
  scale_shape_manual(values = c(21, 23, 25)) +
  geom_jitter(aes(y = Sorensen + delta, fill = sample.type.year,
                  shape = Season), position = position_dodge(0.7), size = 3) +
  scale_colour_manual(values = colvec) +
  labs(y= "DNA-RNA distance", x = "Habitat Type") +
  lims(y = c(min(plot.df$Sorensen, na.rm = T),
             max(plot.df$Sorensen + plot.df$delta, na.rm = T))) +
  theme(
    axis.text.x = element_blank(),
    axis.text = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 10),
    strip.background = element_rect(fill = "gray20"),
    strip.text = element_text(colour = "white", size = 10)) +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))


s <- ggplot(plot.df[panels == "side",], aes(x = sample.type.year)) +
  theme_cust(base_theme = "pubr") +
  facet_grid(.~"Meta-community") +
  geom_linerange(aes(ymin = Sorensen, ymax = Sorensen + delta,
                     colour = sample.type.year, group = Season),
                 colour = "gray30", linetype = "solid", position = position_dodge(0.7)) +
  geom_jitter(aes(y = Sorensen, colour = sample.type.year,
                  shape = Season), colour = "gray20",
              fill = "white", position = position_dodge(0.7), size = 3) +
  scale_fill_manual(values = colvec) +
  scale_shape_manual(values = c(21, 23, 25)) +
  geom_jitter(aes(y = Sorensen + delta, fill = sample.type.year,
                  shape = Season), position = position_dodge(0.7), size = 3) +
  scale_colour_manual(values = colvec) +
  labs(y= "DNA-RNA distance", x = "Habitat Type") +
  lims(y = c(min(plot.df$Sorensen, na.rm = T),
             max(plot.df$Sorensen + plot.df$delta, na.rm = T))) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = "gray20"),
    strip.text = element_text(colour = "white", size = 10)
  ) + guides(shape = "none", fill = "none")

# make a numeric sample type column to add smooth line
plot.df[sample.type.year == "Soil", num.hab := 1]
plot.df[sample.type.year == "Soilwater", num.hab := 2]
plot.df[sample.type.year == "Stream", num.hab := 3]
plot.df[sample.type.year == "Upriver", num.hab := 4]
plot.df[sample.type.year == "Reservoirs", num.hab := 5]
plot.df[sample.type.year == "Downriver", num.hab := 6]
plot.df[sample.type.year == "Estuary", num.hab := 7]

m.down <- ggplot(plot.df[panels == "main",], ) +
  theme_cust(base_theme = "pubr") +
  facet_grid(.~Season, scales = "free_x", space = "free") +
  geom_line(aes(x = num.hab,
                  y = delta),
              stat = "smooth", method = "loess", span = 0.7, se = F,
            alpha = 0.5, colour= "gray10", size = 0.8) +
  geom_jitter(aes(x = num.hab, y = delta, fill = sample.type.year,
                  shape = Season), 
              position = position_dodge(0.7), size = 3) +
  scale_fill_manual(values = colvec, name = "Habitat type") +
  scale_linetype_manual(values = c("solid","dashed","dotted")) +
  scale_shape_manual(values = c(21, 23, 25)) +
  labs(y= expression(paste(Delta, " Distances")), x = "Habitat Type") +
  scale_x_continuous(breaks = 1:length(unique(as.character(plot.df[panels == "main",]$sample.type.year))),
                     labels = unique(as.character(plot.df[panels == "main",]$sample.type.year))) +
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        axis.text = element_text(size = 8),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2))) +
  scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(min(plot.df$delta, na.rm = T)-0.03,
                                                       max(plot.df$delta, na.rm = T)+0.01))

s.down <-ggplot(plot.df[panels == "side",], ) +
  theme_cust(base_theme = "pubr") +
  facet_grid(.~"Meta-community") +
  geom_jitter(aes(x = sample.type.year, y = delta, fill = sample.type.year,
                  shape = Season), 
              position = position_dodge(0.7), size = 3) +
  scale_fill_manual(values = colvec, name = "Habitat type") +
  scale_shape_manual(values = c(21, 23, 25)) +
  labs(y= expression(paste(Delta, " Distances")), x = "Habitat Type") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(min(plot.df$delta, na.rm = T)-0.03,
                                                         max(plot.df$delta, na.rm = T)+0.01)) +
  guides(fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

down <- ggarrange(m.down,s.down, ncol = 2, widths = c(3,1), legend = "none")

up <- ggarrange(m,s, ncol = 2, nrow = 1,  widths = c(3,1), legend = "none")

(dist.p <- ggarrange(up, down, nrow = 2, ncol = 1, common.legend = T, legend.grob = legs, legend = "right",
          heights = c(2,1), align = "hv"))

ggsave("./Figures/Final/distance_lollipop.png", dist.p, 
       width = 18, height = 12, units = "cm", dpi = 300)
ggsave("./Figures/Final/distance_lollipop.tiff", dist.p, 
       width = 18, height = 12, units = "cm", dpi = 300)


# What is behind these patterns? --------------------------------------------------------

# Examine what the proportion of active vs inactive OTUs is in each habitat
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
       c("DR.names", "DnaType","sample.type.year","Season") := list(i.DR.names, i.DnaType,
                                                                    i.sample.type.year, i.Season),
       on = .(Sample)]

# keep only the ones that have both DNA and RNA
#commat <- commat[duplicated(DR.names),]
# doesn't make a difference

# cast so that we can calculate difference between DNA and RNA of each OTU
temp <- dcast(commat, DR.names + sample.type.year + OTU ~ DnaType, value.var = c("reads"))
temp[, diff := DNA - RNA]
temp <- temp[!is.na(diff),]

# calculate the number of positive observations per sample
all <- temp[DNA != 0| RNA != 0, .(n.obs.all = .N), .(DR.names)]
# calculate the number of DNA observations per sample
dna.n <- temp[DNA > 0, .(n.dna = .N), .(DR.names)]
# calculate the number of RNA observations per sample
rna.n <- temp[RNA > 0, .(n.rna = .N), .(DR.names)]
all[dna.n, n.dna := i.n.dna, on = .(DR.names)]
all[rna.n, n.rna := i.n.rna, on = .(DR.names)]

# remove samples that don't have both DNA and RNA
#all <- all[!is.na(n.dna) & !is.na(n.rna),]

# calculate percentage
all[, active := n.rna * 100 / n.obs.all]
all[, inactive := 100 - active]

# melt
all <- melt.data.table(all, id.vars = c("DR.names"),
                measure.vars = c("active", "inactive"),
                variable.name = "Fraction", value.name = "Perc")

# add Seasons
all[met.df, c("sample.type.year","Season") := list(i.sample.type.year, i.Season), on =.(DR.names)]

# summarise
sum.perc <- all[, .(mean = mean(Perc, na.rm = T),
        sd = sd(Perc, na.rm = T)), by = .(sample.type.year, Season, Fraction)]
sum.perc[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater",
                                                                   "Stream","Upriver","Reservoirs",
                                                                   "Downriver", "Estuary",
                                                                   "Sediment", "Tributary",
                                                                   "Riverine \nLakes", "Lake"))]
sum.perc[, Fraction := factor(Fraction, levels = c("active","inactive"))]
# Rename level for plotting
levels(sum.perc$sample.type.year)[levels(sum.perc$sample.type.year) == "Riverine \nLakes"] <- "Riv. Lakes"

# add panel ID
sum.perc[, panels := "main"]
sum.perc[sample.type.year == "Tributary" |
    sample.type.year == "Lake" |
    sample.type.year == "Riv. Lakes" |
    sample.type.year == "Sediment", panels := "side"]
#sum.perc[, sample.type.year := factor(sample.type.year, levels = c("Soil","Sediment",
#                                                            "Soilwater",
#                                                            "Stream", "Tributary",
#                                                            "Riv. Lakes", "Headwater \nPonds", "Lake",
#                                                            "Upriver",
#                                                            "Reservoirs","Downriver",
#                                                            "Estuary"))]

(prop.samp <- ggplot(sum.perc[panels == "main",], aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = mean, fill = Fraction), colour = "gray20"
           ) +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  scale_fill_viridis_d(option = "cividis", direction = -1) +
  labs(x = "Habitat Type", y = "% of bacterial OTUs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 10),
        axis.title.x = element_blank()))

(prop.side <- ggplot(sum.perc[panels == "side",], aes(x = sample.type.year)) +
    theme_cust("pubr") +
    facet_grid(.~Season, scale = "free_x", space = "free") +
    geom_col(aes(y = mean, fill = Fraction), colour = "gray20",
    ) +
    scale_fill_viridis_d(option = "cividis", direction = -1) +
    labs(x = "Habitat Type", y = "% of bacterial OTUs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          strip.background = element_rect(fill = "gray20"),
          strip.text = element_text(colour = "white", size = 10),
          axis.title = element_blank(),
          axis.text.y = element_blank()))

(act.prop <- ggarrange(prop.samp, prop.side, ncol = 2, common.legend = T, widths = c(2,1), labels = "auto",
          vjust =0))

ggsave("./Figures/Final/inactive_active_prop.png", act.prop, width = 8, height =3)

# Calculate percentage of activity of taxa first found in...
# Soil > Soilwater > Stream > Upriver > Reservoir > Downriver > Estuary
# Focus only on the direct continuum
# where were the OTUs first observed (based on DNA)
temp <- temp[(sample.type.year == "Soil" |
        sample.type.year == "Soilwater" |
        sample.type.year == "Stream" |
        sample.type.year == "Upriver" |
        sample.type.year == "Reservoirs" |
        sample.type.year == "Downriver" |
        sample.type.year == "Estuary"),]

# all OTUs in...
soil <- unique(as.character(temp[DNA > 0 & sample.type.year == "Soil",]$OTU))
soilwater <- unique(as.character(temp[DNA > 0 & sample.type.year == "Soilwater",]$OTU))
stream <- unique(as.character(temp[DNA > 0 & sample.type.year == "Stream",]$OTU))
upriver <- unique(as.character(temp[DNA > 0 & sample.type.year == "Upriver",]$OTU))
reservoir <- unique(as.character(temp[DNA > 0 & sample.type.year == "Reservoirs",]$OTU))
downriver <- unique(as.character(temp[DNA > 0 & sample.type.year == "Downriver",]$OTU))
estuary <- unique(as.character(temp[DNA > 0 & sample.type.year == "Estuary",]$OTU))

# now, remove OTUs from vectors that were found before
soilwater <- soilwater[!soilwater %in% soil]
stream <- stream[!stream %in% unique(c(soil,soilwater))]
upriver <- upriver[!upriver %in% unique(c(soil,soilwater,stream))]
reservoir <- reservoir[!reservoir %in% unique(c(soil,soilwater,stream,upriver))]
downriver <- downriver[!downriver %in% unique(c(soil,soilwater,stream, upriver, reservoir))]
estuary <- estuary[!estuary %in% unique(c(soil,soilwater,stream, upriver, reservoir, downriver))]

# sanity check
all.otus <- unique(as.character(temp[DNA > 0 & (sample.type.year == "Soil" |
                                                   sample.type.year == "Soilwater" |
                                 sample.type.year == "Stream" |
                                 sample.type.year == "Upriver" |
                                sample.type.year == "Reservoirs" |
                                sample.type.year == "Downriver" |
                                sample.type.year == "Estuary"),]$OTU))
length(c(soil, soilwater, stream, upriver, reservoir, downriver, estuary)) == length(all.otus)

# make OTU character
temp[, OTU := as.character(OTU)]
# create new column with first observed category for each OTU
temp[OTU %in% soil, first.obs := "Soil"]
temp[OTU %in% soilwater, first.obs := "Soilwater"]
temp[OTU %in% stream, first.obs := "Stream"]
temp[OTU %in% upriver, first.obs := "Upriver"]
temp[OTU %in% reservoir, first.obs := "Reservoirs"]
temp[OTU %in% downriver, first.obs := "Downriver"]
temp[OTU %in% estuary, first.obs := "Estuary"]

# make factor
temp[, first.obs := factor(first.obs, levels = c("Soil", "Soilwater", "Stream",
                                                 "Upriver", "Reservoirs", "Downriver",
                                                 "Estuary"))]
# remove all OTUs that are found outside the direct continuum
first.df <- temp[!is.na(first.obs),]
# all OTUs are characterized
any(!(first.df$OTU %in% c(soil,soilwater,stream, upriver,reservoir,downriver,estuary)) == T)

# sanity check
first.df[is.na(DNA) & RNA >0,]
first.df[DNA == 0 & RNA >0,]

# merge with some meta data
first.df[met.df, c("sample.type.year","Season") := list(i.sample.type.year,
                                                       i.Season), on = .(DR.names)]

# remove any NAs in the dataset
first.df[is.na(DNA),]
first.df[is.na(RNA),]

# calculate number of observations per habitat type and season = Total
inac.all <- first.df[DNA > 0 & RNA == 0,
         .(otu.n.all.dna = .N,
         sum.reads.all.dna = sum(DNA, na.rm = T)), by = .(sample.type.year, Season)] #DR.names
ac.all <- first.df[RNA > 0,
                    .(otu.n.all.rna = .N,
                      sum.reads.all.rna = sum(RNA, na.rm = T)), by = .(sample.type.year, Season)]  # DR.names

# calculate how many OTUs and reads are allocated to the first observed categories
inac <- first.df[DNA > 0 & RNA == 0, .(otu.n = .N,
                               sum.reads = sum(DNA, na.rm = T)), .(first.obs, sample.type.year, Season)] # DR.names
ac <- first.df[RNA > 0, .(otu.n = .N,
                               sum.reads = sum(RNA, na.rm = T)), .(first.obs, sample.type.year, Season)]  # DR.names

# merge together
inac[inac.all, c("otu.n.all.dna","sum.reads.all.dna") :=
          list(i.otu.n.all.dna, i.sum.reads.all.dna), on = .(sample.type.year, Season)] #DR.names
ac[ac.all, c("otu.n.all.rna","sum.reads.all.rna") :=
          list(i.otu.n.all.rna, i.sum.reads.all.rna), on = .(sample.type.year, Season)]  #DR.names

inac[ac, c("otu.n.all.rna", "sum.reads.all.rna","otu.n.rna","sum.reads.rna") :=
          list(i.otu.n.all.rna, i.sum.reads.all.rna, i.otu.n, i.sum.reads),
     on = .(first.obs, sample.type.year, Season)] #DR.names

# calculate percentages
inac[,c("perc.otu", "perc.reads",
           "perc.otu.rna", "perc.reads.rna") := list(otu.n * 100 / otu.n.all.dna,
                                                 sum.reads * 100 / sum.reads.all.dna,
                                                 otu.n.rna * 100 / otu.n.all.rna,
                                                 sum.reads.rna * 100 / sum.reads.all.rna)]
# Sanity check
inac[,.(sum.otu = sum(perc.otu),
        sum.reads = sum(perc.reads),
        sum.otu.rna = sum(perc.otu.rna),
        sum.reads.rna = sum(perc.reads.rna)),, by =.(sample.type.year, Season)]

# re-assign factor levels
plot.df <- inac
plot.df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                  "Upriver", "Reservoirs", "Downriver",
                                                                  "Estuary"))]

# Plot
d <- ggplot(plot.df, aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = perc.otu, fill = first.obs), colour = "gray20") +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 9)) +
  labs(x = "Habitat type", y = "% of OTUs \nin inactive fraction") +
  scale_fill_manual(values = colvec[names(colvec) %in% plot.df$sample.type.year],
                    name = "First \ndetected in")

r <- ggplot(plot.df, aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = perc.otu.rna, fill = first.obs), colour = "gray20") +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 9),
        #strip.background = element_blank(),
        #strip.text = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "Habitat type", y = "% of OTUs \nin reactive fraction") +
  scale_fill_manual(values = colvec[names(colvec) %in% plot.df$sample.type.year])

dr <- ggplot(plot.df, aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = perc.reads, fill = first.obs), colour = "gray20") +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 9)) +
  labs(x = "Habitat type", y = "% of reads in total DNA reads\nof inactive fraction") +
  scale_fill_manual(values = colvec[names(colvec) %in% plot.df$sample.type.year],
                    name = "First observed in")

rr <- ggplot(plot.df, aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = perc.reads.rna, fill = first.obs), colour = "gray20") +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  labs(x = "Habitat type", y = "% of reads in total RNA reads\nof active fraction") +
  scale_fill_manual(values = colvec[names(colvec) %in% plot.df$sample.type.year],
                    name = "First observed in")

first.obs.p <- ggarrange(d,dr,r,rr,
                          ncol = 2, nrow =2, common.legend = T, heights = c(1,1.2))
(first.obs.p <- annotate_figure(first.obs.p, bottom = text_grob("Habitat Type")))


suppl.otu <- ggarrange(d,r, ncol = 2, common.legend = T)
(suppl.otu <- annotate_figure(suppl.otu, bottom = text_grob("Habitat Type")))

frac.first.obs <- ggarrange(dr, rr, ncol = 1, nrow = 2, common.legend = T, heights = c(1,1.2), legend = "none",
                            labels = "auto")
(frac.first.obs <- annotate_figure(frac.first.obs, bottom = text_grob("Habitat Type")))

ggsave("./Figures/Final/inactive_active_firstobs_all.png",
       first.obs.p, height = 12, width = 18, units = "cm")

#ggarrange(prop.samp, ggarrange(d,r,ncol = 1, nrow = 2, common.legend = T, heights = c(1,1.2),
#                               legend = "right"),
#          ncol = 1, nrow = 2, heights = c(1,2))




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
# we do this exercise on the original reads, not the one corrected for phantom taxa
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


# What abundance groups are those in first observed groups, what is their proportion in each habitat?
# Add abundance group to each observation
first.df[DNA >= 47, abg.dna := "Abundant"]
first.df[DNA < 47 & DNA >= 5, abg.dna := "Moderate"]
first.df[DNA < 5, abg.dna := "Rare"]
#first.df[DNA >= 66, abg.rna := "Abundant"]
#first.df[DNA < 66 & DNA >= 9, abg.rna := "Moderate"]
#first.df[DNA < 9, abg.rna := "Rare"]

# Add spatial abundance group to each OTU
for(i in 1:length(grp.names)){
  first.df[OTU %in% ls.asvs[[grp.names[i]]], spab.names :=  grp.names[i]]
}

# calculate percentage of each first observed group's abundance category in habitat type
# n = number of OTUs for each abg category within habitat type/Season
# n.all = number of OTUs in all categories
# sum.reads = sum of reads belonging to OTUs in each category
# sum.reads.all = sum of reads of all OTUs

# there are some taxa without a spatial abundance group
View(first.df[is.na(spab.names) & (DNA > 0 | RNA > 0),])
# these are all phantom taxa, assign them to "rare"
first.df[is.na(spab.names), spab.names := "rare"]

#-----------------------------------

# calculate how many OTUs and reads within each allocated first observed category are
# abundant, moderate or rare
inac.ab <- first.df[DNA > 0 & RNA == 0,
         .(otu.n.ab = .N,
           sum.reads.ab = sum(DNA, na.rm = T)), 
         by = .(first.obs, abg.dna, sample.type.year, Season)] #DR.names
ac.ab <- first.df[RNA > 0,
                  .(otu.n.ab = .N,
                    sum.reads.ab = sum(RNA, na.rm = T)), 
                  by = .(first.obs, abg.dna, sample.type.year, Season)] #DR.names

inac.ab[ac.ab, c("otu.n.ab.rna", "sum.reads.ab.rna") :=
          list(i.otu.n.ab, i.sum.reads.ab), on = .(first.obs, abg.dna, sample.type.year, Season)]

# add original percentage data
plot.ab <- inac.ab[plot.df, on = .(sample.type.year, Season, first.obs)]

# calculate percentage of each abudance group within the first observed fraction
plot.ab[, perc.otu.dna.ab := otu.n.ab * perc.otu / otu.n]
plot.ab[, perc.otu.rna.ab := otu.n.ab.rna * perc.otu.rna / otu.n.rna]

plot.ab[, perc.reads.dna.ab := sum.reads.ab * perc.reads / sum.reads]
plot.ab[, perc.reads.rna.ab := sum.reads.ab.rna * perc.reads.rna / sum.reads.rna]

# make a numeric sample type column to add smooth line
plot.ab[sample.type.year == "Soil", num.hab := 1]
plot.ab[sample.type.year == "Soilwater", num.hab := 2]
plot.ab[sample.type.year == "Stream", num.hab := 3]
plot.ab[sample.type.year == "Upriver", num.hab := 4]
plot.ab[sample.type.year == "Reservoirs", num.hab := 5]
plot.ab[sample.type.year == "Downriver", num.hab := 6]
plot.ab[sample.type.year == "Estuary", num.hab := 7]

ratio.sum <- first.df[ratio > 0 & !is.infinite(ratio), .(mean.ratio = mean(ratio, na.rm = T),
                                                    sd = sd(ratio, na.rm = T)),
                 by = .(first.obs, abg.dna, sample.type.year, Season)]

ratio.sum <- ratio.sum[plot.ab, on = .(first.obs, abg.dna, sample.type.year, Season)]

ggplot(ratio.sum) +
  theme_cust("pubr") +
  facet_grid(sample.type.year~Season) +
  #geom_line(aes(x = log(mean.ratio), 
  #               y = log(perc.reads.rna.ab), colour = first.obs)) +
  geom_point(aes(x = log(mean.ratio), 
                y = log(perc.reads.rna.ab / otu.n.ab.rna), fill = first.obs, shape = abg.dna), size = 1.5) +
  geom_errorbarh(aes(xmax = log(mean.ratio) + log(sd), xmin = log(mean.ratio) - log(sd),
                     y = log(perc.reads.rna.ab / otu.n.ab.rna),
                     colour = first.obs, group = abg.dna)) +
  scale_shape_manual(values = c(21,23,25), name = "Local DNA\nabundance") +
  scale_fill_manual(values = colvec, name = "First \ndetected in") +
  scale_colour_manual(values = colvec) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray20") +
  labs(x = "Mean DNA:RNA ratio\n(log-scale)", y = "% of CSS reads in reactive fraction \n(normalized by number of OTUs, log-scale)") +
  theme(panel.background = element_rect(fill = "gray90"),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 9),
        legend.position = "top") +
  guides(shape = guide_legend(order = 1, override.aes=list(size = 3)),
         fill = guide_legend(order = 2, override.aes=list(shape=21, size = 3)),
         colour = "none")

# Ratio
# NaN = 0 / 0
# Inf DNA > 0 & RNA == 0

# active portion
mean.df <- first.df[RNA > 0, 
                    .(sum.rna.byotu = sum(RNA, na.rm = T),
                      mean.ratio = mean(ratio, na.rm = T)),
         by = .(OTU, sample.type.year, abg.dna, first.obs, Season)]
mean.df[ac.all, "sum.reads.all.rna" := i.sum.reads.all.rna, on = .(sample.type.year,Season)]

mean.df[, perc.sum.rna.byotu := sum.rna.byotu * 100 / sum.reads.all.rna]

# sanity check
mean.df[, .(sum.all = sum(perc.sum.rna.byotu)),
        by = .(sample.type.year, Season)]

# make quartile bins
mean.df[perc.sum.rna.byotu < summary(perc.sum.rna.byotu)[3], quan.bin := 1] #  < 50%
mean.df[perc.sum.rna.byotu >= summary(perc.sum.rna.byotu)[3], quan.bin := 2] # >= 50%

any(is.na(mean.df$quan.bin))
# summarise by bins
mean.df <- mean.df[,.(perc.sum.rna.otu = mean(perc.sum.rna.byotu, na.rm = T),
           mean.ratio = mean(mean.ratio, na.rm = T)),
        by = .(sample.type.year, abg.dna, first.obs, Season, quan.bin)]

mean.df[sample.type.year == "Soil", num.hab := 1]
mean.df[sample.type.year == "Soilwater", num.hab := 2]
mean.df[sample.type.year == "Stream", num.hab := 3]
mean.df[sample.type.year == "Upriver", num.hab := 4]
mean.df[sample.type.year == "Reservoirs", num.hab := 5]
mean.df[sample.type.year == "Downriver", num.hab := 6]
mean.df[sample.type.year == "Estuary", num.hab := 7]

# overwrite soil colour
new.col <- colvec
new.col[names(new.col) == "Soil"] <- "gold"
ann_df <- data.frame(mean.ratio = -1.5,
                     perc.sum.rna.otu = as.numeric(summary(log(mean.df$perc.sum.rna.otu))[3]),
                     Season = factor("Spring",levels = c("Spring","Summer","Autumn")),
                     sample.type.year = factor("Soil"))
line_df <- data.frame(Season = c("Spring","Summer","Spring", "Summer", "Spring", "Summer",
                                 "Spring", "Summer", "Autumn","Spring", "Summer", "Autumn",
                                 "Spring", "Summer", "Autumn",
                                 "Summer"),
                      sample.type.year = c("Soil","Soil","Soilwater","Soilwater", "Stream","Stream",
                                           "Upriver","Upriver","Upriver",
                                           "Downriver","Downriver","Downriver",
                                           "Reservoirs","Reservoirs","Reservoirs",
                                           "Estuary"))
line_df$x <- 0 ; line_df$y <- summary(log(mean.df$perc.sum.rna.otu))[3]
line_df$Season <- factor(line_df$Season, levels = c("Spring","Summer","Autumn"))
setDT(line_df)
line_df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                  "Upriver", "Reservoirs", "Downriver",
                                                                  "Estuary"))]

mean.df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                  "Upriver", "Reservoirs", "Downriver",
                                                                  "Estuary"))]


(perc.cont <- ggplot(mean.df) +
  theme_cust("pubr")+
  facet_grid(sample.type.year~Season, scale = "free_x", space = "free") +
  geom_vline(data = line_df,
             aes(xintercept = x), linetype = "solid", colour = "gray90") +
  geom_hline(data = line_df,
             aes(yintercept = y),
             linetype = "solid", colour = "gray40") +
  geom_point(aes(x = log(mean.ratio), y = log(perc.sum.rna.otu), colour = first.obs,
                 shape = abg.dna, fill = as.character(quan.bin)), alpha = 0.9, size =2) +
  geom_text(data = ann_df,
            aes(y = perc.sum.rna.otu +1.5, x = mean.ratio), label = "Median", size =3) +
  scale_colour_manual(values = new.col,
                      name = "First \ndetected in") +
  scale_shape_manual(values = c(21,23,25),
                     name = "Local DNA\nabundance") +
  scale_fill_manual(values = c("gray80","white"), labels = c("< Median",
                                                    "> Median"),
                    name = "Bins of \n% RNA reads") +
  scale_y_continuous(limits = c(min(log(mean.df$perc.sum.rna.otu)), 
                                max(log(mean.df$perc.sum.rna.otu)) + 1)) +
  theme(panel.grid = element_blank(),
        legend.position = "left",
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 9)) +
  labs(x = "Mean DNA:RNA ratio\n(log-scale)",
       y = "Mean % OTU contribution to\ntotal RNA reads (log-scale)") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 3)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 3)),
         colour = guide_legend(order = 3, override.aes=list(size = 3))))

(perc.cont <- annotate_figure(perc.cont, right = text_grob("Habitat Type", rot = 270, hjust = 0.7)))

(first.obs.col <- ggarrange(frac.first.obs, perc.cont, nrow = 1, ncol = 2, labels = c("","c"), hjust = -10, widths = c(1,2)))

ggsave("./Figures/Final/first.observed.contribution.png", first.obs.col, width = 25, height = 15, units = "cm", dpi = 300)
#----------------------------------------------------------------------

# What's the proportion of first observed OTUs' spatial abundance group?
first.spag <- first.df[,.(spab.names = unique(spab.names)),
                       by = .(OTU, first.obs)][, n := .N, by = .(first.obs)][, n.spag := .N, by = .(first.obs, spab.names)]

#calculate percentage
spag.perc <- first.spag[, .(n = unique(n),
                            n.spag = unique(n.spag)),
                        by = .(first.obs, spab.names)][, perc.spag := n.spag * 100 / n, ]

# sanity check
spag.perc[, .(sum = sum(perc.spag)), by = .(first.obs)]

spag.perc$spab.names <- factor(spag.perc$spab.names, levels = c('abundant','moderate', 'rare',
                                                             'specialist',
                                                             'cosmopolitan',
                                                             'shifter',
                                                             'upper.shifter', 'lower.shifter'),
                              labels = c('Abundant','Moderate', 'Rare',
                                         'Specialist', 'Cosmopolitan', 'Shifter',
                                         'Upper shifter', 'Lower shifter'))

#library(thematic)
#okabe.ito.sel <- okabe_ito(n = 8)
okabe.ito.sel <- c("#000000", "gray70",
                   "#0072B2", "#009E73", "#CC79A7", "#56B4E9","#F0E442", "#E69F00")

# plot this
frac.spag <- ggplot(spag.perc) +
  theme_cust("bw") +
  geom_col(aes(x = first.obs, y = perc.spag, fill = spab.names)) +
  scale_fill_manual(values = okabe.ito.sel, name = "Spatial\nabundance group") +
  labs(x = "First detected in", y = "% of\nbacterial OTUs") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
  

# so, how do they behave along the continuum?
# active portion
mean.df <- first.df[RNA > 0, 
                    .(sum.rna.byotu = sum(RNA, na.rm = T),
                      mean.ratio = mean(ratio, na.rm = T),
                      spab.names = unique(spab.names)),
                    by = .(OTU, sample.type.year, abg.dna, first.obs, Season)]
mean.df[ac.all, "sum.reads.all.rna" := i.sum.reads.all.rna, on = .(sample.type.year,Season)]

mean.df[, perc.sum.rna.byotu := sum.rna.byotu * 100 / sum.reads.all.rna]

mean.df$spab.names <- factor(mean.df$spab.names, levels = c('abundant','moderate', 'rare',
                                                                'specialist',
                                                                'cosmopolitan',
                                                                'shifter',
                                                                'upper.shifter', 'lower.shifter'),
                               labels = c('Abundant','Moderate', 'Rare',
                                          'Specialist', 'Cosmopolitan', 'Shifter',
                                          'Upper shifter', 'Lower shifter'))

spag.box <- ggplot(mean.df) +
  theme_cust("bw") +
  geom_boxplot(aes(x = first.obs, y = log(perc.sum.rna.byotu), fill = spab.names),
               outlier.size = .5,  position = position_dodge(0.7), varwidth = F) +
  scale_fill_manual(values = okabe.ito.sel) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "First detected in", y = "% OTU contribution to total RNA\n(log-scale)")

(spag.col <- ggarrange(frac.spag +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank()), spag.box, nrow = 2, common.legend = T, align = "v", legend = "right",
          heights = c(0.7, 2), labels = "auto"))

ggsave("./Figures/Final/spag_contribution.png", spag.col,
       height = 13, width = 14, units = "cm")

#-------------------------------------------------------------------------------------------------

ggplot(mean.df) +
  theme_cust("pubr")+
  facet_grid(.~Season, scale = "free_x", space = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray10") +
  geom_jitter(aes(y = log(mean.ratio), x = num.hab, colour = first.obs,
                 shape = abg.dna, size = as.character(quan.bin))) +
  scale_colour_manual(values = colvec,
                      name = "First \ndetected in") +
  scale_shape_manual(values = c(21,23,25),
                     name = "Local DNA\nabundance") +
  scale_size_manual(values = c(1,2,5,6), labels = c("< Q1","Q1 - Median",
                                                    "Median - Q3", "> Q3"),
                    name = "IQR of \n% RNA reads")
  


plot.ab <- plot.ab %>%
  mutate(abg.dna = fct_relevel(abg.dna, "Rare","Moderate","Abundant"))

ggplot(plot.ab) +
  theme_cust("pubr")+
  facet_grid(first.obs~Season, scale = "free_x", space = "free") +
  geom_line(aes(x = num.hab, y = log(perc.otu.dna.ab / perc.otu.rna.ab), colour = abg.dna))

ggplot(plot.ab) +
  theme_cust("pubr")+
  facet_grid(first.obs~Season, scale = "free_x", space = "free") +
  scale_fill_viridis_d(direction = 1) +
  scale_alpha_manual(values = c(1,0.5,0.5)) +
  geom_line(aes(x = num.hab, y = perc.reads.rna.ab * -1, colour = abg.dna)) +
  geom_line(aes(x = num.hab, y = perc.reads.dna.ab, colour = abg.dna)) +
  geom_ribbon(aes(x = num.hab,
                  ymin =  perc.reads.rna.ab * -1,
                  ymax = perc.reads.dna.ab, fill = abg.dna, alpha = abg.dna)) +
  geom_point(data = plot.ab[(first.obs == "Downriver" & 
                               (Season == "Spring" | Season == "Autumn")) |
                              first.obs == "Estuary",],
             aes(x = num.hab, y = perc.reads.dna.ab, colour = abg.dna), size = 0.8) +
  geom_point(data = plot.ab[(first.obs == "Downriver" & 
                               (Season == "Spring" | Season == "Autumn")) |
                              first.obs == "Estuary",],
             aes(x = num.hab, y = perc.reads.rna.ab * -1, colour = abg.dna), size = 0.8) +
  scale_colour_viridis_d(direction = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray10")

ggplot(plot.ab) +
  theme_cust("pubr")+
  facet_grid(first.obs~Season, scale = "free_x", space = "free") +
  scale_fill_viridis_d(direction = -1) +
  scale_alpha_manual(values = c(0.3,0.5,0.8)) +
  geom_line(aes(x = num.hab, y = perc.otu.rna.ab * -1, colour = abg.dna)) +
  geom_line(aes(x = num.hab, y = perc.otu.dna.ab, colour = abg.dna)) +
  geom_ribbon(aes(x = num.hab,
                  ymin =  perc.otu.rna.ab * -1,
                  ymax = perc.otu.dna.ab, fill = abg.dna, alpha = abg.dna)) +
  geom_point(data = plot.ab[(first.obs == "Downriver" & 
                               (Season == "Spring" | Season == "Autumn")) |
                              first.obs == "Estuary",],
             aes(x = num.hab, y = perc.otu.dna.ab, colour = abg.dna), size = 0.8) +
  geom_point(data = plot.ab[(first.obs == "Downriver" & 
                               (Season == "Spring" | Season == "Autumn")) |
                              first.obs == "Estuary",],
             aes(x = num.hab, y = perc.otu.rna.ab * -1, colour = abg.dna), size = 0.8) +
  scale_colour_viridis_d(direction = -1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray10")
# geom_line(aes(x = num.hab, y = perc.otu.dna.ab, 
#colour = abg.dna)) +
#  geom_line(aes(x = num.hab, y = perc.otu.rna.ab * -1, 
                colour = abg.dna)) 
# merge with original percentage dataframe

inac[ac, c("otu.n.all.rna", "sum.reads.all.rna","otu.n.rna","sum.reads.rna") :=
       list(i.otu.n.all.rna, i.sum.reads.all.rna, i.otu.n, i.sum.reads),
     on = .(first.obs, sample.type.year, Season)] #DR.names

# calculate percentages
inac[,c("perc.otu", "perc.reads",
        "perc.otu.rna", "perc.reads.rna") := list(otu.n * 100 / otu.n.all.dna,
                                                  sum.reads * 100 / sum.reads.all.dna,
                                                  otu.n.rna * 100 / otu.n.all.rna,
                                                  sum.reads.rna * 100 / sum.reads.all.rna)]
# Sanity check
inac[,.(sum.otu = sum(perc.otu),
        sum.reads = sum(perc.reads),
        sum.otu.rna = sum(perc.otu.rna),
        sum.reads.rna = sum(perc.reads.rna)),, by =.(sample.type.year, Season)]

ratio.df <- first.df[DNA > 0 & RNA > 0, .(mean = mean(ratio, na.rm = T),
                                        ci = conf.int(ratio)),
                   by = .(sample.type.year, Season, spab.names, first.obs)]

ggplot(ratio.df) +
  theme_cust("pubr") +
  facet_grid(spab.names~Season, scale = "free_x", space = "free") +
  geom_jitter(aes(x = log(mean), fill = first.obs,
                  y= sample.type.year), shape = 21,
              position = position_dodge(0.7)) +
  geom_errorbar(aes(xmin = log(mean) - log(ci),
                xmax = log(mean) + log(ci),
                y= sample.type.year, colour = first.obs),
                position = position_dodge(0.7)) +
  geom_vline(xintercept = 0, linetype = "dashed")

geom_smooth(aes(y = log(mean / mean.ac), x = sample.type.year, colour = abg.dna,
                group = abg.dna
), method = "loess", se = F) +

# calculate the mean abundance of each first observed in habitats along the continuum
m.inac <- first.df[DNA > 0 & RNA > 0, .(mean = mean(DNA, na.rm = T),
                                        ci = conf.int(DNA)),
         by = .(sample.type.year, Season, abg.dna, first.obs)]
m.ac <- first.df[DNA > 0 & RNA > 0, .(mean = mean(RNA, na.rm = T),
                                      ci = conf.int(RNA)),
                   by = .(sample.type.year, Season, abg.dna, first.obs)]

test <- m.inac[m.ac, mean.ac := i.mean, on = .(sample.type.year, Season, first.obs)]

ggplot(test) +
  theme_cust("pubr") +
  facet_grid(first.obs~Season, scale = "free_x", space = "free") +
  #geom_jitter(aes(y = log(mean / mean.ac), fill = abg.dna,
  #                x= sample.type.year), shape = 21,
  #            position = position_dodge(0.7)) +
  geom_smooth(aes(y = log(mean / mean.ac), x = sample.type.year, colour = abg.dna,
                  group = abg.dna
                  ), method = "loess", se = F) +
  geom_hline(yintercept = 0, linetype = "dashed")


sum.abg.dna <- first.df[, .(n.otu = .N,
             sum.reads = sum(DNA, na.rm = T)), by = .(sample.type.year, Season, abg.dna, first.obs)]
sum.abg.rna <- first.df[, .(n.otu = .N *-1,
                            sum.reads = sum(RNA, na.rm = T) * -1), by = .(sample.type.year, Season, abg.dna, first.obs)]

test <- rbind(sum.abg.dna, sum.abg.rna)

ggplot() +
  theme_cust("pubr") +
  facet_grid(sample.type.year~Season, scale = "free_y", space = "free") +
  geom_jitter(data = test[sum.reads > 0,], aes(x = log(sum.reads), y = abg.dna, fill = first.obs), shape = 21,
              position = position_dodge(0.7)) +
  geom_jitter(data = test[sum.reads < 0,], aes(x = log(sum.reads *-1)*-1, y = abg.dna, fill = first.obs), shape = 21,
              position = position_dodge(0.7)) +
  geom_vline(xintercept = 0, linetype = "dashed")

ggplot(sum.abg, aes(x = abg.dna, y = log(sum.reads))) +
  geom_col(aes(fill = first.obs)) +
  facet_grid(Season~sample.type.year)

pdf <- first.df[, .(mean = mean(DNA, na.rm = T),
                    sd = sd(DNA, na.rm = T)),  by = .(sample.type.year, Season, first.obs, spab.names)]

ggplot(pdf, aes(x = log(mean), y = sample.type.year)) +
  facet_grid(spab.names~Season) +
  geom_point(aes(fill = first.obs), shape = 21)

# calculate number of observations per habitat type and season = Total
inac.all <- first.df[DNA > 0 & RNA == 0,
                     .(otu.n.all.dna = .N,
                       sum.reads.all.dna = sum(DNA, na.rm = T)), by = .(sample.type.year, Season)] #DR.names
ac.all <- first.df[RNA > 0,
                   .(otu.n.all.rna = .N,
                     sum.reads.all.rna = sum(RNA, na.rm = T)), by = .(sample.type.year, Season)]  # DR.names

inac <- first.df[DNA > 0 & RNA == 0, .(otu.n = .N,
                               sum.reads = sum(DNA, na.rm = T)),
         by = .(first.obs, abg.dna, sample.type.year, Season)]

ac <- first.df[RNA > 0, .(otu.n.rna = .N,
                            sum.reads.rna = sum(RNA, na.rm = T)),
                 by = .(first.obs, abg.dna, sample.type.year, Season)]

# merge together
inac[inac.all, c("otu.n.all.dna","sum.reads.all.dna") :=
       list(i.otu.n.all.dna, i.sum.reads.all.dna), on = .(sample.type.year, Season)] #DR.names
ac[ac.all, c("otu.n.all.rna","sum.reads.all.rna") :=
     list(i.otu.n.all.rna, i.sum.reads.all.rna), on = .(sample.type.year, Season)]  #DR.names

inac[ac, c("otu.n.all.rna", "sum.reads.all.rna","otu.n.rna","sum.reads.rna") :=
       list(i.otu.n.all.rna, i.sum.reads.all.rna, i.otu.n.rna, i.sum.reads.rna),
     on = .(first.obs, abg.dna, sample.type.year, Season)] #DR.names

# calculate percentages
inac[,c("perc.otu", "perc.reads",
        "perc.otu.rna", "perc.reads.rna") := list(otu.n * 100 / otu.n.all.dna,
                                                  sum.reads * 100 / sum.reads.all.dna,
                                                  otu.n.rna * 100 / otu.n.all.rna,
                                                  sum.reads.rna * 100 / sum.reads.all.rna)]
# Sanity check
inac[,.(sum.otu = sum(perc.otu, na.rm = T),
        sum.reads = sum(perc.reads, na.rm = T),
        sum.otu.rna = sum(perc.otu.rna, na.rm = T),
        sum.reads.rna = sum(perc.reads.rna, na.rm = T)),, by =.(sample.type.year, Season)]

# re-assign factor levels
plot.df <- inac
plot.df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                  "Upriver", "Reservoirs", "Downriver",
                                                                  "Estuary"))]

d <- ggplot(plot.df, aes(x = first.obs)) +
  theme_cust("pubr") +
  facet_grid(sample.type.year~Season, space = "free") +
  geom_col(aes(y = perc.otu, fill = abg.dna), colour = "gray20")
  #facet_grid(.~Season, scale = "free_x", space = "free") 
  #theme(axis.text.x = element_blank(),
  #      axis.title.x = element_blank(),
  #      strip.background = element_rect(fill = "gray20"),
  #      strip.text = element_text(colour = "white", size = 10)) +
  #labs(x = "Habitat type", y = "% of OTUs \nin inactive fraction") #+
  #scale_fill_manual(values = colvec[names(colvec) %in% plot.df$sample.type.year])

r <- ggplot(plot.df, aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = perc.otu.rna, fill = first.obs), colour = "gray20") +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  labs(x = "Habitat type", y = "% of OTUs \nin reactive fraction") +
  scale_fill_manual(values = colvec[names(colvec) %in% plot.df$sample.type.year])

dr <- ggplot(plot.df, aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = perc.reads, fill = first.obs), colour = "gray20") +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 10)) +
  labs(x = "Habitat type", y = "% of CSS reads \nin inactive fraction") +
  scale_fill_manual(values = colvec[names(colvec) %in% plot.df$sample.type.year],
                    name = "First observed in")


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
