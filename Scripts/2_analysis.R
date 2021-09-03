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
otu.tab <- select_newest("./Output", "201520162017_fin_css_otu99_table_paper1_")
otu.tab <- as.matrix(read.csv(
  paste0("./Output/", otu.tab),
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))
otu.tab[1:4,1:3]
#otu.tab <- otu.tab[, 1:1000]

# row orders need to match between tax.tab and otu.tab
otu.tab <- otu.tab[order(row.names(otu.tab)),]

# Taxonomy table
tax.tab <- select_newest("./Output", "201520162017_tax_otu99_table_paper1_")
tax.tab <-
  as.matrix(read.csv(
    paste0("./Output/", tax.tab),
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ))
tax.tab[1:4,]
# orders need to match between tax.tab and otu.tab
tax.tab <- tax.tab[order(row.names(tax.tab)),]
#tax.tab <- tax.tab[,1:1000]

# Meta data
met.df <-
  select_newest(path = "./Output", file.pattern = "201520162017_meta_otu99_paper1_")
met.df <-
  read.csv(
    paste0("./Output/", met.df),
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  )

met.df[1:4,]
# Checked by PCoA THW1D is a soilwater sample (actually Hyporheic)
met.df[met.df$seq_name == "THW1D", ]$sample.type.year <- "Hyporheic"
# Some downriver samples were wrongly assigned to "Downriver"
met.df[grep("SWLR", met.df$seq_name),]$sample.type.year <- "Soilwater"
met.df[grep("SLR", met.df$seq_name),]$sample.type.year <- "Soil"


# merge some sample types
met.df$sample.type.year <- factor(met.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheic", 
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                                      "Marine"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Soilwater", 
                                             "Groundwater","Stream", "Tributary",
                                             "Riverine Lakes", "Headwater Ponds", "Lake", "Lake",
                                             "Upriver",# "RO3", "RO2", "RO1","Deep",
                                             "Reservoirs","Reservoirs", "Reservoirs","Reservoirs",
                                             "Downriver",
                                             "Estuary"))

met.df$Season <- factor(met.df$Season, levels = c("Spring","Summer","Autumn"),
                        labels = c("Spring","Summer","Autumn"))

met.df$dna_type <- factor(met.df$dna_type, levels = c("DNA","cDNA"),
                        labels = c("DNA","RNA"))

# Construct phyloseq object
pb <- phyloseq(otu_table(otu.tab, taxa_are_rows = F),
               sample_data(met.df),
               tax_table(tax.tab))

#pb <- prune_taxa(!taxa_sums(pb) == 0, pb)
pb <- prune_samples(!sample_sums(pb) == 0, pb)

# remove THW1D - no meta data
pb <- subset_samples(pb, seq_name != "THW1D")

# check rarefaction curve
#t <- round(otu_table(pb), digits = 0)
#rarecurve(t, step = nrow(t), cex = 0.5, label = F)

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
  "orchid", #"#471063FF", #Reservoir,
  #'salmon',
  #"pink",
  #'navy',
  "#375A8CFF", #Downriver,
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
# custom theme
# Custom theme with smaller legends
theme_cust <- function(base_theme = "bw", base_size = 11, half_line = base_size/2,
                       border = F){
  if(base_theme == "bw"){
    t <- theme_bw(base_size = base_size)
  }
  
  if(base_theme == "pubr"){
    if(border == T){
      t <- theme_pubr(base_size = base_size, border = T)
    } else {
      t <- theme_pubr(base_size = base_size)
    }
  }
  
  t %+replace% theme(legend.position = "right", 
                     legend.spacing = unit(half_line / 2, "pt"), # spacing between legends
                     legend.key.size = unit(0.8, "lines"), # size of legend symbol box
                     legend.box.spacing = unit(1.5 * half_line, "pt"), # spacing between legend and plot
                     legend.text = element_text(size = unit(base_size - 4, "pt")), 
                     legend.title = element_text(size = unit(base_size - 3, "pt"), hjust = 0)
  )
  
}

#theme_set(theme_pubr())

# only save lat, long, of DNA for GIS
# gis <- met.df %>% 
#   filter(dna_type == 'DNA') %>%
#   select(dr_match_name, lat, long) %>%
#   unique()
# write.table(gis, "./Output/gis_samplepnts_miccom.csv", row.names = F, sep = ",")


# Phantom taxa ---------------------------------------------------------------------------------------
# otu.tab <- setDT(as.data.frame(otu_table(pb), strings.As.Factors = F), keep.rownames = "Sample")
# # melt OTU table into long format
# otu.tab <- melt.data.table(
#   otu.tab,
#   id.vars = "Sample",
#   measure.vars = patterns("^OTU_"),
#   variable.name = "OTU",
#   value.name = "reads"
# )
# 
# # Join OTU and meta table
# sumdf <- left_join(
#   otu.tab,
#   met.df[met.df$seq_name %in% otu.tab$Sample,] %>%
#     dplyr::select(
#       seq_name,
#       dr_match_name,
#       year,
#       Season,
#       dna_type,
#       replicate,
#       sample.type,
#       sample.type.year,
#     ),
#   by = c("Sample" = "seq_name")
# )
# 
# # set back to data.table, order data.table by catchment.area
# setDT(sumdf)
# 
# # check for duplicates
# any(duplicated(sumdf[dna_type == "DNA" & OTU == "OTU_10002",]$dr_match_name) == T) # no
# any(duplicated(sumdf[dna_type == "RNA" & OTU == "OTU_10002",]$dr_match_name) == T)
# 
# # CSS reads creates 0.5 reads overwrite all with 1 (integer, makes them 0)
# # Everything below 1 and above 0 is 1
# range(sumdf[reads > 0,]$reads, na.rm = T)
# sumdf[reads > 0 & reads <= 1, reads := 1] # 0.1:0.9 = 1
# 
# # round to counts
# sumdf[, reads := round(reads, digits = 0)]
# range(sumdf[reads > 0,]$reads, na.rm = T)
# 
# # cast into wide format
# castdf <- dcast(sumdf, dr_match_name + OTU ~ dna_type, value.var = "reads")
# 
# # fill NAs with 0, those are taxa not found in RNA or DNA
# castdf[is.na(DNA), DNA := 0]
# castdf[is.na(RNA), RNA := 0]
# 
# # how many samples have an OTU with RNA > 0 and DNA == 0? (phantom taxa)
# castdf[DNA > 0 | RNA > 0, n.all := .N, .(dr_match_name)] # number of any observations above 0 by sample
# temp <- castdf[RNA > 0 & DNA == 0,]
# temp <- temp[, .(n = .N, n.all = unique(n.all)), .(dr_match_name)][, prop := n * 100 / n.all]
# 
# setDT(met.df)
# # which samples have more than 50% phantom taxa
# View(met.df[dr_match_name %in% temp[prop != 100 & prop > 50,]$dr_match_name,] %>% select(dr_match_name,dna_type,LibrarySize))
# 
# # overwrite all DNA observations == 0, where RNA > 0 with 1
# castdf[RNA > 0 & DNA == 0, DNA := 1]
# castdf[RNA > 0 & DNA == 0,] # check if overwrite was successful
# 
# # format back to long format
# temp <- melt.data.table(castdf,
#                 id.vars = c("dr_match_name","OTU"),
#                 measure.vars = c("DNA","RNA"),
#                 variable.name = "dna_type",
#                 value.name = "reads")
# 
# # original and casted data frame have not the same length
# # possible we added rows for samples that didn't have any RNA, but with RNA = 0
# nrow(temp) == nrow(sumdf)
# 
# # do row merge/overwrite
# # we have dealt with the replicates already in the earlier script
# mer <- sumdf[temp, cor.reads := i.reads, on = .(dr_match_name, OTU, dna_type)]
# 
# # cast into wide format
# fin.df <- dcast(mer, Sample ~ OTU, value.var = "cor.reads")
# # make matrix and give row.names
# fin.df <- as.matrix(setDF(fin.df, rownames = fin.df$Sample)[,-1]) # remove first row with "Sample" move into rownames
# 
# # same dimensions as original df?
# dim(fin.df)
# dim(otu_table(pb)) # yes
# 
# write.table(fin.df, paste0("./Output/201520162017_fin_css_otu99_phantomcor_paper1_",Sys.Date(),".csv"),
#             sep = ";", dec = ".", row.names = T)
# rm(temp, castdf, sumdf, otu.tab)


# 3. Analysis ----------------------------------------------------------------------------------

# read in
fin.df <- select_newest("./Output", "201520162017_fin_css_otu99_phantomcor_paper1_")
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

# remove THW1D - no meta data
pb <- subset_samples(pb, seq_name != "THW1D")
pb <- prune_taxa(taxa_sums(pb) != 0, pb)
pb <- prune_samples(sample_sums(pb) != 0, pb) # remove samples with no reads

# Summary stats
dim(otu_mat(pb))
df <- sample_df(pb); setDT(df)
df[which.min(LibrarySize),]
df[which.max(LibrarySize),]
means <- df[, .(mean.lib = mean(LibrarySize)), by = .(sample.type.year)]
#means %>% arrange(mean.lib)

rm(means, df)

## Figure 1: DNA -----------------------------------------------------------------------------------
# Q: Is the microbial assemblage different between habitat types and seasons?

#--------------------------#
#- Multivariate analysis -#
#--------------------------#
# We are dealing with large environmental gradients, thus we expect a high proportion of zeros
# Therefore, we need to select an asymmetrical similarity distance coefficient
# If we're dealing with samples from fairly homogeneous environmental conditions (short envir. gradients)
# and we expect few zeros and symmetric association coefficients we can use the Euclidean distance.
# Transformations can be used to make asymmetric species distributions more symmetric so that e.g.
# Euclidean distances can be used
# E.g. Logarithm transformation log2(x + 1) (in case of microbial data) can be used to make asymmetric species
# distributions less asymmetric

# Several methods were compared:
# PCA
# CA
# PCoA
# NMDS - no convergence

# PCoA was chosen for all ordinations for consistency

# subset only DNA samples
  dna <- subset_samples(pb, dna_type == "DNA")
  
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

ncol(pb.mat) # 16322 OTUs
nrow(pb.mat) # 389 Samples

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray) # 390 registers
# plot with custom function (= made to avoid repetitive code)
# custom function is in ./Functions/custom_fun.R
# extract scores and variance explained
dna.pcoa <- plot_pcoa(pb.bray.pcoa, physeq = dna,  colours = colvec, output = T)

# main PCoA
p <- dna.pcoa$plot + theme_cust(base_theme = "pubr",
                                border = T) + guides(colour = "none")

# zoom into terrestrial part
# different colouring by Season, habitat type as shape
# zoom.terr <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
#     theme_cust() +
#     geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
#     geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
#     geom_point(size = 2.5, alpha = 0.2, fill = "gray20") +
#     new_scale_fill() +
#     geom_point(data = dna.pcoa$df[dna.pcoa$df$sample.type.year == "Soil" | 
#                                     dna.pcoa$df$sample.type.year == "Soilwater" |
#                                     dna.pcoa$df$sample.type.year == "Sediment",], 
#                aes(x = Axis.1, y = Axis.2, fill = Season, shape = sample.type.year), size = 3) +
#     scale_fill_viridis_d(option = "cividis", name = "Season") +
#     scale_shape_manual(values = c(21,23, 25), "Habitat type") +
#     coord_cartesian(ylim = c(-0.11, 0.04), xlim = c(-0.39, -0.17), expand = F) + 
#     labs(x = paste("PC1"), 
#          y = paste("PC2")) +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()) +
#   guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
#          alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
#          fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
# 
# # zoom into estuaries
# # different colouring by distance from mouth, shape is Season
# zoom.est <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
#   theme_cust() +
#   geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
#   geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
#   geom_point(size = 2.5, alpha = 0.2, fill = "gray20") +
#   new_scale_fill() +
#   geom_point(data = dna.pcoa$df[dna.pcoa$df$sample.type.year == "Estuary",], 
#              aes(x = Axis.1, y = Axis.2, fill = abs(distance.from.mouth), shape = Season), size = 3) +
#   scale_fill_viridis_c(name = "Distance from \nmouth (km)", direction = -1) +
#   scale_shape_manual(values = c(21,23,25)) +
#   scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
#   coord_cartesian(ylim = c(-0.17, 0.3), xlim = c(-0.1, 0.28), expand = F) + # ensure aspect ratio
#   labs(x = paste("PC1"), 
#        y = paste("PC2")) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
  #guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
  #       alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
  #       fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

# zoom into rivers (upriver vs. downriver)
# different colouring by distance from mouth, shape is Season
# zoom.river <- ggplot(dna.pcoa$df, aes(x = Axis.1, y = Axis.2)) +
#   theme_cust() +
#   geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
#   geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
#   geom_point(size = 2.5, alpha = 0.2, shape = 21, fill = "gray20") +
#   #scale_fill_manual(values = colvec, name = "Habitat Type") +
#   new_scale_fill() +
#   geom_point(data = dna.pcoa$df[(dna.pcoa$df$sample.type.year == "Downriver" |
#                                   dna.pcoa$df$sample.type.year == "Upriver") & Season == "Spring",],
#              aes(x = Axis.1, y = Axis.2,
#                  fill = distance.from.mouth, shape = sample.type.year), size = 3) +
#   stat_ellipse(data = subset(dna.pcoa$df, (sample.type.year == "Upriver" | 
#                                              sample.type.year == "Downriver")  &
#                                Season == "Spring"),
#                  aes(x = Axis.1, y = Axis.2, group = year), 
#                linetype = "dashed", colour = "grey50") +
#   #stat_ellipse(data = subset(dna.pcoa$df, (sample.type.year == "Upriver" | sample.type.year == "Downriver")),
#   #               aes(x = Axis.1, y = Axis.2, group = paste(year,Season), linetype = as.character(year))) +
#   annotate(geom = "text", x = c(0.2,0.2), y = c(-0.3,-0.19), label = c("2015", "2016"), 
#            size = 3, colour = "grey50") +
#   scale_fill_continuous(type = "viridis", name = "Distance from \nmouth (km)", direction = -1) +
#   #scale_fill_viridis_b(name = "Distance from mouth", direction = -1) +
#   scale_shape_manual(values = c(21,23), name = "Habitat Type") +
#   #coord_cartesian(xlim = c(0.03, 0.32), ylim = c(-0.2,0.3), expand = F) + # ensure aspect ratio
#   labs(x = paste("PC1"), 
#        y = paste("PC2")) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
#   #guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
#   #       alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
#   #       fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))


# annotate collage boxes within main plot
# p <- p + 
#   annotate(geom = "rect", xmin = -0.39, ymin = -0.11, 
#            xmax = -0.17, ymax = 0.04, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
#   annotate(geom = "rect", xmin = -0.1, xmax = 0.28,
#            ymin = - 0.17, ymax = 0.3, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
#   #annotate(geom = "rect", xmin = -0.1, xmax = 0.28,
#   #         ymin = - 0.12, ymax = 0.245, fill = NA, colour = "grey50", linetype = "dashed", alpha = 0.8) +
#   annotate(geom = "text", x = c(-0.41,-0.12), 
#            y = c(0.035,0.295), label = c("b","c"), alpha = 0.7, colour = "grey50") #"d", x = -0.085, y = 0.23
# 
# # combine all plots into one
# (collage <- ggarrange(p, 
#           ggarrange(zoom.terr,  zoom.est, ncol = 2, labels = c("b","c"), align = "hv"), #zoom.river,
#           nrow = 2, labels = c("a"), heights = c(0.6, 0.4))
# )

# save
ggsave(paste0("./Figures/Submission/Figure2_new.tiff"), p,
       width = 14, height = 10, unit = "cm", dpi = 300)

ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_SampleType_n.tiff"), p,
       width = 12, height = 10, unit = "cm", dpi = 300)
ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_SampleType_n.png"),  p,
       width = 12, height = 10, unit = "cm", dpi = 300)

ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_collage_n.tiff"), collage,
       width = 18, height = 15, unit = "cm", dpi = 300)
ggsave(paste0("./Figures/Final/PCoA_hellin_DNA_collage_n.png"),  collage,
       width = 18, height = 15, unit = "cm", dpi = 300)

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

# Results
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# sample.type.year  12    45.490  3.7909  17.500 0.34533  1e-04 ***
#   Season             2     5.222  2.6110  12.053 0.03964  1e-04 ***
#   Residuals        374    81.017  0.2166         0.61503           
# Total            388   131.729                 1.00000  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

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
# habitat alone = F[12] = 38.875, p = < 2.2e-16
# Response: Distances
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Groups     12 1.3397 0.111638  38.875 < 2.2e-16 ***
#   Residuals 376 1.0798 0.002872    

# season alone = F[2] = 58.361, p = < 2.2e-16
# Response: Distances
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Groups      2 0.89849 0.44925  58.361 < 2.2e-16 ***
#   Residuals 386 2.97133 0.00770  

# combined = F[27] = 20.717 , p = < 2.2e-16
# Response: Distances
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Groups     27 2.3975 0.088797  20.717 < 2.2e-16 ***
#   Residuals 358 1.5344 0.004286     

options(scipen = 999)
(perm.dna <- data.frame(Dataset = "DNA", Group = c("Habitat","Season", "Combined"),
                       df = c(perm.mod$aov.tab$Df[1:2],NA),
                       F = c(round(perm.mod$aov.tab$F.Model[1:2], digits = 2), NA),
                       Rsq = c(round(perm.mod$aov.tab$R2[1:2], digits = 2), NA),
                       pval = c(paste( "<", perm.mod$aov.tab$`Pr(>F)`[1:2]), NA),
                       df.p = c(perm1$Df[1], perm2$Df[1], perm3$Df[1]),
                       F.p = c(round(perm1$`F value`[1], digits = 2),
                               round(perm2$`F value`[1], digits = 2),
                               round(perm3$`F value`[1], digits = 2)),
                       pval.p = c(abbrev.p(perm1$`Pr(>F)`[1])[1],
                                  abbrev.p(perm2$`Pr(>F)`[1])[1],
                                  abbrev.p(perm3$`Pr(>F)`[1])[1])))

# look at dispersion
mod # habitats
mod2 # seasons

# individual comparisons (exploring differences in dispersion among groups)
tuk <-TukeyHSD(mod2)
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
#pb.mat <- log1p(pb.mat)
#pb.mat <- decostand(pb.mat, "hellinger")
# creates strong horse-shoe

# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray)
# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)


ncol(pb.mat) # 16322 OTUs
nrow(pb.mat) # 590 registers

ncol(pb.bray.pcoa$vectors) # 589 axes

# Sanity check
# back calculate BC dissimilarity matrix from PCoA
#coor.dist <- as.matrix(dist(pb.bray.pcoa$vectors)) # same as original pb.bray matrix

# plot with custom function
all.pcoa <- plot_pcoa(pb.bray.pcoa, physeq = pb, axes = c(1:3), colours = colvec, output = T)
# plot second and third axes
pcoa.23 <- plot_pcoa(pb.bray.pcoa, physeq = pb, plot.axes = c(3,2), colours = colvec, output = T)

# extract legends
# add grey box to colour legend
#leg1 <- get_legend(all.pcoa$plot + theme_cust(base_theme = "pubr") + 
#                     theme(legend.key = element_rect(fill = "gray50")) +
#                                           guides(colour = guide_legend(order = 3, size = 2.5,
#                                                                 override.aes = list(fill = "grey20",
#                                                                                     shape = 21)),
#                                                  fill = F, shape = F))
#leg2 <- get_legend(all.pcoa$plot + theme_cust(base_theme = "pubr") + guides(colour = F) +
#                     theme(legend.margin = margin(3,5.5,3,5.5)))
# change layout to remove gap between legends
#legs <- gtable_rbind(leg2, leg1)
#legs$layout[4,c(1,3)] <- c(9,9)

find_hull <- function(x, axes){x[chull(x[,paste0("Axis.",axes[1])], x[,paste0("Axis.",axes[2])]),]}
hulls <- ddply(all.pcoa$df, "dna_type", find_hull, axes = c(1,2))
                           
pcoa.plot <- ggplot(all.pcoa$df, aes(x =Axis.1, 
                                           y = Axis.2)) +
  theme_cust(base_theme = "pubr",
             border = T) +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
  scale_fill_manual(values = colvec[names(colvec) %in% as.character(levels(all.pcoa$df$sample.type.year))],
                      name = "Habitat Type") +
  geom_point(aes(fill = sample.type.year, shape = Season), size=2.5) + #colour = dna_type,
  geom_polygon(data = hulls, alpha = 0, aes(linetype = dna_type), fill = "white", colour = "black") +
  scale_linetype_manual(values = c("solid","dotted"), name = "Nucleic Acid Type") +
  coord_fixed(1) + # ensure aspect ratio
  scale_shape_manual(values = c(21,23,25)) +
  #scale_size_manual(values = c(2.5, 2.6), name = "Nucleic Acid \nType") +
  labs(x = paste0("PC1 (", round(all.pcoa$var[1,"var"],
                                                  digits = 1),"%)"), 
       y = paste0("PC2 (", round(all.pcoa$var[2,"var"],
                                                  digits = 1),"%)")) +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)),
         size = FALSE)

# 2nd dimension
hulls <- ddply(pcoa.23$df, "dna_type", find_hull, axes = c(3,2))

pcoa.plot2 <- ggplot(pcoa.23$df, aes(x =Axis.3, 
                                     y = Axis.2)) +
  theme_cust(base_theme = "pubr",
             border = T) +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
  scale_fill_manual(values = colvec[names(colvec) %in% as.character(levels(pcoa.23$df$sample.type.year))],
                    name = "Habitat Type") +
  geom_point(aes(fill = sample.type.year, shape = Season), size=2.5) + #colour = dna_type,
  scale_linetype_manual(values = c("solid","dotted"), name = "Nucleic Acid Type") +
  geom_polygon(data = hulls, alpha = 0, aes(linetype = dna_type), fill = "white", colour = "black") +
  coord_fixed(1) + # ensure aspect ratio
  scale_shape_manual(values = c(21,23,25)) +
  
  #scale_size_manual(values = c(2.5, 2.6), name = "Nucleic Acid \nType") +
  labs(x = paste0("PC1 (", round(pcoa.23$var[1,"var"],
                                 digits = 1),"%)"), 
       y = paste0("PC2 (", round(pcoa.23$var[2,"var"],
                                 digits = 1),"%)")) +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)),
         size = FALSE)
  
pcoa.plot2 <- pcoa.23$plot +  + theme(legend.position = "none")

(p <- ggarrange(pcoa.plot, pcoa.plot2, ncol = 2,
               align = "hv", labels = "auto", common.legend = T, # legend.grob = legs,
               legend = "right"))

# save
ggsave(paste0("./Figures/Submission/Figure3_new.tiff"), p,
       width = 21, height = 11, unit = "cm", dpi = 300)

ggsave(paste0("./Figures/Final/PCoA_all_SampleType_n.tiff"), p,
       width = 20, height = 11, unit = "cm", dpi = 300)
ggsave(paste0("./Figures/Final/PCoA_all_SampleType_n.png"),  p,
       width = 20, height = 11, unit = "cm", dpi = 300)

# check screeplot
ggplot(pb.bray.pcoa$values[pb.bray.pcoa$values$Eigenvalues > 1,],
       aes(x = as.numeric(row.names(pb.bray.pcoa$values[pb.bray.pcoa$values$Eigenvalues >= 1,])), y = Eigenvalues)) +
         geom_col()

# how many axes to have 75% of variance captured?
nrow(pb.bray.pcoa$values[pb.bray.pcoa$values$Cumul_eig <= 0.75,])
# 192 axes

# explore other axes
#for(i in 2:ncol(pb.bray.pcoa$vectors)){
#  temp <- plot_pcoa(pb.bray.pcoa, physeq = pb, plot.axes = c(1,i), colours = colvec, output = T)
#  ggsave(paste0("./Figures/PCoA/PCoA_DR_ax1_",i,".png"), temp$plot,
#         width = 12, height = 10, unit = "cm")
#}


# PERMANOVA -------------------------------------------------------------------------------------
# sensitive towards unbalanced sampling designs = bias.adjust
ord.df <- all.pcoa[["df"]]
ord.df$groups <- paste(ord.df$sample.type.year, ord.df$Season, ord.df$dna_type, sep = "_")

setDT(ord.df)
# Test for significant difference between factors
perm.mod <- adonis(pb.mat ~ sample.type.year + Season + dna_type, 
                   permutations = 9999, data = ord.df, sqrt.dist = T, method = "bray",
                   parallel = cl)
#-- results
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# sample.type.year  12    60.088  5.0073  19.977 0.27534  1e-04 ***
#   Season             2     7.842  3.9208  15.642 0.03593  1e-04 ***
#   dna_type           1     6.426  6.4263  25.637 0.02945  1e-04 ***
#   Residuals        574   143.878  0.2507         0.65929           
# Total            589   218.234                 1.00000  

# habitat = F[12] = 19.977, R^2 = 0.27534, p = 1e-04
# season = F[2] = 15.642, R^2 = 0.03593, p = 1e-04
# nucleic acid type = F[1] = 25.637, R^2 = 0.02945, p = 1e-04
# iterations = 9999

# calculate multivariate dispersions
mod <- betadisper(pb.bray, 
                  group = ord.df$sample.type.year,
                  bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
mod2 <- betadisper(pb.bray, 
                   group = ord.df$Season,
                   bias.adjust = T) # do not set sqrt.dist = T, as pb.bray is already sqrt transformed
mod3 <- betadisper(pb.bray, 
                   group = ord.df$dna_type,
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
# Response: Distances
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Groups     12 1.3880 0.115666  28.881 < 2.2e-16 ***
#   Residuals 577 2.3108 0.004005     
# habitat alone = F[12] = 28.881, p =  < 2.2e-16

# Response: Distances
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Groups      2 0.8522 0.42611  62.033 < 2.2e-16 ***
#   Residuals 587 4.0321 0.00687  
# season alone = F[2] = 62.033 p = < 2.2e-16

# Response: Distances
# Df Sum Sq   Mean Sq F value Pr(>F)
# Groups      1 0.0090 0.0090468  1.1992 0.2739
# Residuals 588 4.4357 0.0075438   
# nucleic acid type alone = F[1] = 1.1992 p = 0.2739
# 
# Response: Distances
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Groups     49 3.8349 0.078263  9.6002 < 2.2e-16 ***
#   Residuals 534 4.3533 0.008152    
# combined = F[49] = 9.6002, p = < 2.2e-16

(perm.dr <- data.frame(Dataset = "DNA and RNA", Group = c("Habitat","Season", "Nucleic Acid Type","Combined"),
                        df = c(perm.mod$aov.tab$Df[1:3],NA),
                        F = c(round(perm.mod$aov.tab$F.Model[1:3], digits = 2), NA),
                        Rsq = c(round(perm.mod$aov.tab$R2[1:3], digits = 2), NA),
                        pval = c(paste( "<", perm.mod$aov.tab$`Pr(>F)`[1:3]), NA),
                        df.p = c(perm1$Df[1], perm2$Df[1], perm3$Df[1],perm4$Df[1]),
                        F.p = c(round(perm1$`F value`[1], digits = 2),
                                round(perm2$`F value`[1], digits = 2),
                                round(perm3$`F value`[1], digits = 2),
                                round(perm4$`F value`[1], digits = 2)),
                        pval.p = c(abbrev.p(perm1$`Pr(>F)`[1])[1],
                                   abbrev.p(perm2$`Pr(>F)`[1])[1],
                                   abbrev.p(perm3$`Pr(>F)`[1])[3],
                                   abbrev.p(perm4$`Pr(>F)`[1])[1])))

tuk <-TukeyHSD(mod)
t <- data.frame(Groups = rownames(tuk$group),unlist(tuk$group))
t <- t %>% separate(Groups, into = c("habitat","season","nucacid",
                                "habitat_y","season_y","nucacid_y"), sep = "_|-")
row.names(t) <- NULL
# all pair-wise comparisons that are significant
t <- t[t$p.adj < 0.05,]

# combine tables into one
permanova.df <- rbind(perm.dna, perm.dr)

# replace NAs with ""
permanova.df <- permanova.df %>% replace(is.na(.), "")

knitr::kable(permanova.df, "latex", booktabs = T) %>%
  kable_styling(position = "center", full_width = F) %>%
  add_header_above(c(" " = 2, "PERMANOVA" = 4, "PERMDISP" = 3)) %>%
  collapse_rows(columns = 1, valign = "middle")

write.xlsx(permanova.df, file = "./Manuscript/Tables/table1.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE)

# \begin{table}[H]
# \centering
# \begin{tabular}{llllllrrl}
# \toprule
# \multicolumn{2}{c}{ } & \multicolumn{4}{c}{PERMANOVA} & \multicolumn{3}{c}{PERMDISP} \\
# \cmidrule(l{3pt}r{3pt}){3-6} \cmidrule(l{3pt}r{3pt}){7-9}
# Dataset & Group & df & F & Rsq & pval & df.p & F.p & pval.p\\
# \midrule
# & Habitat & 12 & 17.5 & 0.35 & < 0.0001 & 12 & 38.87 & < 0.0001\\
# \cmidrule{2-9}
# & Season & 2 & 12.05 & 0.04 & < 0.0001 & 2 & 58.36 & < 0.0001\\
# \cmidrule{2-9}
# \multirow{-3}{*}{\raggedright\arraybackslash DNA} & Combined &  &  &  &  & 27 & 20.72 & < 0.0001\\
# \cmidrule{1-9}
# & Habitat & 12 & 19.98 & 0.28 & < 0.0001 & 12 & 28.88 & < 0.0001\\
# \cmidrule{2-9}
# & Season & 2 & 15.64 & 0.04 & < 0.0001 & 2 & 62.03 & < 0.0001\\
# \cmidrule{2-9}
# & Nucleic Acid Type & 1 & 25.64 & 0.03 & < 0.0001 & 1 & 1.20 & 0.27\\
# \cmidrule{2-9}
# \multirow{-4}{*}{\raggedright\arraybackslash DNA and RNA} & Combined &  &  &  &  & 49 & 9.60 & < 0.0001\\
# \bottomrule
# \end{tabular}
# \end{table}

## 3D plot ---------------------------------------------------------------------------------------
# Save 3D plot for all PCoA, third axis separates seasons
plot.df <- all.pcoa[["df"]]
p <- plot_ly(type = "scatter", mode = "markers")

p <- plot_ly(plot.df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, color = ~sample.type.year, colors = colvec,
             stroke = ~dna_type, strokes = c("black", "white"),
             size = 5, symbol = ~Season, symbols = c(21,23,25), text = ~Sample)
p <- p %>% add_markers()
p <- p %>% layout(scene = list(xaxis = 
                                 list(title = paste("PC1 [", unique(plot.df$x), "%]")),
                               yaxis = 
                                 list(title = paste("PC2 [", unique(plot.df$y), "%]")),
                               zaxis = 
                                 list(title = paste("PC3 [", unique(plot.df$z), "%]"))))
p

# for(i in seq(0,6.3,by=0.1)){
#   outfile <- paste("plot",round(i,digits=2), sep = "_")
#   cam.zoom = 2
#   ver.angle = 0
#   graph <- p %>%
#     layout(scene=list(camera = list(eye = list(x = cos(i)*cam.zoom,y = sin(i)*cam.zoom, z=0.2),
#                                     center = list(x = 0,
#                                                   y = 0,
#                                                   z = 0
#                                     )
#                       )
#     )
#     )
#   graph
#   
#   
#   cat("Now rendering iteration:", i,"\n")
#   plotly::orca(graph, file = paste(outfile,"png", sep="."),
#                  width = 1200,
#                  height = 1050)
# }

htmlwidgets::saveWidget(as_widget(p), "PCoA_DNARNA_Season_3D.html")

## Figure 4: Distance -----------------------------------------------------------------------------------
# Q: How different are the the DNA and RNA assemblages of the same sample?

# Calculate incidence based dissimilarity
# PCoA with Sorensen, incidence based equivalent of Bray Curtis
pb.mat <- otu_mat(pb)
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
# 204 axes

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
                   sample_df(pb) %>% dplyr::select(dr_match_name, dna_type, year, Season, sample.type.year), 
                   stringsAsFactors = F)

# add meta data for both sample x and sample y
dissim.df <- merge(dissim.df, meta, by.x =  "Sample.x", by.y = "Sample")
dissim.df <- merge(dissim.df, meta %>% select(dna_type, Sample, year, dr_match_name), by.x =  "Sample.y", by.y = "Sample")

# omit all samples of 2015 (no RNA was taken, and sample name strategy changed -> creates duplicates)
dissim.df <- dissim.df[dissim.df$year.x != 2015,]
dissim.df <- dissim.df[dissim.df$year.y != 2015,]

# keep all rows where Sample.x and Sample.y are the same
dissim.df <- dissim.df[dissim.df$dr_match_name.x == dissim.df$dr_match_name.y,]
dissim.df <- dissim.df[!(dissim.df$dna_type.x == dissim.df$dna_type.y),] # omit all distances between same dna_type

dissim.dr <- dissim.df %>% select(Metric, ID = dr_match_name.x, year = year.x, Season, sample.type.year, dist)
setDT(dissim.dr)

# add new column to split plot into main and side panel
dissim.dr[, panels := "main"]
dissim.dr[sample.type.year == "Tributary" |
          sample.type.year == "Lake" |
          sample.type.year == "Riverine \nLakes" |
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
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, year, 
                                                   dna_type, distance.from.mouth, dr_match_name), 
                   stringsAsFactors = F)
data <- merge(pb.scores, meta, by = "Sample")
data$Sample <- as.character(data$Sample)

setDT(data)
setcolorder(data, c("Sample","Metric",
                    "sample.type.year","Season","year", "dna_type","distance.from.mouth","dr_match_name",
                    colnames(data)[!(colnames(data) %in% c("Sample","Metric",
                                                           "sample.type.year","Season","year", "dna_type","distance.from.mouth","dr_match_name"))]))
# melt datatable
temp <- melt.data.table(data, id.vars = c("dr_match_name","dna_type", "Metric"), measure.vars = patterns("^Axis."),
             variable.name = "Axis", value.name = "Coordinates")
temp <- dcast(temp, dr_match_name + Axis + Metric ~ dna_type, value.var = c("Coordinates"))
# remove NAs
temp <- na.omit(temp)

# Calculate distance of all axes
temp[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
temp <- temp[, .(sum.dist = sum(pnt.dist)), by = .(Metric, dr_match_name)] # sum temp axes
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
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, year, 
                                                   dna_type, distance.from.mouth, dr_match_name), 
                   stringsAsFactors = F)
bray.df <- setDT(merge(meta, bray.df, by = "Sample"))
soren.df <- setDT(merge(meta, soren.df, by = "Sample"))

# melt and combine bray-curtis and sorensen results into one data frame
temp.75 <- rbind(melt.data.table(bray.df, id.vars = c("dr_match_name","dna_type", "Metric"), measure.vars = patterns("^Axis."),
                      variable.name = "Axis", value.name = "Coordinates"),
      melt.data.table(soren.df, id.vars = c("dr_match_name","dna_type", "Metric"), measure.vars = patterns("^Axis."),
                      variable.name = "Axis", value.name = "Coordinates"))

temp.75 <- dcast(temp.75, dr_match_name + Axis + Metric ~ dna_type, value.var = c("Coordinates"))
# remove NAs
temp.75 <- na.omit(temp.75)
temp.75[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
temp.75 <- temp.75[, .(sum.dist = sum(pnt.dist)), by = .(Metric, dr_match_name)] # sum temp axes
temp.75 <- temp.75[, dist := sqrt(sum.dist)] # take sqrt

# combine back with categories
# all axes
dist.dr <- temp[data, c("sample.type.year",
                        "year", "Season") := list(i.sample.type.year,
                                                  i.year, i.Season), on = .(dr_match_name)]
# 75% variance
dist.75 <- temp.75[data, c("sample.type.year",
                           "year", "Season") := list(i.sample.type.year,
                                                     i.year, i.Season), on = .(dr_match_name)]

# Calculate delta of the two metrics
diff.df <- dcast(dist.75, dr_match_name + sample.type.year + Season ~ Metric, value.var = "dist")
diff.df <- diff.df[, delta := Bray - Sorensen]

diff.df <- melt(diff.df, id.vars = c("dr_match_name","sample.type.year","Season"),
                measure.vars = c("delta","Sorensen","Bray"),
                variable.name = "Metric", value.name = "dist")
sum.delta <- diff.df[, .(mean =  mean(dist, na.rm = T),
                        conf.int = conf.int(dist),
                        stdev = sd(dist, na.rm = T),
                        n = .N),
                    by = .(sample.type.year, Season, Metric)]

# check summary values for results
diff.df[,.(mean =  mean(dist, na.rm = T),
           stdev = sd(dist, na.rm = T),
           n = .N),
        by = .(Season, Metric)] %>% arrange(Metric, mean)


# calculate confidence interval and means of sample type and season combinations
# All axes
sum.dist <- dist.dr[, .(mean =  mean(dist, na.rm = T),
                       conf.int = conf.int(dist),
                       stdev = sd(dist, na.rm = T),
                       n = .N),
                   by = .(Metric, sample.type.year, Season)]
dist.dr[, .(mean =  mean(dist, na.rm = T),
                        stdev = sd(dist, na.rm = T),
                        n = .N),
                    by = .(Metric,Season)]

# only 75%
sum.dist75 <- dist.75[, .(mean =  mean(dist, na.rm = T),
                          conf.int = conf.int(dist),
                          stdev = sd(dist, na.rm = T),
                          n = .N),
                      by = .(Metric, sample.type.year, Season)]

dist.75[, .(mean =  mean(dist, na.rm = T),
            stdev = sd(dist, na.rm = T),
            n = .N),
        by = .(Metric,Season)]

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
                by.x = "ID", by.y = "dr_match_name")

#all.ax <- merge(dissim.dr[dissim.dr$Metric == "Sorensen",], 
#                dist.dr[dist.dr$Metric == "Sorensen",], 
#                by.x = "ID", by.y = "dr_match_name") # works for both metrics

(p <- ggplot(all.ax, aes(x = dist.x, y = dist.y, fill = sample.type.year.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  labs(x = expression(paste("Pairwise ", italic("D")["BC"])),
       y = expression(paste(italic("m")["BC"]^"n = 100%")))
)

ggsave("./Figures/General/allaxesdist_dissim_cor.png", p)

#------------------------------------------------------------------------------------------
# Extract dissimilarity


##########################
# Lollipop plots of distance
temp <- sum.ls[[3]][!is.nan(sum.ls[[3]]$mean),]
# Rename level for plotting
levels(temp$sample.type.year)[levels(temp$sample.type.year) == "Riverine \nLakes"] <- "Riv. Lakes"
names(colvec)[7] <- "Riv. Lakes"

plot.df <- dcast(temp, sample.type.year + Season + panels ~ Metric, value.var = "mean")

temp$Metric <- factor(temp$Metric, levels = c('delta','Bray','Sorensen'), labels = c("Delta","Abundance",
                                                                                  "Incidence"))


# Make a fake plot with all the legend units, habitat type, season and metric
legs <- get_legend(ggplot(temp[!temp$Metric == "Delta",], aes(x = sample.type.year)) +
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
  geom_linerange(aes(ymin = Sorensen, ymax = Bray,
                     colour = sample.type.year, group = Season),
                 colour = "gray30", linetype = "solid", position = position_dodge(0.7)) +
  geom_jitter(aes(y = Sorensen, colour = sample.type.year,
                  shape = Season), colour = "black",
              fill = "white", position = position_dodge(0.7), size = 3) +
  scale_fill_manual(values = colvec, name = "Habitat type") +
  scale_shape_manual(values = c(21, 23, 25)) +
  geom_jitter(aes(y = Bray, fill = sample.type.year,
                  shape = Season), position = position_dodge(0.7), size = 3) +
  scale_colour_manual(values = colvec) +
  labs(y= "Mean DNA-RNA distance", x = "Habitat Type") +
  lims(y = c(min(plot.df$Sorensen, na.rm = T),
             max(plot.df$Bray, na.rm = T))) +
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
  geom_linerange(aes(ymin = Sorensen, ymax = Bray,
                     colour = sample.type.year, group = Season),
                 colour = "gray30", linetype = "solid", position = position_dodge(0.7)) +
  geom_jitter(aes(y = Sorensen, colour = sample.type.year,
                  shape = Season), colour = "gray20",
              fill = "white", position = position_dodge(0.7), size = 3) +
  scale_fill_manual(values = colvec) +
  scale_shape_manual(values = c(21, 23, 25)) +
  geom_jitter(aes(y = Bray, fill = sample.type.year,
                  shape = Season), position = position_dodge(0.7), size = 3) +
  scale_colour_manual(values = colvec) +
  labs(y= "DNA-RNA distance", x = "Habitat Type") +
  lims(y = c(min(plot.df$Sorensen, na.rm = T),
             max(plot.df$Bray, na.rm = T))) +
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
          heights = c(2,1), align = "v", labels = "auto", vjust = .9))

ggsave("./Figures/Submission/Figure4.tiff", dist.p, 
       width = 18, height = 13, units = "cm", dpi = 300)

ggsave("./Figures/Final/distance_lollipop_n.png", dist.p, 
       width = 18, height = 13, units = "cm", dpi = 300)
ggsave("./Figures/Final/distance_lollipop_n.tiff", dist.p, 
       width = 18, height = 13, units = "cm", dpi = 300)


# What is behind these patterns? ------------------------------------------------------------------------
# Figure 5.: Where were OTUs first observed and what is their contribution to DNA and RNA pools? --------

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
       c("dr_match_name", "dna_type","sample.type.year","Season") := list(i.dr_match_name, i.dna_type,
                                                                    i.sample.type.year, i.Season),
       on = .(Sample)]

# cast so that we can calculate difference between DNA and RNA of each OTU
temp <- dcast(commat, dr_match_name + sample.type.year + OTU ~ dna_type, value.var = c("reads"))
temp[, diff := DNA - RNA]
temp <- temp[!is.na(diff),]

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
length(c(soil, soilwater, stream, upriver, reservoir, downriver, estuary)) == length(all.otus) # need to be TRUE

# 12795 of original 16368
# some OTUs are unique to ecosystems outside of the continuum

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
any(!(first.df$OTU %in% c(soil,soilwater,stream, upriver,reservoir,downriver,estuary)) == T) # need to be FALSE

# sanity check
first.df[is.na(DNA) & RNA >0,]
first.df[DNA == 0 & RNA >0,]

# Calculate ratio as a comarison
first.df[, ratio := DNA / RNA]

# merge with some meta data
first.df[met.df, c("sample.type.year","Season") := list(i.sample.type.year,
                                                       i.Season), on = .(dr_match_name)]

# remove any NAs in the dataset
first.df[is.na(DNA),]
first.df[is.na(RNA),]


# Who within the rank abundance curve is most reactive? --------------------------------------------
# Define abundance groups---------------------------------------------------------------------------

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
rel.df <- select_newest("./Objects", "201520162017_css_otu99_paper1_")
rel.df <- readRDS(
  paste0("./Objects/", rel.df))

setDT(rel.df)
# calculate means by sample type
means <- rel.df[, .(mean.css = mean(css.reads, na.rm = T),
                    sd.css = sd(css.reads, na.rm = T)), by = .(sample.type.year, dna_type, OTU)]
# order the abundances to make ranks
means <- means[mean.css != 0,] # remove all 0 observations
means <- means[order(mean.css, decreasing = T)]
means[, rank.abun := 1:.N, by = .(dna_type, sample.type.year)]
means[, log.mean := log1p(mean.css), by = .(dna_type, sample.type.year)]

# smooth and get derivative
classif.thres <- ddply(means, .(dna_type, sample.type.year), function(x){
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
                  sd.min = sd(min)), by = .(dna_type, ab.group)]
#saveRDS(classif.thres, "./Objects/abundance.classification.threshold.rds")


# Apply abundance groups -----------------------------------------------------------------------

# After consulting the means and deviations of the maximum and minimum thresholds across samples
# (see previous script) we settle with:
# Abundant >= 72 css reads
# Medium < 72 & >= 10 css reads
# Rare < 10 css reads

# Add abundance group to each observation
first.df[DNA >= 72, abg.dna := "Abundant"]
first.df[DNA < 72 & DNA >= 10, abg.dna := "Moderate"]
first.df[DNA < 10, abg.dna := "Rare"]
#first.df[DNA >= 66, abg.rna := "Abundant"]
#first.df[DNA < 66 & DNA >= 9, abg.rna := "Moderate"]
#first.df[DNA < 9, abg.rna := "Rare"]

# Extract reactive fraction: RNA > 0
rna.df <- first.df[RNA > 0, 
                   .(sum.rna.byotu = sum(RNA, na.rm = T),
                     mean.ratio = mean(ratio, na.rm = T)),
                   by = .(OTU, sample.type.year, abg.dna, first.obs, Season)]

# Extract all (DNA > 0 and RNA >0 + DNA > 0 and RNA = 0)
dna.df <- first.df[DNA > 0, 
                   .(sum.dna.byotu = sum(DNA, na.rm = T)),
                   by = .(OTU, sample.type.year, abg.dna, first.obs, Season)]
mean.df <- rna.df[dna.df, , on = .(OTU, sample.type.year, abg.dna, first.obs, Season)]

# Calculate the number of reads only in reactive fraction
ac.all <- first.df[RNA > 0,
                   .(otu.n.all.rna = .N,
                     sum.reads.all.rna = sum(RNA, na.rm = T)), by = .(sample.type.year, Season)]  # dr_match_name
mean.df[ac.all, "sum.reads.all.rna" := i.sum.reads.all.rna, on = .(sample.type.year,Season)]

# Sum of all reads of all observations of DNA (not only DNA >0 and RNA = 0, like before)
all <- first.df[DNA > 0,
                .(otu.n.all.dna = .N,
                  sum.reads.all.dna = sum(DNA, na.rm = T)), by = .(sample.type.year, Season)] #dr_match_name
mean.df[all, "sum.reads.all.dna" := i.sum.reads.all.dna, on = .(sample.type.year,Season)]

# Calculate percentage
mean.df[, perc.sum.rna.byotu := sum.rna.byotu * 100 / sum.reads.all.rna]
mean.df[, perc.sum.dna.byotu := sum.dna.byotu * 100 / sum.reads.all.dna]

# sanity check
mean.df[, .(sum.all = sum(perc.sum.rna.byotu, na.rm = T)),
        by = .(sample.type.year, Season)]
mean.df[, .(sum.all = sum(perc.sum.dna.byotu, na.rm = T)),
        by = .(sample.type.year, Season)] # all 100%

# make quartile bins to represent variation better
# overall mean simplifies pattern too much. Tried 4 bins first by:
# 1. < 25th quartile, 2. 25th quartile - Median, 3. Median - 75th quartile, 4. > 75th quartile
# 1+2 and 3+4 were relatively similar, and hence, decide on only two bins.
mean.df[perc.sum.rna.byotu < exp(-5), quan.bin := 1] #  < 50%
mean.df[perc.sum.rna.byotu >= exp(-5), quan.bin := 2] # >= 50%
# mean.df[perc.sum.rna.byotu < summary(perc.sum.rna.byotu)[3], quan.bin := 1] #  < 50%
# mean.df[perc.sum.rna.byotu >= summary(perc.sum.rna.byotu)[3], quan.bin := 2] # >= 50%

# Any OTUs not assigned a quartile bin?
any(is.na(mean.df$quan.bin))
# only OTUs that don't have RNA, that's fine

# save by OTU data frame for later
byotu.df <- mean.df

# summarise by bins
mean.df <- mean.df[,.(perc.sum.rna.otu = mean(perc.sum.rna.byotu, na.rm = T),
                      perc.sum.dna.otu = mean(perc.sum.dna.byotu, na.rm = T),
                      mean.ratio = mean(mean.ratio, na.rm = T)),
                   by = .(sample.type.year, abg.dna, first.obs, Season, quan.bin)]

mean.mean <- mean.df[quan.bin == 2, .(mean = mean(perc.sum.rna.otu),
                                      sd = sd(perc.sum.rna.otu)), by = .(first.obs)]

lins <- dlply(mean.df[quan.bin == 2,], .(sample.type.year, Season), function(x){
  lm(log(perc.sum.dna.otu) ~ log(perc.sum.rna.otu), data = x)
})


ggplot() +
  geom_point(aes(x = lins[[2]]$model$`log(perc.sum.rna.otu)`, y = lins[[2]]$model$`log(perc.sum.dna.otu)`)) +
  geom_line(aes(x = lins[[2]]$model$`log(perc.sum.rna.otu)`, y = lins[[2]]$fitted.values))

(slopes.df <- ldply(lins, coef))
mean(slopes.df$`log(perc.sum.rna.otu)`)

mean.mean %>% arrange(mean)

# Prepare for plotting
# overwrite soil colour, original beige is too hard to see
new.col <- colvec
new.col[names(new.col) == "Soil"] <- "gold"
# make df for annotation of Median
ann_df <- data.frame(perc.sum.dna.otu = -8,
                     perc.sum.rna.otu = as.numeric(summary(log(mean.df$perc.sum.rna.otu))[3]),
                     Season = factor("Summer",levels = c("Spring","Summer","Autumn")),
                     sample.type.year = factor("Soil"))
# leave panels empty where we don't have any data
line_df <- data.frame(Season = c("Summer","Spring", "Summer", "Spring", "Summer",
                                 "Spring", "Summer", "Autumn","Spring", "Summer", "Autumn",
                                 "Spring", "Summer", "Autumn",
                                 "Summer"),
                      sample.type.year = c("Soil","Soilwater","Soilwater", "Stream","Stream",
                                           "Upriver","Upriver","Upriver",
                                           "Downriver","Downriver","Downriver",
                                           "Reservoirs","Reservoirs","Reservoirs",
                                           "Estuary"))
#line_df$x <- 0 ; line_df$y <- summary(log(mean.df$perc.sum.rna.otu))[3]
line_df$x <- 0 ; line_df$y <- -5
line_df$Season <- factor(line_df$Season, levels = c("Spring","Summer","Autumn"))
setDT(line_df)
line_df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                  "Upriver", "Reservoirs", "Downriver",
                                                                  "Estuary"))]
# Make factor
mean.df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                  "Upriver", "Reservoirs", "Downriver",
                                                                  "Estuary"))]

# Categorize median bins
mean.df[, quan.bin := factor(quan.bin, levels = c(2,1),
                           labels = c("> 0.0067","< 0.0067"))]


(perc.cont.dna.rna <- ggplot(mean.df[!is.na(quan.bin),]) +
    theme_cust("pubr")+
    facet_grid(sample.type.year~Season) +
    geom_hline(data = line_df,
               aes(yintercept = y),
               linetype = "solid", colour = "gray70") +
    geom_smooth(aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), group = as.character(quan.bin)),
                method = "lm", se = F, colour = "gray30", size = 0.5, linetype = "solid",
                data = mean.df[!is.na(quan.bin) & quan.bin == "> 0.0067",]) +
    #ggforce::geom_mark_ellipse(data = subset(mean.df[!is.na(quan.bin) & Season == "Summer" 
    #                                        & sample.type.year == "Soil",], (quan.bin == 2)),
    #                  aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), group = as.character(quan.bin)),
    #                  linetype = "dashed", colour = "gray30", expand = unit(1.5, "mm")) +
    geom_point(aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), colour = first.obs,
                   shape = abg.dna, fill = quan.bin), alpha = 0.9, size =2) +
    geom_text(data = ann_df,
              aes(y = perc.sum.rna.otu +1, x = perc.sum.dna.otu+1.2), label = "Threshold", size =2.5,
              colour = "gray50") +
    scale_colour_manual(values = new.col,
                        name = "First \ndetected in") +
    scale_shape_manual(values = c(21,23,25),
                       name = "Local DNA\nabundance") +
    scale_fill_manual(values = c("white","gray80"),
                      name = "Bins of \n% RNA reads") +
    scale_y_continuous(limits = c(min(log(mean.df$perc.sum.rna.otu), na.rm = T) - 1, 
                                  max(log(mean.df$perc.sum.rna.otu), na.rm = T) + 1)) +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          strip.background = element_rect(fill = "gray20"),
          strip.text = element_text(colour = "white", size = 9)) +
    labs(x = "Mean % OTU contribution to\nlocal DNA pool (log-scale)",
         y = "Mean % OTU contribution to\nlocal RNA pool (log-scale)") +
    guides(shape = guide_legend(order = 2, override.aes=list(size = 3)),
           fill = guide_legend(order = 1, override.aes=list(shape=21, size = 3)),
           colour = guide_legend(order = 3, override.aes=list(size = 3))))

#
#stat_ellipse(data = subset(mean.df[!is.na(quan.bin),], (quan.bin == 2)),
#             aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), group = as.character(quan.bin)),
#             linetype = "dashed", colour = "grey50") +

# removed observations are OTUs with no RNA
# Annotate second x-axis (panels)
(ann.perc.cont.dna.rna <- annotate_figure(perc.cont.dna.rna + theme(legend.position = "none"),
                                      right = text_grob("Habitat Type", rot = 270, hjust = 0.7)))
# merge with earlier percentage first detected bar plots
#(first.obs.dna.rna <- ggarrange(frac.first.obs, perc.cont.dna.rna, nrow = 1, ncol = 2, labels = c("","c"), hjust = -10, widths = c(1,2)))
# Save
#ggsave("./Figures/Final/first.observed.contribution.dna.rna.png", 
#       first.obs.dna.rna, width = 25, height = 15, units = "cm", dpi = 300)
#ggsave("./Figures/Final/first.observed.contribution.dna.rna.tiff", 
#       first.obs.dna.rna, width = 25, height = 15, units = "cm", dpi = 300)

# Calculate the percentage they occupy in DNA ---------------------------------------------------------------

# Calculate perc of first observed of < Median, > Median, and no RNA OTUs
subdf <- byotu.df[, .(sum.first.medlow = sum(sum.dna.byotu, na.rm = T)), .(first.obs, Season, sample.type.year, quan.bin)]
all.sum <- byotu.df[,.(sum.dna = sum(sum.dna.byotu, na.rm = T)), .(Season, sample.type.year, quan.bin)]

subdf <- subdf[all.sum, , on = .(Season, sample.type.year, quan.bin)]
subdf[, perc.dna := sum.first.medlow * 100 / sum.dna]

# Categorize median bins
subdf[, quan.bin := factor(quan.bin, levels = c(2,1),
                              labels = c("> Threshold\n(reactive)","< Threshold\n(unreactive)"))]
# Add a new bin for the non RNA OTUs
subdf[is.na(quan.bin), quan.bin := "RNA = 0\n(unreactive)"]

# Level sample types
subdf[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                "Upriver", "Reservoirs", "Downriver",
                                                                "Estuary"))]
(un.perc.dna <- ggplot(subdf[quan.bin == "RNA = 0\n(unreactive)",], aes(x = sample.type.year)) +
    theme_cust("pubr") +
    geom_col(aes(y = perc.dna, fill = first.obs), colour = "gray20", size = 0.3) +
    facet_grid(quan.bin~Season, scale = "free_x", space = "free") +
    theme(axis.text.x = element_blank(), #axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = "gray20"),
          strip.text = element_text(colour = "white", size = 9)) +
    labs(x = "Habitat type", y = "% of reads in\nlocal DNA pool") +
    scale_fill_manual(values = colvec[names(colvec) %in% subdf$sample.type.year],
                      name = "First observed in") +
    guides(fill = "none"))

(perc.dna <- ggplot(subdf[!quan.bin == "RNA = 0\n(unreactive)",], aes(x = sample.type.year)) +
  theme_cust("pubr") +
  geom_col(aes(y = perc.dna, fill = first.obs), colour = "gray20", size = 0.3) +
  facet_grid(quan.bin~Season, scale = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        #axis.title.y = element_blank(),
        strip.background.x = element_blank(),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 9)) +
  labs(x = "Habitat type", y = "% of reads in\nlocal DNA pool") +
  scale_fill_manual(values = colvec[names(colvec) %in% subdf$sample.type.year],
                    name = "First observed in") +
  guides(fill = "none"))



(perc.cont.col <- 
  ggarrange(ggarrange(un.perc.dna, perc.dna, nrow = 2, heights = c(1,2.1), labels = c("a","c")),
          ann.perc.cont.dna.rna, ncol = 2, labels = c("","b"),
          widths = c(1.4,2),
          legend.grob = get_legend(perc.cont.dna.rna), legend = "right"))

#(perc.cont.col <- ggarrange(ggarrange(ann.perc.cont.dna.rna,
#                    legend.grob = get_legend(perc.cont.dna.rna), legend = "right"), 
#          perc.dna, ncol = 2, labels = "auto", widths = c(2,1.1)))
#(perc.cont.col <- ggarrange(ann.perc.cont.dna.rna, 
#                            perc.dna, ncol = 2, labels = "auto", widths = c(2,1.1),
#                            legend.grob = get_legend(perc.cont.dna.rna), legend = "bottom"))

ggsave("./Figures/Submission/Figure5.tiff", 
       perc.cont.col, width = 25, height = 15, units = "cm", dpi = 300)

ggsave("./Figures/Final/first.observed.contribution.thresbins_n.png", 
       perc.cont.col, width = 25, height = 15, units = "cm", dpi = 300)
ggsave("./Figures/Final/first.observed.contribution.thresbins_n.tiff", 
       perc.cont.col, width = 25, height = 15, units = "cm", dpi = 300)

subdf[, .(sum = sum(perc.dna)), by = .(quan.bin)]

# Calculate % OTU for supplementary -------------------------------------------------------------------------

# What's the proportion of reactive vs unreactive? -----------------------------------------------------------------------------------

byotu.df[, un.reactive := "unreactive"]
byotu.df[quan.bin == 2, un.reactive := "reactive"]

first.df[byotu.df, "un.reactive" := i.un.reactive, on = .(OTU, sample.type.year, Season)]

# remove NAs
subdf <- first.df[!is.na(un.reactive),]

# count the number of OTUs
all <- subdf[DNA > 0,
                .(otu.n.all.dna = .N,
                  sum.reads.all.dna = sum(DNA, na.rm = T)), by = .(sample.type.year, Season)] #dr_match_name

# count number of OTUs in each group
group.count <- subdf[DNA > 0,.(otu.n.un.reactive = .N), by = .(sample.type.year, Season, un.reactive)]

# calculate percentage
group.count[all, "otu.n.all.dna" := i.otu.n.all.dna, on = .(sample.type.year, Season)]
group.count[, perc.otu := otu.n.un.reactive * 100 / otu.n.all.dna]

# sanity check
group.count[, .(sum = sum(perc.otu)), by = .(sample.type.year, Season)]

(unreac.pie <- ggplot(group.count, aes(y = perc.otu, fill = un.reactive)) +
  theme_void() +
  facet_grid(sample.type.year~Season) +
  geom_col(aes(x = ""), colour = "gray20", size = 0.3) +
  coord_polar("y", start = 0) +
  scale_fill_viridis_d(direction = -1, option = "magma", name = "Fraction", labels = c("Reactive","Unreactive"))+
  theme(legend.position = "top")) 

ggsave("./Figures/Final/reactive.unreative.pies_n.png", 
       unreac.pie, width = 10, height = 15, units = "cm", dpi = 300)


(unreac.bar <- ggplot(group.count, aes(y = perc.otu, fill = un.reactive)) +
  theme_cust("pubr") +
  facet_grid(.~Season, scale = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 9)) +
  geom_col(aes(x = sample.type.year),colour = "gray20", size = 0.3) +
  scale_fill_viridis_d(direction = -1, option = "cividis", name = "Fraction", labels = c("Reactive","Unreactive")) +
  labs(y = "% of bacterial OTUs", x = "Habitat type"))

ggsave("./Figures/Final/reactive.unreative.barplot_n.png", 
       unreac.bar, width = 12, height = 6, units = "cm", dpi = 300)

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
