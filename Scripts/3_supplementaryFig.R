#-- Script for the publication:
#-- Title
#-- Responsible author: Masumi Stadler

# This script is the last of a series of scripts that were used to analyse the data
# used in the publication.

# In this script we create additional supplementary figures
# that are not included directly in R markdown for the Supplements or the previous scripts,
# as data is too big to be directly knitted within R markdown

###----------------###
#-   Supplements   - #
###----------------###

# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("phyloseq", # wrangling
              "tidyverse", "data.table", # wrangling
              "ggpubr", # plotting
              "vegan", "ade4") # stats


### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)


### Functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")

# Set seed for session and reproducibility of permutations
# (just for consistency, no random iteration in this script)
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
            "#050416FF") #Estuary)

names(colvec) <- as.character(sample.factors) # assign names for later easy selection

# set theme for plotting
theme_set(theme_bw())

# Fig. S2 -------------------------------------------------------------------------------------

alpha.df <- select_newest("./Output/", "alpha_div_otu99_summary_")
alpha.df <- read.csv(paste0("./Output/", alpha.df), sep = ",", stringsAsFactors = F)

# make to data table
setDT(alpha.df)

alpha.df$Data <- factor(alpha.df$Data, levels = c("css", "lib15000", "lib25000", "lib50000"),
                          labels = c("CSS", "Rarefied: Lib15000", "Rarefied: Lib25000","Rarefied: Lib50000"))
alpha.df$DnaType <- factor(alpha.df$DnaType, levels = c("DNA", "cDNA"),
                             labels =  c("DNA", "RNA"))

melt.alpha <- melt(alpha.df, id.vars = c("DR.names", "Data", "DnaType", "sample.type.year","Season"),
                   measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
                   variable.name = "Diversity",
                   value.name = "Index")


# make plot with ggpubr to include pearson's correlation outputs directly in the plot
(rar.box <- ggplot(melt.alpha, aes(x = DnaType, y = Index)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(Diversity ~ Data, scales = "free") +
  labs(x = "Nucleic acid type", y = expression(paste(alpha," - Diversity"))))

ggsave("./Figures/Final/alphadiv_rar_comp_boxplots.png", rar.box,
       width = 15, height = 12, unit = "cm")

cast.alpha <- melt.alpha %>%
  group_by(DR.names, Data, Diversity, DnaType) %>%
  summarise(Season = unique(Season),
            Index = mean(Index, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = DnaType, values_from = Index)

cast.alpha$Season <- factor(cast.alpha$Season, levels = c("spring", "summer", "autumn"),
                           labels =  c("Spring", "Summer", "Autumn"))

(rar.lin <- ggplot(cast.alpha,
                  aes(x = RNA, y = DNA, fill = Season)) +
  geom_point(shape = 21) +
  facet_wrap(Diversity~Data, scales = "free") +
  scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00"))
) # colour-blind friendly

ggsave("./Figures/Final/alphadiv_rar_comp_scatter.png", rar.lin,
       width = 18, height = 15, unit = "cm")


# Fig. S6 -------------------------------------------------------------------------------------
# Abundance classification
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
(p <- annotate_figure(p, bottom = text_grob("OTU Rank")))
ggsave("./Figures/Final/abundance_class_ex_otu99.png", p,
       width = 22, height = 9, unit = "cm")


# Fig. S7 -------------------------------------------------------------------------------------
# Taxonomic composition
# We want to show the taxonomic composition of our samples plus the abundance
# As we have too many samples, best would probably be to calcualte the mean abundance for each group
# Groups are: sample.type.year + DnaType + Season

# melt OTU table
pb.df <- as.data.frame(otu_mat(pb)) # make data.frame
setDT(pb.df, keep.rownames = "Sample") # make data.table

# melt
pb.df <- melt.data.table(pb.df, 
                         id.vars = "Sample", # skip measure.var, takes all columns
                         variable.name = "OTU",
                         value.name = "css")

# add taxonomy data to the abundance data
tax.df <- as.data.frame(tax.tab) # make data.frame
setDT(tax.df, keep.rownames = "OTU") # make data.table
pb.df <- pb.df[tax.df, on = .(OTU)] # merge

# add meta data
meta <- sample_df(met.df) %>%
  dplyr::select(DnaType, Year, Season, sample.type.year, soilorwater, LibrarySize)
setDT(meta, keep.rownames = "Sample") # make data.table

pb.df <- pb.df[meta, on = .(Sample)] # merge

# calculate mean abundance (css) for each category
mean.pb <- pb.df[, .(css.mean = mean(css, na.rm = T),
                     css.sd = sd(css, na.rm = T)),
                 by = .(OTU, sample.type.year, Season, DnaType)][
                   , css.sum := sum(css.mean, na.rm = T),
                   by = .(sample.type.year, Season, DnaType)
                 ][, css.rel := css.mean * 1 / css.sum]
mean.pb <- mean.pb[tax.df, on = .(OTU)] # merge with taxonomy data
mean.pb[, ID := paste(DnaType, Season, sample.type.year, sep = "_")] # add plot ID

# mean of means for results section
temp <- mean.pb[, .(sum = sum(css.rel, na.rm = T)), .(sample.type.year,phylum, DnaType, Season)]
temp <- temp[, .(mean = mean(sum) *100, sd = sd(sum)), .(phylum)]
temp[order(mean, decreasing = T),]

# calculate mean library size per category
lib.size <- pb.df[, .(lib.mean = mean(LibrarySize, na.rm = T),
                      lib.sd = sd(LibrarySize, na.rm = T)),
                  by = .(sample.type.year, DnaType, Season)]
lib.size[, ID := paste(DnaType, Season, sample.type.year,sep = "_")] # add plot ID

# Overwrite IDs as factors to set a order for plotting
mean.pb$ID <- factor(mean.pb$ID, 
                     levels = c("DNA_Spring_Soil", "DNA_Spring_Sediment",
                                "DNA_Spring_Soilwater",
                                "DNA_Spring_Groundwater", "DNA_Spring_Stream",
                                "DNA_Spring_Tributary", "DNA_Spring_Riverine \nLakes",
                                "DNA_Spring_Headwater \nPonds", "DNA_Spring_Lake",
                                "DNA_Spring_Upriver", "DNA_Spring_Downriver",
                                "DNA_Spring_Reservoirs", 
                                "DNA_Spring_Estuary",
                                "DNA_Summer_Soil", "DNA_Summer_Sediment",
                                "DNA_Summer_Soilwater",
                                "DNA_Summer_Stream",
                                "DNA_Summer_Tributary", "DNA_Summer_Riverine \nLakes",
                                "DNA_Summer_Headwater \nPonds", "DNA_Summer_Lake",
                                "DNA_Summer_Upriver", "DNA_Summer_Downriver",
                                "DNA_Summer_Reservoirs", 
                                "DNA_Summer_Estuary",
                                "DNA_Autumn_Tributary", "DNA_Autumn_Riverine \nLakes",
                                "DNA_Autumn_Lake", "DNA_Autumn_Upriver", "DNA_Autumn_Downriver",
                                "DNA_Autumn_Reservoirs",
                                "RNA_Spring_Soil",
                                "RNA_Spring_Soilwater",
                                "RNA_Spring_Stream",
                                "RNA_Spring_Tributary", "RNA_Spring_Riverine \nLakes",
                                "RNA_Spring_Lake",
                                "RNA_Spring_Upriver", "RNA_Spring_Downriver",
                                "RNA_Spring_Reservoirs",
                                "RNA_Summer_Soil", "RNA_Summer_Sediment",
                                "RNA_Summer_Soilwater",
                                "RNA_Summer_Stream",
                                "RNA_Summer_Tributary", "RNA_Summer_Riverine \nLakes",
                                "RNA_Summer_Lake",
                                "RNA_Summer_Upriver", "RNA_Summer_Downriver",
                                "RNA_Summer_Reservoirs", 
                                "RNA_Summer_Estuary",
                                "RNA_Autumn_Tributary", "RNA_Autumn_Riverine \nLakes",
                                "RNA_Autumn_Lake", "RNA_Autumn_Upriver", "RNA_Autumn_Downriver",
                                "RNA_Autumn_Reservoirs"))

lib.size$ID <- factor(lib.size$ID, 
                      levels = c("DNA_Spring_Soil", "DNA_Spring_Sediment",
                                 "DNA_Spring_Soilwater",
                                 "DNA_Spring_Groundwater", "DNA_Spring_Stream",
                                 "DNA_Spring_Tributary", "DNA_Spring_Riverine \nLakes",
                                 "DNA_Spring_Headwater \nPonds", "DNA_Spring_Lake",
                                 "DNA_Spring_Upriver", "DNA_Spring_Downriver",
                                 "DNA_Spring_Reservoirs", 
                                 "DNA_Spring_Estuary",
                                 "DNA_Summer_Soil", "DNA_Summer_Sediment",
                                 "DNA_Summer_Soilwater",
                                 "DNA_Summer_Stream",
                                 "DNA_Summer_Tributary", "DNA_Summer_Riverine \nLakes",
                                 "DNA_Summer_Headwater \nPonds", "DNA_Summer_Lake",
                                 "DNA_Summer_Upriver", "DNA_Summer_Downriver",
                                 "DNA_Summer_Reservoirs", 
                                 "DNA_Summer_Estuary",
                                 "DNA_Autumn_Tributary", "DNA_Autumn_Riverine \nLakes",
                                 "DNA_Autumn_Lake", "DNA_Autumn_Upriver", "DNA_Autumn_Downriver",
                                 "DNA_Autumn_Reservoirs",
                                 "RNA_Spring_Soil",
                                 "RNA_Spring_Soilwater",
                                 "RNA_Spring_Stream",
                                 "RNA_Spring_Tributary", "RNA_Spring_Riverine \nLakes",
                                 "RNA_Spring_Lake",
                                 "RNA_Spring_Upriver", "RNA_Spring_Downriver",
                                 "RNA_Spring_Reservoirs",
                                 "RNA_Summer_Soil", "RNA_Summer_Sediment",
                                 "RNA_Summer_Soilwater",
                                 "RNA_Summer_Stream",
                                 "RNA_Summer_Tributary", "RNA_Summer_Riverine \nLakes",
                                 "RNA_Summer_Lake",
                                 "RNA_Summer_Upriver", "RNA_Summer_Downriver",
                                 "RNA_Summer_Reservoirs", 
                                 "RNA_Summer_Estuary",
                                 "RNA_Autumn_Tributary", "RNA_Autumn_Riverine \nLakes",
                                 "RNA_Autumn_Lake", "RNA_Autumn_Upriver", "RNA_Autumn_Downriver",
                                 "RNA_Autumn_Reservoirs"))


# colour blind friendly, derived from https://medialab.github.io/iwanthue/
col_vector<-c("#350070","#dcd873","#06a644","#b9ce40","#003499","#6b8aff","#f5d249","#ec83f6","#7beb8b",
              "#c947b1","#01b072","#df3587","#006f26","#ff83d0","#215a00","#99a1ff","#668200",
              "#015db8","#e77e28","#019cf8","#b5221d","#67b8ff","#bd0a35","#b2e294","#840066",
              "#314800","#ffb0ed","#954700","#3d0e52","#ff9b61","#59003b","#ff6e83","#aa7dbf","#620009")

# normal
#col_vector<-c("#6893ff","#bbce1a","#6b63ed","#18cc58","#b22bb2",
#"#71dd6e","#ff6fea","#428700","#005eca","#f1bf45","#36508f","#d78600","#4bb8ff",
#"#ce5300","#5bd5f6","#da0053","#01823f","#c80074","#97d68c","#a5006f","#bccf63",
#"#85307a","#848a00","#ff94d5","#006845","#ff6f91","#00b5be","#a8041e","#9ac89a",
#"#ff7f64","#707c48","#ff9b3a","#635b00","#9c5700")

x.types.labs<-c("Soil","Sediment","Soilwater","Groundwater","Stream",
                "Tributary","Riverine Lakes","Headwater Ponds","Lake",
                "Upriver","Downriver", "Reservoirs","Estuary",
                "Soil","Sediment","Soilwater","Stream","Tributary","Riverine Lakes",
                "Upstream Ponds","Lake","Upriver","Downriver", "Reserviors", "Estuary",
                "Tributary","Riverine Lakes","Lake","Upriver","Downriver", 'Reservoirs',
                "Soil","Soilwater","Stream","Tributary","Riverine Lakes",
                "Lake","Upriver","Downriver", "Reservoirs",
                "Soil","Sediment","Soilwater",
                "Stream","Tributary","Riverine Lakes","Lake","Upriver",
                "Downriver", "Reservoirs", "Estuary",
                "Tributary","Riverine Lakes","Lake","Upriver",
                "Downriver", "Reservoirs")

(lib.bar <- 
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
  coord_cartesian(clip="off"))

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
(labelled.tax <- tax + 
  coord_cartesian(ylim=c(0,1), clip="off") +
  annotate("segment", x = 0, xend = 13.5, y = -0.41, yend = -0.41) +
  annotate("segment", x = 13.7, xend = 25.5, y = -0.41, yend = -0.41) +
  annotate("segment", x = 25.7, xend = 31.5, y = -0.41, yend = -0.41) +
  annotate("segment", x = 32, xend = 40.5, y = -0.41, yend = -0.41) +
  annotate("segment", x = 40.7, xend = 51.5, y = -0.41, yend = -0.41) +
  annotate("segment", x = 51.7, xend = 57.5, y = -0.41, yend = -0.41) +
  annotate("text", x = 6.75, y = -0.44, label = "Spring") +
  annotate("text", x = 19.6, y = -0.44, label = "Summer") +
  annotate("text", x = 28.6, y = -0.44, label = "Autumn") +
  annotate("text", x = 36.25, y = -0.44, label = "Spring") +
  annotate("text", x = 46.1, y = -0.44, label = "Summer") +
  annotate("text", x = 54.6, y = -0.44, label = "Autumn") +
  annotate("segment", x = 0, xend = 31.5, y = -0.48, yend = -0.48) +
  annotate("segment", x = 32, xend = 57.5, y = -0.48, yend = -0.48) +
  annotate("text", x = 15.75, y = -0.52, label = "DNA") +
  annotate("text", x = 44.75, y = -0.52, label = "RNA")
)

combo  <- ggarrange(lib.bar, labelled.tax, nrow = 2,
                    align = "v", heights = c(0.2, 0.8), legend.grob = tax.leg, legend = "right")

#tax.leg, nrow = 2, heights = c(0.9, 0.1), widths = c(0.8, 0.2)
ggsave("./Figures/Final/Tax_LibSiz_Phyla_otus.png", combo,
       width = 300, height = 203, unit = "mm", dpi = 400)


# Fig. S8 -------------------------------------------------------------------------------------
# Only terrestrial PCoA

ter <- subset_samples(pb, sample.type.year == "Soil" |
                        sample.type.year == "Sediment" |
                        sample.type.year == "Soilwater" |
                        sample.type.year == "Groundwater")

# remove ASVs that do not appear in this dataset
ter <- prune_taxa(taxa_sums(ter) != 0, ter)
ter <- prune_samples(sample_sums(ter) != 0, ter)

pb.mat <- otu_mat(ter)
pb.mat <- log2(pb.mat + 1)

nrow(pb.mat)
ncol(pb.mat)
# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # TRUE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)
# plot with custom function
dna.pcoa <- plot_pcoa(pb.bray.pcoa, 
                      physeq = ter, colours = colvec, output = T)

ggsave("./Figures/General/terr_PCoA.png", dna.pcoa$plot,
       width = 12, height = 10, unit = "cm")

# Fig. S9 -------------------------------------------------------------------------------------
# RNA PCoA

# subset only RNA samples
rna <- subset_samples(pb, DnaType == "RNA")

# remove ASVs that do not appear in this dataset
rna <- prune_taxa(taxa_sums(rna) != 0, rna)
rna <- prune_samples(sample_sums(rna) != 0, rna)

# extract ASV matrix
pb.mat <- otu_mat(rna)
pb.mat <- log2(pb.mat + 1)

nrow(pb.mat)
ncol(pb.mat)
# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray) # TRUE

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)
# plot with custom function
rna.pcoa <- plot_pcoa(pb.bray.pcoa, 
                      physeq = rna, colours = colvec, output = T)

p <- rna.pcoa$plot + guides(alpha = "none")

# save
ggsave(paste0("./Figures/Final/PCoA_log_RNA_SampleType.png"),  p,
       width = 12, height = 10, unit = "cm")
