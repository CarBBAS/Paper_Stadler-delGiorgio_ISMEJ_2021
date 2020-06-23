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
library(doMC) # parallel computing
library(vegan) #vegdist, diversity
library(ape) #pcoa
library(ade4) #is.euclid

#-----------#
# FUNCTIONS #
#-----------#
source("./Functions/big_data.R")
source("./Functions/phyloseq_index.R")
source("./Functions/peak_identification.R")


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

############
# Analysis #
############
# make ordinations
# tried NMDS stress does not reach convergence

####################
# Both DNA and RNA #
####################
# extract species table with species in columns
pb.mat <- t(otu_mat(pb))

# 1. PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray.pcoa <- ape::pcoa(pb.bray)

# 2. PCoA with Jaccard (presence-absence)
pb.jac <- vegdist(pb.mat, method = "jaccard", binary = T)
is.euclid(pb.jac) # TRUE
pb.jac.pcoa <- ape::pcoa(pb.jac)

#cl <-
# Mantel test to see if matrices are different
#mantel(pb.bray, pb.jax, method = "pearson", permutations = 999, parallel = )

# extract scores and variance explained
pb.scores <- rbind(data.frame(Sample = as.character(row.names(pb.bray.pcoa$vectors)),
                        Metric = "Bray", pb.bray.pcoa$vectors[,1:2], stringsAsFactors = F),
                   data.frame(Sample = as.character(row.names(pb.jac.pcoa$vectors)),
                              Metric = "Jaccard", pb.jac.pcoa$vectors[,1:2], stringsAsFactors = F))  # get first two axes
pb.var <- rbind(data.frame(x = round(100 * pb.bray.pcoa$values$Eigenvalues[1] / sum(pb.bray.pcoa$values$Eigenvalues), 2),
                           y = round(100 * pb.bray.pcoa$values$Eigenvalues[2] / sum(pb.bray.pcoa$values$Eigenvalues), 2),
                           Metric = "Bray", stringsAsFactors = F),
                data.frame(x = round(100 * pb.jac.pcoa$values$Eigenvalues[1] / sum(pb.jac.pcoa$values$Eigenvalues), 2),
                           y = round(100 * pb.jac.pcoa$values$Eigenvalues[2] / sum(pb.jac.pcoa$values$Eigenvalues), 2),
                           Metric = "Jaccard", stringsAsFactors = F))
pb.scores <- merge(pb.scores, pb.var, by = "Metric")

# merge with a selection of meta data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, DnaType, LibrarySize,
                                                   bact.abundance, bact.production,catchment.area, lat, long), 
                   stringsAsFactors = F)
pdataframe <- merge(pb.scores, meta, by = "Sample")
pdataframe$Sample <- as.character(pdataframe$Sample)

#View(data.frame(pb.scores$Sample, str_replace(pb.scores$Sample, "D$", ""),
#           str_replace(pb.scores$Sample, "R$", ""), pb.scores$DnaType))

# overwrite factor levels
pdataframe$sample.type.year <- factor(pdataframe$sample.type.year, levels = c("Soil","Sediment",
                                                                            "Soilwater","Hyporheicwater", 
                                                                            "Wellwater","Stream", "Tributary",
                                                                            "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                            "Upriver","RO3", "RO2", "RO1","Deep",
                                                                            "Downriver",
                                                                            "Marine"),
                                     labels = c("Soil","Sediment",
                                                "Soilwater","Hyporheicwater", 
                                                "Wellwater","Stream", "Tributary",
                                                "Headwater \nLake", "Upstream \nPond", "Lake", "Lake",
                                                "Upriver","RO3","RO2", "RO1","Hypolimnion","Downriver",
                                                "Estuary"))
pdataframe$Season <- factor(pdataframe$Season, levels = c("spring", "summer", "autumn"), 
                           labels = c("Spring", "Summer","Autumn"))
pdataframe$DnaType <- factor(pdataframe$DnaType, levels = c("DNA", "cDNA"), labels = c("DNA", "RNA"))

# create colour vector for plotting
colvec <- c("red4","chocolate3","orangered2","orange3",
            "khaki","cadetblue","darksalmon",
            "darkolivegreen","darkolivegreen3","gold",
            "royalblue","mediumorchid4", "violet","palevioletred2","navy","skyblue",
            "seagreen3")
#"chocolate4","seagreen3","azure","gold"
# separate legends for plotting
(
  bot.leg <-
    get_legend(
      ggplot(pdataframe %>% filter(Metric == "Bray"), aes(x = Axis.1, y = Axis.2)) +
        theme_bw() +
        geom_point(aes(shape = Season, alpha = DnaType),
                   size = 2.5) +
        theme(legend.position = "bottom") +
        scale_shape_manual(values = c(21, 23, 25)) +
        scale_alpha_manual(values = c(1, 0.5), name = "Nucleic Acid Type") +
        guides(shape = guide_legend(order = 1),
               alpha = guide_legend(order = 2))
    )
)

# main plot
(pcoa.bray <- ggplot(pdataframe %>% filter(Metric == "Bray"), aes(x = Axis.1, y = Axis.2)) +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
               size = 2.5) +
    #scale_fill_viridis_d(name = "Sample Type") +
    scale_fill_manual(values = colvec, name = "Sample Type") +
    scale_shape_manual(values = c(21,23,25)) +
    scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
    labs(x = paste("PCoA 1 [", pdataframe[pdataframe$Metric == "Bray",]$x,"%]"), 
         y = paste("PCoA 2 [", pdataframe[pdataframe$Metric == "Bray",]$y,"%]")) +
    theme(legend.key.size = unit(1.5, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(override.aes=list(shape=21)),
           alpha = "none", shape = "none"))

(pcoa.jac <- ggplot(pdataframe %>% filter(Metric == "Jaccard"), aes(x = Axis.1, y = Axis.2)) +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
               size = 2.5) +
    #scale_fill_viridis_d(name = "Sample Type") +
    scale_fill_manual(values = colvec, name = "Sample Type") +
    scale_shape_manual(values = c(21,23,25)) +
    scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
    labs(x = paste("PCoA 1 [", pdataframe[pdataframe$Metric == "Jaccard",]$x,"%]"), 
         y = paste("PCoA 2 [", pdataframe[pdataframe$Metric == "Jaccard",]$y,"%]")) +
    theme(legend.key.size = unit(1.5, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(override.aes=list(shape=21)),
           alpha = "none", shape = "none"))

# arrange main plot and legend
(pcoa.arr <- ggarrange(ggarrange(pcoa.bray, pcoa.jac,
                       ncol = 2,
                       common.legend = T,
                       labels = c("a", "b"),
                       legend = "right"
                      ), bot.leg,
                      nrow = 2, heights = c(3, 0.2)))

ggsave("./Figures/Final/All_DNARNA_SampleType_PCoA.tiff", pcoa.arr,
       width = 30, height = 15, unit = "cm")
#############

rm(pcoa.bray, pcoa.jac, pcoa.arr,bot.leg)
####################################################################################
#-----------------------------------------------#
# Calculate distance between DNA and RNA points #
#-----------------------------------------------#

#View(data.frame(pb.scores$Sample, pb.scores$sample.type.year, pb.scores$Year))

# correct a few wrong sample names for matching DNA and RNA
pdataframe[pdataframe$Sample == "RO2R52R", "Sample"] <- "RO2.52R"
pdataframe[pdataframe$Sample == "SWR34R", "Sample"] <- "SW34R"
pdataframe[pdataframe$Sample == "RO2.36pD", "Sample"] <- "RO2.36D"
pdataframe[pdataframe$Sample == "RO2.36pR", "Sample"] <- "RO2.36R"
pdataframe[pdataframe$Sample == "RO2111.60mD", "Sample"] <- "RO2111.90mD"
pdataframe[pdataframe$Sample == "RO2.30DPR", "Sample"] <- "RO2.30R" # two DNA
pdataframe[pdataframe$Sample == "RO301.HypoR", "Sample"] <- "RO31.HypoR"
pdataframe[pdataframe$Sample == "RO301R", "Sample"] <- "RO31R" 
pdataframe[pdataframe$Sample == "RO304R", "Sample"] <- "RO34R" 
pdataframe[pdataframe$Sample == "RO307R", "Sample"] <- "RO37R" 

pdataframe[pdataframe$Sample == "L230R", "Sample"] <- "L330R" # L230 does not exist
#pdataframe[pdataframe$Sample == "TR49", "Sample"] <- # two DNA

# remove Ds and Rs to match counterpart samples
pdataframe$ID[pdataframe$DnaType == "DNA"] <- str_replace(pdataframe$Sample[pdataframe$DnaType == "DNA"], "D$", "")
pdataframe$ID[pdataframe$DnaType == "RNA"] <- str_replace(pdataframe$Sample[pdataframe$DnaType == "RNA"], "R$", "")

# export table to look at point positions in GIS
#write.table(pb.scores, "./Output/BrayCurtis_scores_withmeta.csv", sep = ",", dec = ".", row.names = F)

# calculate mean coordinates for duplicates
sum <- pdataframe %>% 
  filter(!Year == 2015) %>% 
  dplyr::group_by(ID, DnaType, Metric) %>%
  dplyr::summarise(x = mean(Axis.1), y = mean(Axis.2),
                   n = n()) %>%
  ungroup()

setDT(sum)
# Calculate distance
temp <- dcast(sum, ID + Metric ~ DnaType, value.var = c("x","y"))
temp <- dcast(sum, ID + Metric ~ DnaType, value.var = c("x","y"))
temp[, distance := sqrt((x_DNA - x_RNA)^2 + (y_DNA - y_RNA)^2)]

# combine back with categories
dist.dr <- temp[pdataframe, c("sample.type.year",
                  "Year", "Season" ) := list(i.sample.type.year,
                                             i.Year, i.Season), on = .(ID)]
# add new column to split plot into main and side panel
dist.dr[, panels := "main"]
dist.dr[sample.type.year == "Tributary" |
          sample.type.year == "Lake" |
          sample.type.year == "Headwater \nLake" |
          sample.type.year == "Sediment", panels := "side"]

write.table(dist.dr, "./Output/bray_pcoa_dnarna_distance.csv",
            sep = ";", dec = ".", row.names = F)

# Bray Curtis
# plot main plot
(
  main.b <-
    ggplot(dist.dr[panels == "main" & Metric == "Bray", ], aes(
      x = sample.type.year, y = distance, fill = Season
    )) +
    theme_pubr() +
    geom_boxplot(width = 0.5, outlier.size = 0.5, size = 0.3) + # outlier.alpha = 0
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    new_scale_fill() +
    #geom_point(aes(fill = as.character(Year)), position = position_jitterdodge(), colour = "black", shape = 21) + #alpha = 0.5,
    #scale_fill_manual(name = "Year", values = c("white","gray20")) +
    labs(x = "Sample type", y = "Distance between DNA and RNA \nin PCoA space (Bray-Curtis)") +
    #facet_grid(.~Year, scales = "free") +
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
    ggplot(dist.dr[panels == "side" & Metric == "Bray", ], aes(
      x = sample.type.year, y = distance, fill = Season
    )) +
    theme_pubr() +
    geom_boxplot(width = 0.5, outlier.size = 0.5, size = .3) + # outlier.alpha = 0
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    new_scale_fill() +
    #geom_point(
    #  aes(fill = as.character(Year)),
    #  position = position_jitterdodge(),
    #  colour = "black",
    #  shape = 21
    #) + #alpha = 0.5,
    #scale_fill_manual(name = "Year", values = c("white", "gray20")) +
    labs(x = "Sample type", y = "Distance between DNA and RNA \nin PCoA space") +
    lims(y = c(0, 0.6)) +
    #facet_grid(.~Year, scales = "free") +
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

# Jaccard
# plot main plot
(
  main.j <-
    ggplot(dist.dr[panels == "main" & Metric == "Jaccard", ], aes(
      x = sample.type.year, y = distance, fill = Season
    )) +
    theme_pubr() +
    geom_boxplot(width = 0.5, outlier.size = 0.5, size = 0.3) + # outlier.alpha = 0
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    new_scale_fill() +
    #geom_point(aes(fill = as.character(Year)), position = position_jitterdodge(), colour = "black", shape = 21) + #alpha = 0.5,
    #scale_fill_manual(name = "Year", values = c("white","gray20")) +
    labs(x = "Sample type", y = "Distance between DNA and RNA \nin PCoA space (Jaccard)") +
    #facet_grid(.~Year, scales = "free") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title = element_text(size = 10)
    )
)

# side panel
(
  side.j <-
    ggplot(dist.dr[panels == "side" & Metric == "Jaccard", ], aes(
      x = sample.type.year, y = distance, fill = Season
    )) +
    theme_pubr() +
    geom_boxplot(width = 0.5, outlier.size = 0.5, size = .3) + # outlier.alpha = 0
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    new_scale_fill() +
    #geom_point(
    #  aes(fill = as.character(Year)),
    #  position = position_jitterdodge(),
    #  colour = "black",
    #  shape = 21
    #) + #alpha = 0.5,
    #scale_fill_manual(name = "Year", values = c("white", "gray20")) +
    labs(x = "Sample type", y = "Distance between DNA and RNA \nin PCoA space") +
    lims(y = c(0, 0.6)) +
    #facet_grid(.~Year, scales = "free") +
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


# combine both plots
(p <- ggarrange(
  ggarrange(main.b, 
          side.b,
          widths = c(3,1),
          ncol = 2, nrow = 1, 
          common.legend = T,
          legend = "top",
          align = "h",
          font.label = list(size = 10)),
  ggarrange(main.j, 
            side.j,
            widths = c(3,1),
            ncol = 2, nrow = 1, 
            common.legend = F,
            legend = "none",
            align = "h",
            font.label = list(size = 10)),
  nrow = 2))
# add x axis title to be in the middle of two panels
(p <- annotate_figure(p, bottom = text_grob("Sample Type")))

# save
ggsave("./Figures/Final/All_log_DNARNA_withinPCoA_distance.tiff", 
       p, width = 18, height = 15, unit = "cm")

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
  # quantifies uncertainty associated with predicted the identity of a new taxa
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
  # the more variation in abundances between different taxa within the community the lower isJ
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

#------------------------------------------------------------------------------------------------#
# Combine with bray distance data
# first get some meta data
alpha.df <- merge(alpha.df , sample_df(pb) %>% dplyr::select(Sample = DadaNames, DnaType))

# correct a few wrong sample names for matching DNA and RNA
alpha.df[alpha.df$Sample == "RO2R52R", "Sample"] <- "RO2.52R"
alpha.df[alpha.df$Sample == "SWR34R", "Sample"] <- "SW34R"
alpha.df[alpha.df$Sample == "RO2.36pD", "Sample"] <- "RO2.36D"
alpha.df[alpha.df$Sample == "RO2.36pR", "Sample"] <- "RO2.36R"
alpha.df[alpha.df$Sample == "RO2111.60mD", "Sample"] <- "RO2111.90mD"
alpha.df[alpha.df$Sample == "RO2.30DR", "Sample"] <- "RO2.30R" # two DNA
#alpha.df[alpha.df$Sample == "TR49", "Sample"] <- # two DNA


# remove Ds and Rs to match counterpart samples
alpha.df$ID[alpha.df$DnaType == "DNA"] <- str_replace(alpha.df$Sample[alpha.df$DnaType == "DNA"], "D$", "")
alpha.df$ID[alpha.df$DnaType == "cDNA"] <- str_replace(alpha.df$Sample[alpha.df$DnaType == "cDNA"], "R$", "")


dna.alpha <- alpha.df[alpha.df$DnaType == "DNA",]
rna.alpha <- alpha.df[alpha.df$DnaType == "cDNA",]

############################

# calculate mean coordinates for duplicates
sum <- alpha.df %>%
  dplyr::group_by(ID, Rarefy, DnaType) %>%
  dplyr::summarise(Shannon = mean(Shannon, na.rm = T),
                   Simpson = mean(Simpson, na.rm = T),
                   Pielou = mean(Pielou, na.rm = T),
                   Chao1 = mean(Chao1, na.rm = T)) %>%
  ungroup()

setDT(sum)
temp <- dcast(sum, ID + Rarefy ~ DnaType, value.var = c("Simpson","Shannon","Pielou","Chao1"))

temp <- melt(temp, id.vars = c("ID","Rarefy"),
     variable.name = "Index",
     value.name = "Diversity") %>%
  separate(Index, into = c("Index","DnaType"), sep = "_")

temp <- dcast(temp, ID + Rarefy + Index ~ DnaType, value.var = "Diversity") 


ggplot(temp, aes(x = cDNA, y = DNA)) +
  theme_pubr() + 
  geom_point() +
  facet_grid(Rarefy~Index, scales = "free")

#

dna.alpha <- melt(setDT(dna.alpha), id.vars = c("ID","Rarefy"),
                  measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
                  variable.name = "Index",
                  value.name = "Diversity")

# merge with distance
plot.df <- merge(dna.alpha, dist.dr %>% dplyr::select(ID, sample.type.year,Year, Season, distance), by = "ID")
plot.df$log.dist <- log(plot.df$distance)
plot.df$log.div <- log(plot.df$Diversity)

colvec <- c("red4","chocolate3","orangered2","orange3",
            "cadetblue", "darksalmon",
            "darkolivegreen","darkolivegreen3",
            "royalblue","mediumorchid4", "violet","palevioletred2","navy","skyblue",
            "seagreen3")

p <- ggscatter(plot.df, x = "distance", y = "Diversity", 
               color = "sample.type.year",
               palette = colvec,
               xlab = "Distance between DNA and RNA \nin PCoA space (Bray Curtis)",
               ylab = "DNA Alpha Diversity",
               legend.title = "Sample Type")
(pf <- facet(p, facet.by = c("Index","Rarefy"), ncol = 3, nrow = 4, scales = "free")+
  stat_cor(method = "pearson", label.x = 0.2, cor.coef.name = "r"))

ggsave("./Figures/General/All_DNARNA_distance_alphadiv.png", pf,
       width = 25, height = 18, unit = "cm")

#(ap <- ggplot(plot.df, 
#              aes(x = log(distance), y = log(Diversity), colour = sample.type.year)) +
#    theme_pubr() +
#    geom_point() + 
#    facet_grid(Index~Rarefy, scales = "free"))

rna.alpha <- alpha.df[alpha.df$DnaType == "cDNA",]

rna.alpha <- melt(setDT(rna.alpha), id.vars = c("ID","Rarefy"),
                  measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
                  variable.name = "Index",
                  value.name = "Diversity")

# merge with distance
plot.df <- merge(rna.alpha, dist.dr %>% dplyr::select(ID, sample.type.year,Year, Season, distance), by = "ID")
plot.df$log.dist <- log(plot.df$distance)
plot.df$log.div <- log(plot.df$Diversity)

colvec <- c("red4","chocolate3","orangered2","orange3",
            "cadetblue", "darksalmon",
            "darkolivegreen","darkolivegreen3",
            "royalblue","mediumorchid4", "violet","palevioletred2","navy","skyblue",
            "seagreen3")

p <- ggscatter(plot.df, x = "distance", y = "Diversity", 
               color = "sample.type.year",
               palette = colvec,
               xlab = "Distance between DNA and RNA \nin PCoA space (Bray Curtis)",
               ylab = "RNA Alpha Diversity",
               legend.title = "Sample Type")
(pf <- facet(p, facet.by = c("Index","Rarefy"), ncol = 3, nrow = 4, scales = "free")+
    stat_cor(method = "pearson", label.x = 0.2, cor.coef.name = "r"))

ggsave("./Figures/General/All_DNARNA_distance_alphadiv.png", pf,
       width = 25, height = 18, unit = "cm")

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
                                                "Wellwater","Stream", "Tributary",
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
                                                 "Wellwater","Stream", "Tributary",
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
