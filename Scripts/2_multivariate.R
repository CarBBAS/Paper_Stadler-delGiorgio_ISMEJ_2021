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

# get number of samples for methods


############
# Analysis #
############
# find abundant, moderate and rare taxa by sample type
# melt into long format
rel.df <- select_newest("./Output", "201520162017_css")
rel.df <- read.csv(
  paste0("./Output/", rel.df),
  sep = "\t",
  dec = ".",
  stringsAsFactors = F
)


molten.asv <- melt.data.table(setDT(as.data.frame(asv.tab), keep.rownames = "ASV"),
                id.vars = "ASV", # skip measure.var, takes all columns
                variable.name = "Sample",
                value.name = "css.reads")

# merge relevant meta data
molten.asv[sample_df(met.df), c("sample.type.year",
                                "Season") := list(i.sample.type.year,
                                                  i.Season), on = c("Sample" = "DadaNames")]

# calculate relative abundance
molten.asv[, total := sum(css.reads), .(Sample)]
molten.asv[, rel.abund := css.reads / total]


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
                                                "Groundwater","Stream", "Tributary",
                                                "Headwater \nLakes", "Upstream \nPonds", "Lake", "Lake",
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
          sample.type.year == "Headwater \nLakes" |
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
    geom_boxplot(width = 0.7, outlier.size = 0.5, size = 0.3) + # outlier.alpha = 0
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
    geom_boxplot(width = 0.7, outlier.size = 0.5, size = .3) + # outlier.alpha = 0
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
    geom_boxplot(width = 0.7, outlier.size = 0.5, size = 0.3) + # outlier.alpha = 0
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
    geom_boxplot(width = 0.7, outlier.size = 0.5, size = .3) + # outlier.alpha = 0
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
ggsave("./Figures/Final/All_DNARNA_withinPCoA_distance.tiff", 
       p, width = 18, height = 13, unit = "cm")

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
