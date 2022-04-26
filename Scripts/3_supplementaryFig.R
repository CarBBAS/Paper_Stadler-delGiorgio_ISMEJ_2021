#-- Script for the publication:
#-- Terrestrial connectivity, upstream aquatic history and seasonality shape bacterial
#-- community assembly within a large boreal aquatic network. The ISME Journal.
#-- Authors: Masumi Stadler & Paul A. del Giorgio
#-- Responsible code author: Masumi Stadler

# This script is the fourth of a series of scripts that were used to analyse the data
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
              "plyr", "tidyverse", "data.table", # wrangling
              "ggpubr", "gridExtra", # plotting
              "vegan", "ade4",
              "doMC") # stats


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

### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# 2. Read and prepare data ---------------------------------------------------------------

# do we have several files per object? -> take newest version
# read in
fin.df <- select_newest("./Output", "201520162017_fin_css_otu99_phantomcor_paper1_")
fin.df <- as.matrix(read.csv(
  paste0("./Output/", fin.df),
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))

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

# Update phyloseq object
pb <- phyloseq(otu_table(fin.df, taxa_are_rows = F),
               sample_data(met.df),
               tax_table(tax.tab))

# remove THW1D - no meta data
pb <- subset_samples(pb, seq_name != "THW1D")
pb <- prune_taxa(taxa_sums(pb) != 0, pb)
pb <- prune_samples(sample_sums(pb) != 0, pb) # remove samples with no reads

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

# set theme for plotting
theme_set(theme_bw())


# Fig. S2: Rarefaction richness comparison -------------------------------------------------------------------
# Do not run this code if you have a slow computer, or if you want to go quickly through the script
# Read in intermediate data frame in line 214 onwards
df <- sample_df(pb)
hist(df$LibrarySize, breaks = 20)

#choose rarefaction thresholds
min_lib <- c(min(df$LibrarySize), 5000, 10000, 20000, 30000)
set.seed(3)
rar.ls <- list()
for(i in 1:length(min_lib)){
  out <- rarefy_even_depth(pb, sample.size = min_lib[i], rngseed = 3)
  rar.ls[[i]] <- out
  names(rar.ls)[i] <- paste("rar",min_lib[i], sep = "_")
}

# amend original matrix
rar.ls[[length(rar.ls)+1]] <- pb
names(rar.ls)[length(rar.ls)] <- "orig"

alpha <- llply(rar.ls, function(x){
  rar <- as.data.frame(otu_mat(x))
  data.table::setDT(rar, keep.rownames = "Sample")
  rar <- setDF(rar)
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
  alpha <- plyr::ddply(rar, ~ Sample, function(z){
    data.frame(Shannon = vegan::diversity(z[,-1], index = "shannon"),
               Simpson = vegan::diversity(z[,-1], index = "simpson"),
               Pielou = vegan::diversity(z[,-1], index = "shannon") / log(sum(z[,-1] > 0))
               )#Chao1 = vegan::estimateR(z[,-1],)["S.chao1",]
  })
  
  return(alpha)
}, .parallel = T)

#saveRDS(alpha, "./Objects/alpha_rarefaction.rds")
alpha <- readRDS("./Objects/alpha_rarefaction.rds")
  
#names(alpha) <- c(paste0("lib", min_lib), "css")
alpha.df <- bind_rows(alpha, .id = "Data")
setDT(alpha.df); setDT(met.df)
# combine with DR.names and meta data
alpha.df[met.df, c("Season","sample.type.year","dr_match_name","dna_type") := 
           list(i.Season, i.sample.type.year, i.dr_match_name, i.dna_type), on = c('Sample' = 'seq_name')]
#sumdf <- readRDS("./Objects/summary.meta.with.oldnames.rds")

write.table(alpha.df, paste0("./Output/alpha_div_otu99_summary_", Sys.Date(),".csv"), sep = ",", dec = ".", row.names = F)

#-----------------------------------------------------------------------------------------------
# Do not run above code again, takes time
# read in processed data frame

alpha.df <- select_newest("./Output/", "alpha_div_otu99_summary_")
alpha.df <- read.csv(paste0("./Output/", alpha.df), sep = ",", stringsAsFactors = F)

# make to data table
setDT(alpha.df)

alpha.df$Data <- factor(alpha.df$Data, levels = c("orig", "rar_1470", "rar_5000", "rar_10000", "rar_20000",
                                                  "rar_30000"),
                          labels = c("CSS", "Rarefied: Lib1470", "Rarefied: Lib5000","Rarefied: Lib10000",
                                     "Rarefied: Lib20000","Rarefied: Lib30000"))

melt.alpha <- melt(alpha.df, id.vars = c("dr_match_name", "Data", "dna_type", "sample.type.year","Season"),
                   measure.vars = c("Shannon","Simpson","Pielou"), #"Chao1"
                   variable.name = "Diversity",
                   value.name = "Index")

melt.alpha$Season <- factor(melt.alpha$Season, levels = c("Spring", "Summer","Autumn"))


# make plot with ggpubr to include pearson's correlation outputs directly in the plot
(rar.box <- ggplot(melt.alpha, aes(x = dna_type, y = Index)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(Diversity ~ Data, scales = "free") +
  labs(x = "Nucleic acid type", y = expression(paste(alpha," - Diversity")))+
  theme(strip.background = element_rect(fill= "white"),
        panel.grid = element_blank()))

ggsave("./Figures/Final/alphadiv_rar_comp_boxplots.png", rar.box,
       width = 20, height = 12, unit = "cm")

cast.alpha <- melt.alpha %>%
  group_by(dr_match_name, Data, Diversity, dna_type) %>%
  dplyr::summarise(Season = unique(Season),
            Index = mean(Index, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = dna_type, values_from = Index)

(rar.lin <- ggplot(cast.alpha,
                  aes(x = RNA, y = DNA, fill = Season)) +
  geom_point(shape = 21) +
  facet_wrap(Diversity~Data, scales = "free", ncol = 6) +
  scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) +
    theme(strip.background = element_rect(fill= "white"),
          panel.grid = element_blank()))
# colour-blind friendly

ggsave("./Figures/Final/alphadiv_rar_comp_scatter.png", rar.lin,
       width = 26, height = 15, unit = "cm")

# Together
p <- ggarrange(rar.box, rar.lin, nrow = 2, labels = "auto", common.legend = T, legend = "bottom")

ggsave("./Figures/Final/alphadiv_comp.png", p,
       width = 26, height = 30, unit = "cm")

# Fig. S7: Abundance classification ---------------------------------------------------------------------------------
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

# Get derivatives for an example and plot for supplementary material
x <- means[sample.type.year == "Soilwater" & dna_type == "DNA"]

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
       width = 25, height = 9, unit = "cm")




# Fig. S8: Taxonomic composition ----------------------------------------------------------------------------------
# We want to show the taxonomic composition of our samples plus the abundance
# As we have too many samples, best would probably be to calcualte the mean abundance for each group
# Groups are: sample.type.year + dna_type + Season

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

length(levels(factor(tax.df$phylum)))
length(levels(factor(tax.df$class)))
length(levels(factor(tax.df$order)))
length(levels(factor(tax.df$family)))
length(levels(factor(tax.df$genus)))

# add meta data
meta <- sample_df(pb) %>%
  dplyr::select(dna_type, year, Season, sample.type.year, LibrarySize)
setDT(meta, keep.rownames = "Sample") # make data.table

pb.df <- pb.df[meta, on = .(Sample)] # merge

# calculate mean abundance (css) for each category
mean.pb <- pb.df[, .(css.mean = mean(css, na.rm = T),
                     css.sd = sd(css, na.rm = T)),
                 by = .(OTU, sample.type.year, Season, dna_type)][
                   , css.sum := sum(css.mean, na.rm = T),
                   by = .(sample.type.year, Season, dna_type)
                 ][, css.rel := css.mean * 1 / css.sum]
mean.pb <- mean.pb[tax.df, on = .(OTU)] # merge with taxonomy data
mean.pb[, ID := paste(dna_type, Season, sample.type.year, sep = "_")] # add plot ID

# mean of means for results section
temp <- mean.pb[, .(sum = sum(css.rel, na.rm = T)), .(sample.type.year, phylum, dna_type, Season)]
temp <- temp[, .(mean = mean(sum) *100, sd = sd(sum)), .(phylum)]
temp[order(mean, decreasing = T),]

# calculate mean library size per category
lib.size <- pb.df[, .(lib.mean = mean(LibrarySize, na.rm = T),
                      lib.sd = sd(LibrarySize, na.rm = T)),
                  by = .(sample.type.year, dna_type, Season)]
lib.size[, ID := paste(dna_type, Season, sample.type.year,sep = "_")] # add plot ID

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

# GTDB has some "sub" phyla separeted with "_"
mean.pb[, phylum := sapply(strsplit(phylum, split = "_"),  "[", 1)]
mean.pb[, phylum := factor(phylum)]

# colour blind friendly, derived from https://medialab.github.io/iwanthue/
col_vector<-rev(c("#350070","#06a644","#dcd873","#b9ce40","#6b8aff","#f5d249","#ec83f6","#7beb8b", "#003499",
              "#c947b1","#01b072","#df3587","#006f26","#ff83d0","#215a00","#99a1ff","#668200",
              "#015db8","#e77e28","#019cf8","#b5221d","#67b8ff","#bd0a35","#b2e294","#840066",
              "#314800","#ffb0ed","#954700","#3d0e52","#ff9b61","#59003b","#ff6e83","#aa7dbf","#620009",
              "#b32d54","#36dee6","#d54a4a","#7cbf76","#c86f2d","#dd89c2","#472f80","#86a73c",
              "#6f81ea", "#beb337","#3f7f25","#df70c7","#4cbe84","#d14184"))
df<- data.frame(col_vector = col_vector, i = seq(1:length(col_vector)))
ggplot(df, aes(y = i, x = i, colour = levels(factor(mean.pb$phylum)), label = levels(factor(mean.pb$phylum)))) +
  geom_text() +
  scale_colour_manual(values = col_vector)

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
      scale_fill_manual(name = "Phylum", values = col_vector) +
      scale_colour_manual(name = "Phylum", values = col_vector) +
      theme(legend.position = "right") +
      guides(fill = guide_legend(nrow = length(levels(mean.pb$phylum))/2, byrow =F))
  )


tax <- ggplot(mean.pb, 
              aes(x = ID, y = css.rel, fill = phylum, colour = phylum)) +
  theme_pubr() +
  geom_bar(stat = "identity", width = 0.8) +
  guides(fill = "none", colour = "none") +
  labs(y = "Average relative abundance") +
  scale_fill_manual(values = col_vector) +
  scale_colour_manual(values = col_vector) +
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
ggsave("./Figures/Final/Tax_LibSiz_Phyla_otus_gtdb.png", combo,
       width = 300, height = 203, unit = "mm", dpi = 400)


# Fig. S9: Only terrestrial PCoA -------------------------------------------------------------------------------------

ter <- subset_samples(pb, sample.type.year == "Soil" |
                        sample.type.year == "Sediment" |
                        sample.type.year == "Soilwater" |
                        sample.type.year == "Groundwater")

# remove ASVs that do not appear in this dataset
ter <- prune_taxa(taxa_sums(ter) != 0, ter)
ter <- prune_samples(sample_sums(ter) != 0, ter)

pb.mat <- otu_mat(ter)
pb.mat <- decostand(pb.mat, "hellinger")
#pb.mat <- log2(pb.mat + 1)

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

# Make hulls around DNA and RNA
find_hull <- function(x, axes){x[chull(x[,paste0("Axis.",axes[1])], x[,paste0("Axis.",axes[2])]),]}
hulls <- ddply(dna.pcoa$df, "dna_type", find_hull, axes = c(1,2))

p <- ggplot(dna.pcoa$df, aes(x = Axis.1, 
                                     y = Axis.2)) +
  theme_cust(base_theme = "pubr",
             border = T) +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
  scale_fill_manual(values = colvec[names(colvec) %in% as.character(levels(dna.pcoa$df$sample.type.year))],
                    name = "Habitat Type") +
  geom_point(aes(fill = sample.type.year, shape = Season), size=2.5) + #colour = dna_type,
  geom_polygon(data = hulls, alpha = 0, aes(linetype = dna_type), fill = "white", colour = "black") +
  scale_linetype_manual(values = c("solid","dotted"), name = "Nucleic Acid Type") +
  coord_fixed(1) + # ensure aspect ratio
  scale_shape_manual(values = c(21,23,25)) +
  #scale_size_manual(values = c(2.5, 2.6), name = "Nucleic Acid \nType") +
  labs(x = paste0("PC1 (", round(dna.pcoa$var[1,"var"],
                                 digits = 1),"%)"), 
       y = paste0("PC2 (", round(dna.pcoa$var[2,"var"],
                                 digits = 1),"%)")) +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)),
         size = FALSE)

ggsave("./Figures/Final/terr_PCoA_gtdb_hulls.png", p,
       width = 12, height = 10, unit = "cm")

# Fig. S10: RNA PCoA -------------------------------------------------------------------------------------
# subset only RNA samples
rna <- subset_samples(pb, dna_type == "RNA")

# remove ASVs that do not appear in this dataset
rna <- prune_taxa(taxa_sums(rna) != 0, rna)
rna <- prune_samples(sample_sums(rna) != 0, rna)

# extract ASV matrix
pb.mat <- otu_mat(rna)
pb.mat <- decostand(pb.mat, "hellinger")
#pb.mat <- log2(pb.mat + 1)

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
p <- rna.pcoa$plot

# save
ggsave(paste0("./Figures/Final/PCoA_RNA_SampleType_gtdb.png"),  p,
       width = 12, height = 10, unit = "cm")

# Extra --------------------------------------------------------------------------------------------
# Fig. Sx -------------------------------------------------------------------------------------
# DNA ~ RNA richness
alpha.df <- select_newest("./Output/", "alpha_div_otu99_summary_")
alpha.df <- read.csv(paste0("./Output/", alpha.df), sep = ",", stringsAsFactors = F)

# make to data table
setDT(alpha.df)

# keep CSS results only
temp <- alpha.df[dr_match_name %in% unique(alpha.df[dna_type == "RNA",]$dr_match_name),][Data == "orig",]

# keep Shannon and Pielou, we're not very interested in the other alpha indices
cast.alpha <- dcast(temp,
                    dr_match_name + Season + sample.type.year ~ dna_type,
                    value.var = c("Shannon", "Pielou"))

# and melt
melt.alpha <- melt(cast.alpha, id.vars = c("dr_match_name"),
                   measure.vars = c("Shannon_DNA",
                                    "Shannon_RNA",
                                    "Pielou_DNA",
                                    "Pielou_RNA"),
                   value.name = "Index") %>%
  separate(col = "variable", into = c("Diversity","dna_type"), sep = "_") %>%
  drop_na() %>%
  dcast(dr_match_name + Diversity ~ dna_type, value.var = "Index")

melt.alpha[cast.alpha, c("sample.type.year", "Season") :=
             list(i.sample.type.year, i.Season), on = .(dr_match_name)]

# Visual inspection
ggscatter(
  melt.alpha, x = "RNA", y = "DNA",
  color = "Season", add = "reg.line"
)+
  facet_wrap("Diversity", scales = "free") +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Season)
  )

# the lack of correlation/dependency in all seasons except summer is interesting

# calculate linear models for alpha diversity indices
alpha.mod <- dlply(melt.alpha, ~ paste(Diversity, Season, sep = "_"), function(df) {
  setDF(df)
  lm(DNA ~ sqrt(RNA), data = df, na.action = na.omit)
})

# check model assumptions
# Homogeneity of variances
(diagn.plots <- llply(alpha.mod, function(list){
  par(mfrow=c(2,2))
  plot(list)
  p <- recordPlot()
  plot.new()
  return(p)
}))

# there is no point in getting assumptions right, as there is probably no relationship
# in autumn and spring

# extract P-value
(val <- ldply(alpha.mod, anova))
(val <- val[!is.na(val$`Pr(>F)`),]) # keep only the rows that have the P-value
colnames(val)[c(1,6)] <- c("Div_Season","Pval") # rename P-value column


# Pielou autumn significant, BUT assumptions not fulfilled.
# Transformations are not helping.... We include it in the results but interpret with caution

melt.alpha[, show := F][Diversity == "Pielou" & (Season == "Summer" | Season == "Autumn"), show := T]
melt.alpha[Diversity == "Shannon" & Season == "Summer", show := T]

melt.alpha$Season <- factor(melt.alpha$Season, levels = c("Spring","Summer","Autumn"))
melt.alpha$Diversity <- factor(melt.alpha$Diversity, levels = c("Shannon","Pielou"))

library(ggpmisc)

(alpha.p <- ggplot() +
    theme_pubr() +
    geom_point(data = melt.alpha, aes(x = RNA, y = DNA, fill = Season), shape = 21) +
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    geom_smooth(data = melt.alpha[show == T,] %>% drop_na(), 
                mapping = aes(x = RNA, y = DNA,
                              colour = Season), method = "lm", se = F) +
    facet_wrap("Diversity", scales = "free") +
    stat_poly_eq(formula = y ~ x,
                 data = melt.alpha[show == T,] %>% drop_na(), 
                 mapping = aes(x = RNA, y = DNA,
                               label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~"),
                               colour = Season), 
                 parse = TRUE,
                 label.x = "right", label.y = "bottom",
                 vstep = 0.05, size = 2.5) +
    labs(x = "RNA") +
    guides(colour = "none") +
    scale_colour_manual(values = c("gold", "#D55E00")) +
    theme(strip.background = element_rect(fill= "white"),
          panel.grid = element_blank()))


ggsave("./Figures/Final/alpha.div.dna.rna.lm.png", alpha.p,
       width = 15, height = 10, units = "cm")

# Fig. Sx3 --------------------------------------------------------------------------------------------
# Check richness along continuum

# merge some sample types
melt.alpha$sample.type.year <- factor(melt.alpha$sample.type.year, levels = c("Soil","Sediment",
                                                                              "Soilwater",
                                                                              "Stream", "Tributary",
                                                                              "Riverine \nLakes", "Lake", 
                                                                              "Upriver", "Reservoirs", "Downriver",
                                                                              "Estuary"),
                                      labels = c("Soil","Sediment",
                                                 "Soilwater",
                                                 "Stream", "Tributary",
                                                 "Riverine Lakes", "Lake", 
                                                 "Upriver", "Reservoirs", "Downriver",
                                                 "Estuary"))

melt.alpha <- melt(melt.alpha, id.vars = c("sample.type.year","Season","Diversity","dr_match_name"),
                   measure.vars = c("DNA","RNA"), variable.name = "dna_type", value.name = "Values")

melt.alpha$dna_type <- factor(melt.alpha$dna_type, levels = c("DNA","RNA"),
                              labels = c("DNA", "RNA"))


(shan <- ggplot(melt.alpha[Diversity == "Shannon" & Season != "Autumn",],
                aes(x = sample.type.year, y = Values, fill = sample.type.year)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(name = "Habitat type", values = colvec) +
    facet_wrap(dna_type~Season, scales = "free_y") +
    labs(x = "Habitat type", y = "Shannon") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill= "white"),
          panel.grid = element_blank(),
          axis.title.x = element_blank()))

(pie <- ggplot(melt.alpha[Diversity == "Pielou" & Season != "Autumn",],
               aes(x = sample.type.year, y = Values, fill = sample.type.year)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(name = "Habitat type", values = colvec) +
    facet_wrap(dna_type~Season, scales = "free_y") +
    labs(x = "Habitat type", y = "Pielou") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill= "white"),
          panel.grid = element_blank()))

(p <- ggarrange(shan, pie,
                nrow = 2,  common.legend = T,
                align = "hv"))

ggsave("./Figures/Final/continuum_alpha.png", p,
       width = 18, height = 23, units = "cm")