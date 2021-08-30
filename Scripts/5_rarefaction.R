# Re-do all analyses with rarefaction thresholds

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
                                             "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
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

# Rarefaction richness comparison -------------------------------------------------------------------
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