# Do some simulations to test our delta-distance approach

# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("phyloseq", # wrangling
              "plyr", "tidyverse", "data.table", # wrangling
              "ggpubr", "gridExtra", # plotting
              "vegan", "ade4",
              "doMC") # stats


### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

### Functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")

# Set seed for session and reproducibility of permutations
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

# 3. Simulation prep/exploration ---------------------------------------------------------------------------------------------------
# Check species abundance distribution of following habitats:

rel.df <- select_newest("./Objects", "201520162017_css_otu99_paper1_")
rel.df <- readRDS(
  paste0("./Objects/", rel.df))

setDT(rel.df)

# merge some sample types
rel.df$sample.type.year <- factor(rel.df$sample.type.year, levels = c("Soil","Sediment",
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

# calculate means by sample type
means <- rel.df[, .(mean.css = mean(css.reads, na.rm = T),
                    sd.css = sd(css.reads, na.rm = T)), by = .(sample.type.year, Season, dna_type, OTU)]
# order the abundances to make ranks
means <- means[mean.css != 0,] # remove all 0 observations
means <- means[order(mean.css, decreasing = T)]
means[, rank.abun := 1:.N, by = .(dna_type, Season, sample.type.year)]
means[, log.mean := log1p(mean.css), by = .(dna_type, Season, sample.type.year)]

# Let's plot
means %>%
  filter(Season != "Autumn" &
           sample.type.year %in% c("Reservoirs","Stream","Soilwater", "Tributary","Riverine \nLakes", "Lake", "Upriver")) %>%
  ggplot(., aes(x = rank.abun, y = log.mean)) +
  geom_line(aes(colour = dna_type), size = 1.5) +
  facet_grid(Season~sample.type.year)
#scale_colour_manual(values = colvec[names(colvec) %in% unique(means$sample.type.year)])

# 4. Simulation ---------------------------------------------------------------------------------------------------------
# From Mario's paper:
# Lennon JT, Muscarella ME, Placella SA, Lehmkuhl BK.
# How, when, and where relic DNA affects microbial diversity. mBio 2018; 9: e00637-18.

# Two scenarios in their paper:
# Main: Intact community
# 1. Less even relic community
# 2. More even relic community

# Let's do a similar approach where we use different SADs to explore the delta-distances

# Create SADs ----------------------------------------------------------------------------------------------------------
# So we will do the following:
# We will create 3 different "regional pools" mimicking different habitat types
# which vary in their evenness
# We hypothesize that:
# Steeper communities exemplify habitats of strong selection (= higher delta distances)
# as only a few taxa become very dominant

# What was the max number of species that we observed in our dataset?
ncol(otu_table(pb))
# We observed a regional species pool of 16322
# But we'll shrink it to 10000, like in Lennon 2018

# 'Normal' community, moderate steepness
gamma1 <- rlnorm(n=10000, meanlog = 1, sdlog = 0.9)
gamma1 <- gamma1[rev(order(gamma1))]

gamma2 <- rlnorm(n=10000, meanlog = 1, sdlog = 1.2)
gamma2 <- gamma2[rev(order(gamma2))]

# Intermediate steepness
gamma3 <- rlnorm(n=10000, meanlog = 1, sdlog = 1.8)
gamma3 <- gamma3[rev(order(gamma3))]

# Very steep community
gamma4 <- rlnorm(n=10000, meanlog = 1, sdlog = 2.7) # Less even, 1.8
gamma4 <- gamma4[rev(order(gamma4))]

# Create OTUs: 10000 OTUs
otus <- paste("OTU_", sprintf("%05d", seq(1:length(gamma1))), sep = "") 

# Initiate Communities ------------------------------------------------------------------------------------------------
# We will sample the DNA from these regional pools
# Then, we duplicate the randomly sampled DNA community, so that RNA is a derivative of the DNA community
# and then, we apply three different degrees of replacement
# Here, replacement means replacing existing OTUs in DNA with
# OTUs that are not in DNA, and we will re-sample the number of removed
# DNA reads from the regional species pool

# what's the mean library size?
mean(sample_df(pb)$LibrarySize)
# So, we will sample 25000 times

# Create community matrix
# We will have three degrees of replacement
# And for each SAD - replacement combination, we will do the exercise 9 times to get replicates
# to correct for randomness

# Create empty data.frame and populate with OTUs
commat <- data.frame(matrix(ncol = length(otus), nrow = 0))
colnames(commat) <- otus

# We will write a function to apply the sampling and replacement procedure to each treatment/replicate
simulate.dnarna <- function(otus, size, rep.size, prob){
  # sample DNA community from SAD
  dna <- data.frame(table(sample(otus, size, replace = T, prob = prob))) %>%
    select(OTU = Var1, Reads = Freq) %>%
    setDT() %>%
    dcast(., ... ~ OTU, value.var = "Reads") %>%
    .[,-1] # remove empty column
  
  # Replace OTUs for RNA
  # Remove a certain amount of taxa, given by replacement size (= rep.size)
  rna <- sample(dna, size = round(ncol(dna) * rep.size), replace = F)
  both <- bind_rows(dna, rna) # combine DNA and RNA
  
  # Re-sample from regional pool, take OTUs that have not been sampled
  not.sampled <- otus[!(otus %in% colnames(both))]
  prob.resample <- gamma1[which(!(otus %in% colnames(both)))]
  
  # how many reads have been removed in the random removal exercise?
  # the difference will be the number of reads, which are re-sampled from regional pool to achieve equal library size
  resample.size <- abs(diff(rowSums(both, na.rm = T)))
  resample <- data.frame(table(sample(not.sampled, size = resample.size, replace = T,
                                      prob = prob.resample))) %>%
    select(OTU = Var1, Reads = Freq) %>%
    setDT() %>%
    dcast(., ... ~ OTU, value.var = "Reads") %>%
    .[,-1] # remove empty column
  
  # create empty df, with all regional OTUs
  regional.df <- data.frame(matrix(ncol = length(otus), nrow = 0))
  colnames(regional.df) <- otus
  # merge with dna-rna preliminary data
  regional.df <- bind_rows(regional.df, both) %>% setDT()
  # add resampled OTUs
  regional.df <- bind_rows(regional.df, resample)
  
  # replace all NAs with 0
  regional.df[is.na(regional.df)] <- 0
  
  # merge the two rows, preliminary RNA and re-sampled RNA
  # then, remove the re-sampled row
  regional.df[2,] <- regional.df[2,] + regional.df[3,]
  regional.df <- regional.df[-3,]
  # Sanity check
  #rowSums(regional.df) #ok, equal library size
  return(regional.df) # return community matrix of DNA and RNA
  
}

# Run function ---------------------------------------------------------------------------------------------------
# create empty list to populate with data
com.ls <- list()

# run loop to populate list with different SAD and replacement values, 9 times
for(j in 1:9){ # replicate loop
  set.seed(j) # always change seed otherwise, all replicate iterations will be the same
  for(i in 1:16){
    # set different SADs for different iterations in 3x steps
    # e.g. 1-3 will be SAD1, 4-6 SAD3, 7-9 SAD3
    if(i < 5){
      prob <- gamma1; label <- "gamma1"
    } else if(i >= 5 & i < 9){
      prob <- gamma2; label <- "gamma2"
    } else if(i >= 9 & i < 13){
      prob <- gamma3; label <- "gamma3"
    } else {
      prob <- gamma4; label <- "gamma4"
    }
    
    # set different replacement values
    # First in each set is high, 2nd is medium and 3rd is low
    if(i %in% seq(from = 1, to = 16, by = 4)){
      rep.size <- (1/2); rep.label <- "1/2" 
    } else if(i %in% seq(from = 2, to = 16, by = 4)){
      rep.size <- (1/3); rep.label <- "1/3" 
    } else if(i %in% seq(from = 3, to = 16, by = 4)){
      rep.size <- (1/6); rep.label <- "1/6" 
    } else {
      rep.size <- (1/9); rep.label <- "1/9"
    }
    # We need unique site numbers for each iteration
    # Hence, we will update the site number with each run
    sites.per.round <- 16
    if(j > 1){
      site.no <- (sites.per.round * j) - (sites.per.round - i)
    } else {
      site.no <- i
    }
    
    # apply function
    com.ls[[site.no]] <- simulate.dnarna(otus, size = 25000, rep.size = rep.size, prob = prob) %>%
      mutate(site = paste0("site",site.no),
             SAD = label, replacement = rep.label, dna_type = c("DNA","RNA")) %>% # add meta data
      select(site:dna_type, everything()) # reorder columns
    
  }
}

# combine list into data frame
sim.df <- bind_rows(com.ls)

# save as rds, a run takes quite some time
# if your computer cannot run this loop, read in the output in the next line
saveRDS(sim.df, "./Objects/simdf.rds")

# Read simulation output
sim.df <- readRDS("./Objects/simdf.rds")

# Correct for phantom taxa
# first, melt OTUs into long format
melt.com <- melt.data.table(
  sim.df,
  id.vars = c("site","SAD","replacement", "dna_type"),
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

# change dna type levels
melt.com[, dna_type := factor(dna_type, levels = c("DNA", "RNA"))]

# cast into wide format so that DNA and RNA reads are next to each other
castdf <- dcast(melt.com, site+SAD+replacement+OTU ~ dna_type, value.var = "reads")

# fill NAs with 0, those are taxa not found in RNA or DNA
castdf[is.na(DNA), DNA := 0]
castdf[is.na(RNA), RNA := 0]

# overwrite all DNA observations == 0, where RNA > 0 with 1
castdf[RNA > 0 & DNA == 0, DNA := 1]
castdf[RNA > 0 & DNA == 0,] # check if overwrite was successful

# calculate some metrics
castdf[DNA > 0 & RNA > 0, n.cooc := .N, .(site, SAD, replacement)] # how many OTUs are shared between DNA and RNA?
castdf[DNA > 0 & RNA == 0, n.uni.dna := .N, .(site, SAD, replacement)] # how many OTUs are unique to DNA?
castdf[, diff := DNA-RNA] # mean abundance difference

# Calculate evenness
castdf[, diversity.dna := diversity(DNA), by = .(site, SAD, replacement)]
castdf[, diversity.rna := diversity(RNA), by = .(site, SAD, replacement)]
castdf[, evenness.dna := diversity.dna / log(specnumber(DNA)), by = .(site, SAD, replacement)]
castdf[, evenness.rna := diversity.rna / log(specnumber(RNA)), by = .(site, SAD, replacement)]

# merge into a meta dataframe and summarize for each site for later
mets <- castdf[, .(evenness.dna = unique(.SD$evenness.dna),
                   evenness.rna = unique(.SD$evenness.rna),
                   n.cooc = unique(.SD[!is.na(n.cooc),]$n.cooc),
                   n.uni.dna = unique(.SD[!is.na(n.uni.dna),]$n.uni.dna),
                   mean.ab.diff = mean(diff, na.rm = T)), by = .(site, SAD, replacement)]

# format back to long format
temp <- melt.data.table(castdf,
                        id.vars = c("site", "SAD","replacement", "OTU"),
                        measure.vars = c("DNA","RNA"),
                        variable.name = "dna_type",
                        value.name = "reads")

# do row merge/overwrite
mer <- melt.com[temp, cor.reads := i.reads, on = .(site, SAD, replacement, dna_type, OTU)]

# create ID
mer[, ID := paste(site, SAD, replacement, dna_type, sep = "_")]

# cast into wide format
fin.df <- dcast(mer, ID ~ OTU, value.var = "cor.reads")
# make final community matrix and give row names
fin.df <- as.matrix(setDF(fin.df, rownames = fin.df$ID)[,-1]) # remove first row with "ID" move into rownames

# save meta
mer[, site.no := as.numeric(str_extract(site, "\\-*\\d+\\.*\\d*"))] # extract site numbers
meta <- mer %>% select(ID, site, SAD, replacement, dna_type, site.no) %>% unique() # shrink to unique entries

# Plot rank abundance distributions
#r.df <- mer[meta, c("sim.ID", "site.no") := list(i.sim.ID, i.site.no), on = .(ID)]
r.df <- mer
r.df <- r.df[order(cor.reads, decreasing = T)]
r.df[, rank.abun := 1:.N, by = .(ID)]
r.df[, log.reads := log1p(cor.reads)]

pj <- mets %>% select(SAD, evenness.dna) %>% unique() %>%
  group_by(SAD) %>%
  summarize(mean.even = mean(evenness.dna))

r.df[, SAD := factor(SAD, levels = c("gamma1","gamma2","gamma3", "gamma4"),
                        labels = c(paste0("SAD1:\nJ = ", round(pj[1,2], digits = 2)),
                                   paste0("SAD2:\nJ = ", round(pj[2,2], digits = 2)),
                                   paste0("SAD3:\nJ = ", round(pj[3,2], digits = 2)),
                                   paste0("SAD4:\nJ = ", round(pj[4,2], digits = 2))))]
# Change to numbers, more meaningful than classifier
r.df[, replacement := factor(replacement, levels = c("1/2", "1/3", "1/6", "1/9"))]

# Let's plot
sad.p <- r.df %>%
  filter(site.no <= 16) %>%
  ggplot(., aes(x = rank.abun, y = log.reads, alpha = dna_type)) +
  theme_pubr() +
  geom_line(aes(colour = dna_type), size = 1.5) +
  facet_grid(replacement~SAD) +
  scale_alpha_manual(values = c(1,0.7)) +
  scale_colour_viridis_d(direction = -1, end = 0.9, option = "B", name = "Nucleic acid type") +
  scale_x_continuous(n.breaks = 3) +
  labs(x = "Rank", y = "Abundance (log-scale)") +
  guides(alpha = "none") +
  theme(strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(colour = "white", size = 10),
        strip.text.y = element_text(angle = 0, size = 12),
        legend.position = "top",
        panel.spacing.x = unit(1, "lines"),)
sad.p <- annotate_figure(sad.p, right = text_grob("Replacement", rot = 270))

ggsave("./Figures/Final/sim_sads.png", sad.p, 
       width = 18, height = 13, units = "cm", dpi = 300)

# Do multivariate analysis -------------------------------------------------------------------------------------------
# PCoA with Bray-Curtis
pb.bray <- vegdist(fin.df, method = "bray")
is.euclid(pb.bray) # FALSE
pb.bray <- sqrt(pb.bray) # make Euclidean
is.euclid(pb.bray)
# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray)

ncol(fin.df) # 10000 OTUs
nrow(fin.df) # 84 registers

ncol(pb.bray.pcoa$vectors) # 83 axes

axes <- c(1,2)
# extract scores and variance explained
pdataframe <- data.frame(Sample = as.character(row.names(pb.bray.pcoa$vectors)),
                         pb.bray.pcoa$vectors[,axes[1:length(axes)]],
                         stringsAsFactors = F)  # extract site scores

pb.var <- data.frame(Axes = axes,
                     var = round(100 * pb.bray.pcoa$values$Eigenvalues[axes] / sum(pb.bray.pcoa$values$Eigenvalues), 2),
                     stringsAsFactors = F)

# merge with a selection of meta data
pdataframe <- merge(pdataframe, meta, by.x = "Sample", by.y = "ID")

# main plot
(p <- ggplot(pdataframe, aes_string(x = paste0("Axis.", axes[1]), 
                                    y = paste0("Axis.", axes[2]))) +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    geom_point(aes(shape = dna_type, colour = replacement), size = 2, alpha = 0.5) + 
    coord_fixed(1) + # ensure aspect ratio
    labs(x = paste0("PC",  axes[1]," (", round(pb.var[pb.var$Axes == axes[1],"var"],
                                               digits = 1),"%)"), 
         y = paste0("PC",  axes[2]," (", round(pb.var[pb.var$Axes == axes[2],"var"],
                                               digits = 1),"%)")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

# Calculate Sorensen -----------------------------------------------------------------------------------------
# Calculate incidence based dissimilarity
# PCoA with Sorensen, incidence based equivalent of Bray Curtis
pb.soren <- vegdist(fin.df, method = "bray", binary = T)
is.euclid(pb.soren) # FALSE
pb.soren <- sqrt(pb.soren) # make Euclidean
is.euclid(pb.soren)
# make PCoA
pb.soren.pcoa <- ape::pcoa(pb.soren)
# how many axes to have 75% of variance captured?
nrow(pb.soren.pcoa$values[pb.soren.pcoa$values$Cumul_eig <= 0.75,])
# 50 axes

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

bray.df <- setDT(merge(meta, bray.df, by.y = "Sample", by.x = "ID"))
soren.df <- setDT(merge(meta, soren.df, by.y = "Sample", by.x = "ID"))

# melt and combine bray-curtis and sorensen results into one data frame
temp.75 <- rbind(melt.data.table(bray.df, id.vars = c("ID","site", "SAD", "replacement",
                                                      "dna_type", "Metric"), measure.vars = patterns("^Axis."),
                                 variable.name = "Axis", value.name = "Coordinates"),
                 melt.data.table(soren.df, id.vars = c("ID","site", "SAD", "replacement",
                                                       "dna_type", "Metric"), measure.vars = patterns("^Axis."),
                                 variable.name = "Axis", value.name = "Coordinates"))
setDT(temp.75)

temp.75 <- dcast(temp.75, site + SAD + replacement + Axis + Metric ~ dna_type, value.var = c("Coordinates"))
# remove NAs
temp.75 <- na.omit(temp.75)
temp.75[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
temp.75 <- temp.75[, .(sum.dist = sum(pnt.dist)), by = .(Metric, site, SAD, replacement)] # sum temp axes
temp.75 <- temp.75[, dist := sqrt(sum.dist)] # take sqrt

dist.75 <- temp.75[meta, c("site.no") := list(i.site.no), on = .(site)]

# Calculate delta of the two metrics
diff.df <- dcast(dist.75, site + SAD + replacement + site.no ~ Metric, value.var = "dist")
diff.df <- diff.df[, delta := Bray - Sorensen]

# merge some more variables
diff.df <- merge(diff.df, mets %>% select(site,evenness.dna:mean.ab.diff), by = "site")

# Calculate mean of replicates
diff.df <- diff.df[, .(Bray = mean(Bray, na.rm = T),
                       Sorensen = mean(Sorensen, na.rm = T),
                       delta = mean(delta, na.rm = T),
                       sd.delta = sd(delta, na.rm = T),
                       evenness = mean(evenness.dna, na.rm = T),
                       n.cooc = mean(n.cooc, na.rm = T),
                       n.uni.dna = mean(n.uni.dna, na.rm = T),
                       mean.ab.diff = mean(mean.ab.diff, na.rm = T),
                       sd.uni.dna = sd(n.uni.dna, na.rm = T),
                       sd.ab.diff = sd(mean.ab.diff, na.rm = T)),
                   by = .(SAD, replacement)]

diff.df[, SAD := factor(SAD, levels = c("gamma1","gamma2","gamma3", "gamma4"),
                     labels = c(paste0("SAD1:\nJ = ", round(pj[1,2], digits = 2)),
                                paste0("SAD2:\nJ = ", round(pj[2,2], digits = 2)),
                                paste0("SAD3:\nJ = ", round(pj[3,2], digits = 2)),
                                paste0("SAD4:\nJ = ", round(pj[4,2], digits = 2))))]
# Change to numbers, more meaningful than classifier
diff.df[, replacement := factor(replacement, levels = c("1/2", "1/3", "1/6", "1/9"))]

# melt
# diff.df <- melt(diff.df, id.vars = c("sample","sim.ID"),
#                 measure.vars = c("delta","Sorensen","Bray"),
#                 variable.name = "Metric", value.name = "dist")
diff.df[,plot.id := paste(SAD, replacement)]
diff.df[, plot.id := factor(plot.id,
                            levels = c(paste0("SAD1:\nJ = ", round(pj[1,2], digits = 2), " 1/2"),
                                       paste0("SAD1:\nJ = ", round(pj[1,2], digits = 2), " 1/3"),
                                       paste0("SAD1:\nJ = ", round(pj[1,2], digits = 2), " 1/6"),
                                       paste0("SAD1:\nJ = ", round(pj[1,2], digits = 2), " 1/9"),
                                       paste0("SAD2:\nJ = ", round(pj[2,2], digits = 2), " 1/2"),
                                       paste0("SAD2:\nJ = ", round(pj[2,2], digits = 2), " 1/3"),
                                       paste0("SAD2:\nJ = ", round(pj[2,2], digits = 2), " 1/6"),
                                       paste0("SAD2:\nJ = ", round(pj[2,2], digits = 2), " 1/9"),
                                       paste0("SAD3:\nJ = ", round(pj[3,2], digits = 2), " 1/2"),
                                       paste0("SAD3:\nJ = ", round(pj[3,2], digits = 2), " 1/3"),
                                       paste0("SAD3:\nJ = ", round(pj[3,2], digits = 2), " 1/6"),
                                       paste0("SAD3:\nJ = ", round(pj[3,2], digits = 2), " 1/9"),
                                       paste0("SAD4:\nJ = ", round(pj[4,2], digits = 2), " 1/2"),
                                       paste0("SAD4:\nJ = ", round(pj[4,2], digits = 2), " 1/3"),
                                       paste0("SAD4:\nJ = ", round(pj[4,2], digits = 2), " 1/6"),
                                       paste0("SAD4:\nJ = ", round(pj[4,2], digits = 2), " 1/9")))]

delta<-ggplot(diff.df, aes(x = SAD)) +
  theme_pubr() +
  geom_errorbar(aes(ymin = delta - sd.delta, ymax = delta + sd.delta,
                    group = plot.id), position = position_dodge(0.7), width = 0.2) +
  geom_jitter(aes(y = delta, fill = replacement, shape = SAD), position = position_dodge(0.7), size = 3) +
  scale_fill_viridis_d(option = "E", name = "Replacement", direction = -1) +
  scale_shape_manual(values = c(21,23,25,24)) +
  labs(x = NULL, y = expression(paste(Delta," Distance"))) +
  guides(fill = guide_legend(order = 2, override.aes=list(shape = 21)),
         shape = guide_legend(order = 1))

# Test a few hypotheses
# delta distance is negatively driven by how many OTUs are unique to DNA
uni.lm <- lm(diff.df$n.uni.dna ~ diff.df$delta)
#plot(uni.lm)
summary(uni.lm)

lm.eq <- function(x){
  lm_coef <- list(a = as.numeric(round(summary(x)$coefficients[1], digits = 2)),
                  b = as.numeric(round(summary(x)$coefficients[2], digits = 2)),
                  r2 = round(summary(x)$r.squared, digits = 2),
                  pval = abbrev.p(anova(x)$`Pr(>F)`[1])[1]);
  if(lm_coef$b < 0){
    lm_coef$b <- abs(lm_coef$b)
    lm_eq <- c(substitute(italic(y) == a - b %.% italic(x), lm_coef),
                             substitute(~~italic(R)^2~"="~r2*","~~italic(p)~pval, lm_coef))
  } else{
    lm_eq <- c(substitute(italic(y) == a + b %.% italic(x), lm_coef),
               substitute(~~italic(R)^2~"="~r2*","~~italic(p)~pval, lm_coef))
  }
  
  return(as.character(as.expression(lm_eq)))
}


a <- ggplot(diff.df, aes(x = delta, n.uni.dna)) +
  theme_pubr() +
  geom_smooth(method = 'lm', colour = "black", size = 0.5) +
  geom_point(aes(fill = replacement, shape = SAD), size = 3) +
  annotate(geom = "text",
           y = 6300,
           x = 0.3,
           label = lm.eq(uni.lm)[1],
           parse = T) +
  annotate(geom = "text",
           y = 5900,
           x = 0.3,
           label = lm.eq(uni.lm)[2],
           parse = T) +
  scale_fill_viridis_d(option = "E", name = "Replacement", direction = -1) +
  scale_shape_manual(values = c(21,23,25,24)) +
  theme(legend.position = "right") +
  labs(x = expression(paste(Delta, " Distance")),
       y = "Number of OTUs not in RNA\n(= unreactive OTUs)") +
  guides(fill = guide_legend(order = 2, override.aes=list(shape = 21)),
         shape = guide_legend(order = 1))

# Delta distance is positively driven by abundance differences between DNA and RNA
ab.lm <- lm(diff.df$mean.ab.diff ~ diff.df$delta)
#plot(ab.lm)
summary(ab.lm)

b <- ggplot(diff.df, aes(x = delta, mean.ab.diff)) +
  theme_pubr() +
  geom_smooth(method = 'lm', colour = "black", size = 0.5) +
  geom_point(aes(fill = replacement, shape = SAD), size = 3) +
  annotate(geom = "text",
           y = 0.25,
           x = 0.3,
           label = lm.eq(ab.lm)[1],
           parse = T) +
  annotate(geom = "text",
           y = 0.22,
           x = 0.3,
           label = lm.eq(ab.lm)[2],
           parse = T) +
  scale_fill_viridis_d(option = "E", name = "Replacement", direction = -1) +
  scale_shape_manual(values = c(21,23,25,24)) +
  theme(legend.position = "right") +
  labs(x = expression(paste(Delta, " Distance")),
       y = "Mean reads difference\nbetween DNA and RNA") +
  guides(fill = guide_legend(order = 2, override.aes=list(shape = 21)),
         shape = guide_legend(order = 1))

out <- ggarrange(delta,a,b, ncol = 3, align = 'hv', common.legend = T, labels = "auto")

ggsave("./Figures/Final/sim_delta.png", out, 
       width = 28, height = 10, units = "cm", dpi = 300)
