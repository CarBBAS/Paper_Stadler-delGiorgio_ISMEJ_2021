###-----------------------------------------###
#-   Pre-processing of data for Chapter 1   - #
###-----------------------------------------###
###-----------------------------------###
#-   Binning and CSS transformation   - #
###-----------------------------------###
# This script processes the meta data of the whole La Romaine project to extract the necessary information
# for the first paper of this thesis.

# First we take all samples except those from Bioassays and 2018 (as RNA data is not avaialable at this moment)

# Second, we split the data set into ID groups (Year, Season, DnaType, ASV) and evaluate whether
# there is an actual observation or if that particular ASV is absent in that ID group.
# To do this, we will split the x-axis which is catchment.area into bins. If within a bin,
# there are several samples but only one actual observation then we will remove the observation
# if the number of reads is below 10. This function consumes a lot of memory and CPU power,
# thus we only execute the bin_qual_control() function on the actually present IDs.

# Third, we will transform the corrected counts with metgenomeSeq to CSS ().

# Finally, we calculate the relative abundances based on the css read numbers for completeness.
#setwd("/media/shared/Documents/University/PhD/Analyses/Molecular/lr.chapter1")
setwd("/home/bioinf/data/Molecular/Masumi/Bioinf.LaRomaine/catchment_microbial_community")
#----------#
# PACKAGES #
#----------#
library(phyloseq)
library(metagenomeSeq)
library(tidyverse)
library(data.table)
library(doMC)
library(plyr)
# allow progress bar
library(svMisc)
library(pbapply)

#-----------#
# FUNCTIONS #
#-----------#
source("./Functions/phyloseq_index.R")
#source("./Functions/peak_identification.R")

#---------------------#
# Raw data processing #
#---------------------#
# Read in microbial data
tax <-
  readRDS("./MotherData/taxtab_idtaxa_silva_v138_2018.rds") # assigned taxonomy
seqtab <- readRDS("./MotherData/nochim_seqtab_2018.rds")

# Read in meta data
meta <-
  read.csv(
    "./MotherData/master_bacmeta.csv",
    sep = ";",
    dec = ".",
    stringsAsFactors = F
  )

# phyloseq needs the sample names of the meta data to be the same as the microbial data
meta <- sample_data(meta)

# Assign rownames to be Sample ID's
rownames(meta) <- meta$DadaNames

# Construct phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = F),
               sample_data(meta),
               tax_table(tax))

# Filter only bacteria, omitting chloroplasts and mitochondria
pb <- ps %>%
  subset_taxa(domain == "Bacteria" &
                family  != "Mitochondria" &
                order   != "Chloroplast")

#----------------------------#
# Rename ASV exact sequences #
#----------------------------#
# Rename exact sequences to ASV_n
asv_id <- data.frame(
  Sequence = taxa_names(pb),
  ID = paste("ASV", seq(length = length(taxa_names(
    pb
  ))), sep = "_"),
  stringsAsFactors = F
)

# create folder for saving
dir.create(file.path("./Output"))
# save for later cross checking
write.table(
  asv_id,
  paste0("./Output/ASV_sequences_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = F
)

taxa_names(pb)[1:5] # check ASVs = exact sequence as colname
# actually overwrite
taxa_names(pb) <-
  paste("ASV", seq(length = length(taxa_names(pb))), sep = "_")
taxa_names(pb)[1:5]

# Filter out ASVs that do not have any abundance in data set
pb <- filter_taxa(pb, function(x)
  sum(x) > 0, TRUE)

#--------------------------------#
# Shrink data to relevant subset #
#--------------------------------#
# Add another column of library size to meta data
sample_data(pb)$LibrarySize <- sample_sums(pb)

# Subset to data set that focuses on PR -> RO2
# Subset only shallow samples
pb <- subset_samples(pb, SeqDepth == "Shallow")

# Subset only 2015 to 2017
pb <- subset_samples(pb, Year == 2015 | Year == 2016 | Year == 2017)

# Omit balnk and Bioassay
pb <- subset_samples(pb, !(sample.type.year == "Bioassay" | sample.type.year == "Blank"))

# remove ASVs that do not appear in this dataset
pb <- prune_taxa(taxa_sums(pb) != 0, pb)

# export shrinked ASV table
asv.tab <- otu_table(pb)
write.table(
  asv.tab,
  paste0("./Output/201520162017_raw_asv_table_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = F
)

# export shrinked data
met.df <- sample_df(pb)
write.table(
  met.df,
  paste0("./Output/201520162017_meta_data_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = F
)

tax.df <- tax_mat(pb)
write.table(
  tax.df,
  paste0("./Output/201520162017_tax_table_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = T
)

# remove unnecessary objects
rm(asv_id, meta, ps, seqtab, tax)

###########################
#- Correct too few reads -#
###########################
# transform ASV table into data.table
asv.tab <- setDT(as.data.frame(asv.tab, strings.As.Factors = F), keep.rownames = "Sample")

# melt asv table into long format
asv.tab <- melt.data.table(
  asv.tab,
  id.vars = "Sample",
  measure.vars = patterns("^ASV_"),
  variable.name = "ASV",
  value.name = "reads"
)

# Join asv and meta table
sumdf <- left_join(
  asv.tab,
  met.df[met.df$DadaNames %in% asv.tab$Sample,] %>%
    dplyr::select(
      DadaNames,
      Year,
      Season,
      DnaType,
      sample.type,
      sample.type.year,
      soilorwater,
      catchment.area,
      distance.from.mouth
    ),
  by = c("Sample" = "DadaNames")
)

# set back to data.table, order data.table by catchment.area
setDT(sumdf)

# adjust catchment area to include soil, soilwater, and hyporheic water
#sumdf[, catchment.area := catchment.area + 45] # add 45 to all

sumdf[sample.type.year == "Soil", catchment.area := -10]
sumdf[sample.type.year == "Soilwater", catchment.area := -20]
sumdf[sample.type.year == "Hyporheicwater", catchment.area := -30]

# add ID column for parallel computing
sumdf[, ID := paste(Year, Season, DnaType, ASV, sep = ".")]

# order by catchment.area
sumdf <- sumdf[order(ID, catchment.area)]

# split data frame in present and absent
# (Otherwise computational power is overwhelmed)
sumdf[, n := .N, by = .(ID)] # number of samples by factorial combination
sumdf[, n.obs := nrow(.SD[reads > 0]), by = .(ID)] # how many of those have an actual observation of ASV?

# initiate presence-absence col
sumdf[, PA := character()]
sumdf[is.na(PA) &
        n.obs == 0, PA := "Absent", by = .(ID)] # if for a factorial combination, there was not a single observation define as absent
sumdf[is.na(PA) & n.obs > 0, PA := "Present", by = .(ID)]
# all observations above 0 are present

#- Actual quality control -#

# register number of cores
detectCores() #24, we do not have 24. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12
# check if cores were allocated correctly

# create folder for saving
dir.create(file.path("./Objects"))

present <- sumdf[PA == "Present",]
absent <- sumdf[PA == "Absent",]

# remove single observations with read numbers lower than 10 reads by sample.type ~ Year ~ Season ~ DnaType combination
present[, bin := paste(sample.type.year, Year, Season, DnaType, ASV, sep = "_")]
# number of observations by bin
present[, n.bin := .N, by = .(bin)]
present[n.bin == 1 & reads < 10, cor.reads := 0] # select all with only one observation by bin, overwrite all reads less than 10 with 0
present[is.na(cor.reads), cor.reads := reads] # fill in rest
present[, bin := NULL][, n.bin := NULL] # remove unncessary columns for downstream analysis

# add 0 column to absent for cor.reads
absent[, cor.reads := 0]

# combine present and absent back together
fin <- bind_rows(present, absent)

# extract corrected reads for metagenomeseq as matrix
cor.reads <- setDF(spread(fin %>% select(Sample, ASV, cor.reads), key = ASV, value = cor.reads))
row.names(cor.reads) <- cor.reads$Sample
cor.reads$Sample <- NULL
cor.reads <- as.matrix(cor.reads)

# join back to phyloseq so that orders of samples and ASVs match across data frames
cor.pb <- phyloseq(otu_table(cor.reads, taxa_are_rows = F),
                   sample_data(met.df),
                   tax_table(tax.df))

######################
# CSS transformation #
######################
# We apply the cumuluative sum scaling transformation (Paulson et al. Nature Methods 2013)
# metagenomeSeq needs, samples in columns, "features" (= ASVs) in rows
# export phyloseq object into MRexperiment
pb.ms <- phyloseq_to_metagenomeSeq_mod(cor.pb) # encountered error with original function use customised function

# inspect data
head(pData(pb.ms), 3) # meta data
head(fData(pb.ms), 3) # taxonomy table
head(MRcounts(pb.ms[, 1:2])) # count table

# we've done all the subsetting in phyloseq, directly move to CSS transformation

# calculate normalisation factor
p <- round(cumNormStatFast(pb.ms), digits = 2)

# apply CSS transformation
pb.ms <- cumNorm(pb.ms, p = p)

# check if normalisation was succesful
t(otu_table(cor.pb))[1:5, 1:5]
MRcounts(pb.ms, norm = T)[1:5, 1:5]

# export normalised count matrices
pb.mat <- MRcounts(pb.ms, norm = T)
exportMat(pb.mat, file = paste0("./Output/201520162017_CSS_asvtab_", Sys.Date(),".tsv"))

# export sample statistics
# create folder for saving
dir.create(file.path("./Output/StatTables"))
exportStats(pb.ms, file = paste0("./Output/StatTables/201520162017_CSS_transf_stats_",Sys.Date(),".tsv"))
# read with read.csv(file, sep = "\t")

# transform data into long format to merge with fin
# melt asv table into long format
# to do this easily we transpose and melt matrix
css <- melt.data.table(
  setDT(as.data.frame(t(pb.mat)), keep.rownames = "Sample"),
  id.vars = "Sample",
  measure.vars = patterns("^ASV_"),
  variable.name = "ASV",
  value.name = "css.reads"
)

fin <- merge(fin, css, by = c("Sample", "ASV"))

#---------------------#
# Relative abundances #
#---------------------#
# !! not recommended to use !! #
# but calculate for completeness
# calculate relative abundances
fin[, library.size := sum(css.reads), .(Sample)]
fin[, rel.abund := css.reads / library.size]
# 273 samples

#sanity check
san <- fin[, .(check = sum(rel.abund)), .(Sample)]
san$check # all 1

#---------------------#
# Who became absent ? #
#---------------------#
# Classify those as absent after correcting by bin
fin[, n := .N, by = .(ID)] # number of samples by factorial combination
fin[, n.obs := nrow(.SD[css.reads > 0]), by = .(ID)] # how many of those have an actual observation of ASV?

# add z-standardised css.reads to compare rare and abundant things later in regressions
fin[, z.css.reads := (css.reads - mean(css.reads, na.rm = T)) / sd(css.reads, na.rm = T)]

#------#
# Save #
#------#
# rearrange data frame before saving
setcolorder(
  fin,
  c(
    "ID",
    "Sample",
    "ASV",
    "Year",
    "Season",
    "DnaType",
    "sample.type",
    "sample.type.year",
    "soilorwater",
    "catchment.area",
    "distance.from.mouth",
    "library.size",
    "n",
    "n.obs",
    "PA",
    "reads",
    "cor.reads",
    "css.reads",
    "rel.abund",
    "z.css.reads"
  )
)

saveRDS(fin,
        paste0("./Objects/201520162017_css_",
               Sys.Date(),
               ".rds"))

##################################################################
###################################################################################################
## Rarefaction
# overwrite library size data with corrected reads
sample_data(cor.pb)$LibrarySize <- sample_sums(cor.pb)
hist(sample_data(cor.pb)$LibrarySize)

# how many samples do we loose with different thresholds?
higher <- subset_samples(cor.pb, LibrarySize < 70000)
high <- subset_samples(cor.pb, LibrarySize < 60000)
med2 <- subset_samples(cor.pb, LibrarySize < 30000)
med <- subset_samples(cor.pb, LibrarySize < 25000)
low <- subset_samples(cor.pb, LibrarySize < 20000)
lower <- subset_samples(cor.pb, LibrarySize < 15000)
sample_data(lower)

data.frame(
  LibrarySize = c(70000, 60000, 30000, 25000, 20000, 15000),
  Count = c(
    length(sample_names(higher)),
    length(sample_names(high)),
    length(sample_names(med2)),
    length(sample_names(med)),
    length(sample_names(low)),
    length(sample_names(lower))
  )
)

# Almost all are sample with low library size are soily samples.
# How many samples do we still have for each factor combination if we remove them?
sample_data(subset_samples(cor.pb, LibrarySize >= 25000)) %>%
  group_by(year, Season, DnaType, soilorwater) %>% dplyr::summarise(n())
# the sample sizes are acceptable for now. We'll go for rarefaction at 25000.

# Because rarefaction is random and we might be introducing biases,
# we will do a resampling technique to estimate the mean of each ASV in the samples.

# make a loop to try different library sizes
trials <- 100  #set number of permutation
min_lib <- c(15000, 25000, 50000) # set minimum library size

permlist <- list() #empty list to fill by loop

# Loop 100 times and rarefy
for (j in 1:length(min_lib)) {
  for (i in 1:trials) {
    # subsample e.g. rarefy
    set.seed(i)
    permlist[[i]] <-
      data.table(
        otu_mat(
          rarefy_even_depth(
            pb,
            sample.size = min_lib[j],
            verbose = F,
            replace = T
          )
        ),
        keep.rownames = T,
        stringsAsFactors = F
      )
    
    progress(i)
    Sys.sleep(0.01)
    if (i == trials)
      cat("Done!\n")
  }
  
  # We have completed the iteration. Now, we will recombine the iteration
  # results into one data frame and calculate the mean and standard deviation for each ASV in each sample.
  permdf <-
    bind_rows(permlist)  # combine all data frames within lists into one data frame
  
  means <-
    permdf[, pblapply(.SD, mean, na.rm = T), by = list(rn)]  #calculate means by row and group
  sd <-
    permdf[, pblapply(.SD, sd, na.rm = T), by = list(rn)]  #calculate standard deviation by row and group
  
  # melt dataframes
  melmean <- melt.data.table(
    means,
    id.vars = c("rn"),
    measure.vars = patterns("^ASV_"),
    variable.name = "ASV",
    value.name = "iter.mean"
  )
  melsd <- melt.data.table(
    sd,
    id.vars = c("rn"),
    measure.vars = patterns("^ASV_"),
    variable.name = "ASV",
    value.name = "iter.sd"
  )
  
  # merge
  perm.rar <- merge(melmean, melsd, by = c("rn", "ASV"), all = T)
  apply(perm.rar, 2, function(x)
    any(is.na(x))) # any NAs? - yes in SD
  colnames(perm.rar)[1] <- "Sample" # rename first ID column
  
  # save
  write.table(
    perm.rar,
    paste0("./Output/perm.rar_lib", min_lib[j], "_", Sys.Date(), ".csv"),
    sep = ";",
    dec = ".",
    row.names = F
  )
}


#---------------------#
#------- Done! -------#
# Move to next script #
#---------------------#
sessionInfo()

