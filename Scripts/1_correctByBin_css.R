#-- Script for the publication:
#-- Title
#-- Responsible author: Masumi Stadler

# This script is the second of a series of scripts that were used to analyse the data
# used in the publication.

###-----------------------------------###
#-   Binning and CSS transformation   - #
###-----------------------------------###
# This script processes the meta data of the whole La Romaine project to extract the necessary information for the paper

# First we take all samples except those from Bioassays and 2018 (as RNA data is not avaialable at this moment)
# For users/viewers: Only parts of the whole data will be published at point of publication
# as other data will be used for further thesis chapters, thus a few lines will have to be skipped as they are only
# relevant for the whole data

# Second, we split the data set into ID groups (Year, Season, DnaType, ASV) and evaluate whether
# there is an actual observation or if that particular ASV is absent in that ID group.
# To do this, we will split the x-axis which is catchment.area into bins. If within a bin,
# there are several samples but only one actual observation then we will remove the observation
# if the number of reads is below 10. This function consumes a lot of memory and CPU power,
# thus we only execute the bin_qual_control() function on the actually present IDs.

# Third, we will transform the corrected counts with metgenomeSeq to CSS ().

# Finally, we calculate the relative abundances based on the css read numbers for completeness.

# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("phyloseq", "metagenomeSeq", # OTU clustering
              "tidyverse", "plyr", "data.table", # wrangling
              "doMC", "pbapply", # parallel computing
              "svMisc") # progress bar in loop


### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

### Functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")

### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# Set seed for session and reproducibility of permutations
# (just for consistency, no random iteration in this script)
set.seed(3)

# 2. Raw data processing ---------------------------------------------------------------
# Read in OTU table and taxonomy information
#tax <- readRDS("./Objects/otu_taxtab_99.rds")
#seqtab <- readRDS("./Objects/otu_seqtab_99.rds")
tax <- read.csv("./Output/OTU_98_taxonomy.csv", 
                sep = ',', row.names = 1, colClasses = "character", stringsAsFactors = F)
tax.names <- rownames(tax) # somehow mutate_if removes rownames, save and add later
tax <- tax %>% mutate_if(is.character, list(~na_if(.,"")))
# some NAs are just an empty string, fill with NA
rownames(tax) <- tax.names; rm(tax.names)
tax <- as.matrix(tax[order(rownames(tax)),])
#seqtab <- read.csv("./Output/OTU_99_table.csv",
#                   sep = ',', row.names = 1, stringsAsFactors = F)
# faster to read in R object
seqtab <- readRDS('./Objects/OTU_98_table.rds')
seqtab <- seqtab[order(rownames(seqtab)),order(colnames(seqtab))]


# Read in meta data
meta <-
  read.csv(
    "./MotherData/master_bacmeta.csv",
    sep = ";",
    dec = ".",
    stringsAsFactors = F
  )

# we have more samples than in the data set needed (aka meta), as new 2018 data were added.
# Omit rows with 2018 samples
seqtab <- seqtab[rownames(seqtab) %in% meta$DadaNames,]

# phyloseq needs the sample names of the meta data to be the same as the microbial data
meta <- sample_data(meta)

# Assign rownames to be Sample ID's
rownames(meta) <- meta$DadaNames
meta <- meta[order(rownames(meta)),]

# Construct phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = F),
               sample_data(meta),
               tax_table(tax))

# How many unclassified OTUs?
t <- ps %>% subset_taxa(taxa_sums(ps) != 0) %>% subset_taxa(is.na(domain))
nrow(tax_mat(t))
rm(t)

# Filter only bacteria, omitting chloroplasts and mitochondria
pb <- ps %>%
  subset_taxa(domain == "Bacteria" &
                family  != "Mitochondria" &
                order   != "Chloroplast")

# 3. Shrink to relevant subset ---------------------------------------------------------------
# Filter out OTUs that do not have any abundance in data set
pb <- filter_taxa(pb, function(x)
  sum(x) > 0, TRUE)

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

# extract individual tables from phyloseq obj
otu.tab <- otu_table(pb)

met.df <- sample_df(pb)

tax.df <- tax_mat(pb)

# remove unnecessary objects
rm(meta, ps, seqtab, tax)

# 4. Correct too few reads ---------------------------------------------------------------
# transform ASV table into data.table
otu.tab <- setDT(as.data.frame(otu.tab, strings.As.Factors = F), keep.rownames = "Sample")

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

# correct one miscategorisation
sumdf[sumdf$Sample == "RO2111.60mD",]$sample.type.year <- "Deep"


# 5. Match DNA and RNA ---------------------------------------------------------------
# new sample name column to identify which DNAs belong to which RNA
sumdf[, DR.names := Sample]
# correct a few wrong sample names for matching DNA and RNA
sumdf[sumdf$DR.names == "RO2R52R", "DR.names"] <- "RO2.52R"
sumdf[sumdf$DR.names == "SWR34R", "DR.names"] <- "SW34R"
sumdf[sumdf$DR.names == "RO2.36pD", "DR.names"] <- "RO2.36D"
sumdf[sumdf$DR.names == "RO2.36pR", "DR.names"] <- "RO2.36R"
sumdf[sumdf$DR.names == "RO2111.60mD", "DR.names"] <- "RO2111.90mD"
sumdf[sumdf$DR.names == "RO2.30DPR", "DR.names"] <- "RO2.30R" # two DNA
sumdf[sumdf$DR.names == "RO301.HypoR", "DR.names"] <- "RO31.HypoR"
sumdf[sumdf$DR.names == "RO301R", "DR.names"] <- "RO31R" 
sumdf[sumdf$DR.names == "RO304R", "DR.names"] <- "RO34R" 
sumdf[sumdf$DR.names == "RO307R", "DR.names"] <- "RO37R" 
sumdf[sumdf$DR.names == "L230R", "DR.names"] <- "L330R" # L230 does not exist

# remove Ds and Rs to match counterpart DR.names
sumdf$DR.names[sumdf$DnaType == "DNA"] <- str_replace(sumdf$DR.names[sumdf$DnaType == "DNA"], "D$", "")
sumdf$DR.names[sumdf$DnaType == "cDNA"] <- str_replace(sumdf$DR.names[sumdf$DnaType == "cDNA"], "R$", "")

# calculate the sum of OTUs per DnaType, omit those OTUs that only appear in RNA
sum.reads <- sumdf[, .(sum.reads = sum(reads)), by = .(DnaType, OTU)]
notindna <- sum.reads[DnaType == "DNA" & sum.reads == 0,]$OTU

# order for overview
sum.reads <- sum.reads[order(OTU, DnaType),]

nrow(tax.df[rownames(tax.df) %in% sum.reads[OTU %in% as.character(notindna) & sum.reads > 100,]$OTU,])
# 162 OTUs only in RNA and that have a high read number (> 100)
# 490 with ASVs and clustering by 99.7% similarity

sumdf <- sumdf[!(OTU %in% notindna),]
length(notindna)
# removing 3353 OTUs

#combine back with some meta Data and sample names
# add ID column for parallel computing
sumdf[, ID := paste(Year, Season, DnaType, OTU, sep = ".")]

# split data frame in present and absent
# (Otherwise computational power is overwhelmed)
sumdf[, n := .N, by = .(ID)] # number of samples by factorial combination
sumdf[, n.obs := nrow(.SD[reads > 0]), by = .(ID)] # how many of those have an actual observation of an OTU?

# initiate presence-absence col
sumdf[, PA := character()]
sumdf[is.na(PA) &
        n.obs == 0, PA := "Absent", by = .(ID)]
# if for a factorial combination, there was not a single observation define as absent
sumdf[is.na(PA) & n.obs > 0, PA := "Present", by = .(ID)]
# all observations above 0 are present

# 6. Quality control ---------------------------------------------------------------
getDoParWorkers() # 12
# check if cores were allocated correctly

# create folder for saving
dir.create(file.path("./Objects"))

present <- sumdf[PA == "Present",]
absent <- sumdf[PA == "Absent",]

# remove single observations with read numbers lower than 10 reads by sample.type ~ Year ~ Season ~ DnaType combination
present[, bin := paste(sample.type.year, Year, Season, DnaType, OTU, sep = "_")]
# number of observations by bin
present[, n.bin := .N, by = .(bin)]
present[n.bin == 1 & reads < 10, cor.reads := 0]
# select all with only one observation by bin, overwrite all reads less than 10 with 0
present[is.na(cor.reads), cor.reads := reads] # fill in rest
present[, bin := NULL][, n.bin := NULL] # remove unncessary columns for downstream analysis

# add 0 column to absent for cor.reads
absent[, cor.reads := 0]

# combine present and absent back together
fin <- bind_rows(present, absent)

# extract corrected reads for metagenomeseq as matrix
cor.reads <- setDF(spread(fin %>% select(Sample, OTU, cor.reads), key = OTU, value = cor.reads))
row.names(cor.reads) <- cor.reads$Sample
cor.reads$Sample <- NULL
cor.reads <- as.matrix(cor.reads)

# one sample has no OTU reads left
cor.reads <- cor.reads[which(rowSums(cor.reads) != 0),]

# join back to phyloseq so that orders of samples and ASVs match across data frames
cor.pb <- phyloseq(otu_table(cor.reads, taxa_are_rows = F),
                   sample_data(met.df),
                   tax_table(tax.df[row.names(tax.df) %in% colnames(cor.reads),]))

# save initial sumdf data frame to be used later as reference (sample ~ corrected sample name)
saveRDS(sumdf, "./Objects/summary.meta.with.oldnames.rds")

# 7. CSS transformation ---------------------------------------------------------------
# We apply the cumuluative sum scaling transformation (Paulson et al. Nature Methods 2013)
# metagenomeSeq needs, samples in columns, "features" (= ASVs) in rows
# export phyloseq object into MRexperiment
pb.ms <- physeq_to_metagenomeSeq_mod(cor.pb) # encountered error with original function use customised function

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
t(otu_table(cor.pb))[1100:1110, 1:5]
MRcounts(pb.ms, norm = T)[1100:1110, 1:5]

# export normalised count matrices
pb.mat <- MRcounts(pb.ms, norm = T)
pb.mat <- round(pb.mat, digits = 0) # export as count data instead of decimals
exportMat(pb.mat, file = paste0("./Output/201520162017_CSS_otu99_otutab_", Sys.Date(),".tsv"))

# export sample statistics
# create folder for saving
dir.create(file.path("./Output/StatTables"))
exportStats(pb.ms, file = paste0("./Output/StatTables/201520162017_CSS_otu99_transf_stats_",Sys.Date(),".tsv"))
# read with read.csv(file, sep = "\t")

# transform data into long format to merge with fin
# melt asv table into long format
# to do this easily we transpose and melt matrix
css <- melt.data.table(
  setDT(as.data.frame(t(pb.mat)), keep.rownames = "Sample"),
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "css.reads"
)

fin <- merge(fin, css, by = c("Sample", "OTU"))

# 8. Calculate relative abundances ---------------------------------------------------------------
# !! not recommended to use !! #
# but calculate for completeness
# calculate relative abundances
fin[, library.size := sum(css.reads), .(Sample)]
fin[, rel.abund := css.reads / library.size]
# 273 samples

#sanity check
san <- fin[, .(check = sum(rel.abund)), .(Sample)]
san$check # all 1

# 9. Identify who became absent after quality control ---------------------------------------------------------------
# Classify those as absent after correcting by bin
fin[, n := .N, by = .(ID)] # number of samples by factorial combination
fin[, n.obs := nrow(.SD[css.reads > 0]), by = .(ID)] # how many of those have an actual observation of ASV?

# add z-standardised css.reads to compare rare and abundant things later in regressions
fin[, z.css.reads := (css.reads - mean(css.reads, na.rm = T)) / sd(css.reads, na.rm = T)]

# calculate mean reads for duplicates
dupl.mean <- fin[, .(reads = mean(reads, na.rm = T),
                  cor.reads = mean(cor.reads, na.rm = T),
                  css.reads = mean(css.reads, na.rm = T), 
                  rel.abund = mean(rel.abund, na.rm = T),
                  z.css.reads = mean(z.css.reads, na.rm = T)), by = .(DR.names, DnaType, OTU)]
out <- dupl.mean[fin, c("Sample","Season","Year","sample.type.year","soilorwater","catchment.area",
              "distance.from.mouth", "n", "n.obs", "PA", "ID", "library.size") := 
            list(i.Sample,i.Season,i.Year,i.sample.type.year,i.soilorwater,i.catchment.area,
                 i.distance.from.mouth, i.n, i.n.obs, i.PA, i.ID, i.library.size), on = .(DR.names, DnaType, OTU)]

# 11. Export ---------------------------------------------------------------
# rearrange data frame before saving
setcolorder(
  out,
  c(
    "ID",
    "Sample",
    "DR.names",
    "OTU",
    "Year",
    "Season",
    "DnaType",
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

saveRDS(out,
        paste0("./Objects/201520162017_css_otu99_",
               Sys.Date(),
               ".rds"))
write.table(out, paste0("./Output/201520162017_css_otu99_",
                      Sys.Date(),
                      ".csv"), sep = ",", dec = ".",  row.names = F)
## Extract final info from all phyloseq objects
# export shrinked ASV table

otu.tab <- setDF(dcast(out, Sample ~ OTU, value.var = "css.reads"))
row.names(otu.tab) <- otu.tab$Sample; otu.tab[, "Sample"] <- NULL
otu.tab <- as.matrix(otu.tab)

# remove empty samples or empty OTUs
otu.tab <- otu.tab[rowSums(otu.tab) > 0,]
otu.tab <- otu.tab[, colSums(otu.tab) > 0]

# row orders need to match between tax.tab and otu.tab
otu.tab <- otu.tab[order(row.names(otu.tab)),]
write.table(
  otu.tab,
  paste0("./Output/201520162017_fin_css_otu99_table_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = T
)

# export shrinked meta data
met.df <- sample_df(pb)

# only take those samples in OTU table, will remove duplicates and one lost sample
met.df <- met.df[met.df$DadaNames %in% levels(factor(row.names(otu.tab))),]
setDT(met.df)
met.df <- met.df[out, c("DR.names", "sample.type.year") := 
                 list(i.DR.names, i.sample.type.year), on = .(DadaNames == Sample)]
setDF(met.df)
# Assign rownames to be Sample ID's
rownames(met.df) <- met.df$DadaNames
write.table(
  met.df,
  paste0("./Output/201520162017_meta_otu99_data_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = T
)

tax.df <- tax_mat(pb)
tax.df <- tax.df[row.names(tax.df) %in% levels(factor(colnames(otu.tab))),]
# orders need to match between tax.tab and otu.tab
tax.df <- tax.df[order(row.names(tax.df)),]

write.table(
  tax.df,
  paste0("./Output/201520162017_tax_otu99_table_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = T
)

#rm(met.df, out, pb, pb.mat, pb.ms, present, san, seqdf, sum.reads, sumdf, tax.df, tax.tab, test,notindna, p, absent,
#   otu.tab, cor.reads, css, dupl.mean, fin)




# 12. Rarefaction ---------------------------------------------------------------
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
            cor.pb,
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
    measure.vars = patterns("^OTU_"),
    variable.name = "OTU",
    value.name = "iter.mean"
  )
  melsd <- melt.data.table(
    sd,
    id.vars = c("rn"),
    measure.vars = patterns("^OTU_"),
    variable.name = "OTU",
    value.name = "iter.sd"
  )
  
  # merge
  perm.rar <- merge(melmean, melsd, by = c("rn", "OTU"), all = T)
  apply(perm.rar, 2, function(x)
    any(is.na(x))) # any NAs? - yes in SD
  colnames(perm.rar)[1] <- "Sample" # rename first ID column
  
  # save
  write.table(
    perm.rar,
    paste0("./Output/perm.rar_lib_otu99_", min_lib[j], "_", Sys.Date(), ".csv"),
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


#R version 4.0.3 (2020-10-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.5 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

#locale:
#  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
#[3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
#[5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
#[7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods  
#[9] base     

#other attached packages:
#  [1] doMC_1.3.7          iterators_1.0.13    foreach_1.5.1       data.table_1.13.2  
#[5] plyr_1.8.6          forcats_0.5.0       stringr_1.4.0       dplyr_1.0.2        
#[9] purrr_0.3.4         readr_1.4.0         tidyr_1.1.2         tibble_3.0.4       
#[13] ggplot2_3.3.2       tidyverse_1.3.0     DECIPHER_2.16.1     RSQLite_2.2.1      
#[17] Biostrings_2.56.0   XVector_0.28.0      IRanges_2.22.2      S4Vectors_0.26.1   
#[21] BiocGenerics_0.34.0

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.5       lubridate_1.7.9  assertthat_0.2.1 digest_0.6.27   
#[5] R6_2.4.1         cellranger_1.1.0 backports_1.1.10 reprex_0.3.0    
#[9] httr_1.4.2       pillar_1.4.6     zlibbioc_1.34.0  rlang_0.4.8     
#[13] readxl_1.3.1     rstudioapi_0.11  blob_1.2.1       bit_4.0.4       
#[17] munsell_0.5.0    broom_0.7.2      compiler_4.0.3   modelr_0.1.8    
#[21] pkgconfig_2.0.3  tidyselect_1.1.0 codetools_0.2-17 fansi_0.4.1     
#[25] crayon_1.3.4     dbplyr_1.4.4     withr_2.3.0      grid_4.0.3      
#[29] jsonlite_1.7.1   gtable_0.3.0     lifecycle_0.2.0  DBI_1.1.0       
#[33] magrittr_1.5     scales_1.1.1     cli_2.1.0        stringi_1.5.3   
#[37] fs_1.5.0         xml2_1.3.2       ellipsis_0.3.1   generics_0.0.2  
#[41] vctrs_0.3.4      tools_4.0.3      bit64_4.0.5      glue_1.4.2      
#[45] hms_0.5.3        colorspace_1.4-1 rvest_0.3.6      memoise_1.1.0   
#[49] haven_2.3.1