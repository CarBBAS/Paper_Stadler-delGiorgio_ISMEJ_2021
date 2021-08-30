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
# If within a bin, there are several samples but only one actual observation then
# we will remove the observation if the number of reads is below 10

# Third, we will transform the corrected counts with metgenomeSeq to CSS reads.

# Finally, we calculate the relative abundances based on the CSS read numbers for completeness.

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
tax <- readRDS("./Objects/OTU_99_gtdb_paper1_taxonomy.rds")
seqtab <- readRDS("./Objects/OTU_99_gtdb_paper1_table.rds")
#tax <- read.csv("./Output/OTU_99_gtdb_taxonomy.csv", 
#                sep = ',', row.names = 1, colClasses = "character", stringsAsFactors = F)
#seqtab <- read.csv("./Output/OTU_99_gtdb_table.csv",
#                   sep = ',', row.names = 1, stringsAsFactors = F)

tax.names <- rownames(tax) # somehow mutate_if removes rownames, save and add later
tax <- tax %>% as.data.frame() %>% mutate_if(is.character, list(~na_if(.,"")))
# some NAs are just an empty string, fill with NA
rownames(tax) <- tax.names; rm(tax.names)
tax <- as.matrix(tax[order(rownames(tax)),])

seqtab <- seqtab[order(rownames(seqtab)),order(colnames(seqtab))]

#seqtab[row.names(seqtab) == "HW24R", "OTU_128"]

# # Read in meta data
# meta <-
#   read.csv(
#     "./MotherData/main_bac_match_2021-02-22.csv", #/master_bacmeta
#     sep = ",",
#     dec = ".",
#     stringsAsFactors = F
#   )
# 
# naming <-  read.csv(
#   "./MotherData/sequence_metadata.csv", #/master_bacmeta
#   sep = ",",
#   dec = ".",
#   stringsAsFactors = F
# )
# 
# # merge the two
# setDT(meta); setDT(naming)
# 
# meta <- naming[meta[,-c("year","lat","long","sampling.date","sample.type.year"), with = T], ,
#                on = c("dr_match_name" = "dna.match")]
# rm(naming)
# 
# # extract columns used in this paper
# meta <- meta %>% dplyr::select(seq_name:sample.type,Season,depth,distance.from.mouth, catchment.area,Comments)
# write.table(meta, "./MotherData/meta_paper1.csv", sep = ",", dec = ".", row.names = F)

meta <- read.csv("./MotherData/meta_paper1.csv", sep = ",", dec = ".", stringsAsFactors = F)

# phyloseq needs the sample names of the meta data to be the same as the microbial data
meta <- sample_data(meta)

# Assign rownames to be Sample ID's
rownames(meta) <- meta$seq_name
meta <- meta[order(rownames(meta)),]

# Construct phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = F),
               sample_data(meta),
               tax_table(tax))

# How many unclassified OTUs?
t <- ps %>% subset_taxa(taxa_sums(ps) != 0) %>% subset_taxa(is.na(domain))
nrow(tax_mat(t)) # 35995
# of a total
nrow(tax_mat(ps)) # 81266
rm(t)

# Filter only bacteria, omitting chloroplasts and mitochondria
pb <- ps %>%
  subset_taxa(domain == "Bacteria" &
                family  != "Mitochondria" &
                order   != "Chloroplast")

# 3. Separate duplicates ------------------------------------------------------------------------
# Filter out OTUs that do not have any abundance in data set
pb <- filter_taxa(pb, function(x)
  sum(x) > 0, TRUE)

# Add another column of library size to meta data
sample_data(pb)$LibrarySize <- sample_sums(pb)

dups <- sample_df(pb)
# extract duplicates
setDT(dups)
t <- ddply(dups, .(dr_match_name, seq_depth, dna_type), function(x){
  if(nrow(x) > 1){
    return(x)
  }
}) %>% select(dr_match_name, seq_depth, dna_type, seq_name, replicate, LibrarySize)

which.bigger <-ddply(t, .(dr_match_name, dna_type, seq_depth), function(x){
  x[which.max(x$LibrarySize),]
})
# seems like mostly replicate 2 has a larger library size. But this is sample dependent.
# extract always the sample with more replicates
# BUT we keep both duplicates for "true" replicates (labelled t2)

# first, make a df with all samples (we will do the PCoA with both datasets)
all <- pb

# extract only one sample per duplicate
t <- t[t$dr_match_name %in% t[grep("s2", t$seq_name),]$dr_match_name,]

which.bigger <-ddply(t, .(dr_match_name, dna_type, seq_depth), function(x){
  x[which.max(x$LibrarySize),]
})

# samples that have s2 duplicates
dupl.names <- t$seq_name

# samples that have higher library size
better.samples <- which.bigger$seq_name

# extract a df without any s2 samples and add "better" duplicates
no.dups <- rbind(dups[!(seq_name %in% dupl.names),],
      dups[seq_name %in% better.samples,])

# sanity check
nrow(no.dups) + length(better.samples) == nrow(sample_df(all))

# extract a phyloseq object
no.dups <- subset_samples(all, seq_name %in% no.dups$seq_name)
# remove taxa that don't appear in this data set
no.dups <- filter_taxa(no.dups, function(x)
  sum(x) > 0, TRUE)

# 4. Shrink to relevant subset -------------------------------------------------------------------------------------
# Subset only shallow samples
#all <- subset_samples(all, seq_depth == "Shallow")
no.dups <- subset_samples(no.dups, seq_depth == "Shallow")

# Subset only 2015 to 2017
#all <- subset_samples(all, year == 2015 | year == 2016 | year == 2017)
no.dups <- subset_samples(no.dups, year == 2015 | year == 2016 | year == 2017)

# Omit blank and Bioassay
#all <- subset_samples(all, !(sample.type.year == "Bioassay" | sample.type.year == "Blank"))
no.dups <- subset_samples(no.dups, !(sample.type.year == "Bioassay" | sample.type.year == "Blank"))

# remove ASVs that do not appear in this dataset
#all <- prune_taxa(taxa_sums(all) != 0, all)
no.dups <- prune_taxa(taxa_sums(no.dups) != 0, no.dups)

# Samples without (technical) duplicates
# run again with all, if desired

# extract individual tables from phyloseq obj
otu.tab <- otu_table(no.dups)

met.df <- sample_df(no.dups)

tax.df <- tax_mat(no.dups)

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
  met.df[met.df$seq_name %in% otu.tab$Sample,] %>%
    dplyr::select(
      seq_name,
      dr_match_name,
      year,
      Season,
      dna_type,
      sample.type,
      sample.type.year,
      catchment.area,
      distance.from.mouth
    ),
  by = c("Sample" = "seq_name")
)

# set back to data.table, order data.table by catchment.area
setDT(sumdf)


# adjust catchment area to include soil, soilwater, and hyporheic water
#sumdf[, catchment.area := catchment.area + 45] # add 45 to all

sumdf[sample.type.year == "Soil", catchment.area := -10]
sumdf[sample.type.year == "Soilwater", catchment.area := -20]
sumdf[sample.type.year == "Hyporheicwater", catchment.area := -30]

# calculate the sum of OTUs per DnaType, omit those OTUs that only appear in RNA
sum.reads <- sumdf[, .(sum.reads = sum(reads)), by = .(dna_type, OTU)]
notindna <- sum.reads[dna_type == "DNA" & sum.reads == 0,]$OTU

# order for overview
sum.reads <- sum.reads[order(OTU, dna_type),]

nrow(tax.df[rownames(tax.df) %in% sum.reads[OTU %in% as.character(notindna) & sum.reads > 100,]$OTU,])
# 147 OTUs only in RNA and that have a high read number (> 100)

sumdf <- sumdf[!(OTU %in% notindna),]
length(notindna)
# removing 2615 OTUs

length(unique(sumdf$OTU))
# we still have 16445 OTUs of...
length(unique(otu.tab$OTU))
# 19060

#combine back with some meta Data and sample names
# add ID column for parallel computing
sumdf[, ID := paste(year, Season, dna_type, OTU, sep = ".")]

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

# remove single observations with read numbers lower than 10 reads by sample.type ~ Year ~ Season ~ dna_type combination
present[, bin := paste(sample.type.year, year, Season, dna_type, OTU, sep = "_")]
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

# remove any OTU with no reads anymore
cor.reads <- cor.reads[,names(which(colSums(cor.reads) != 0))]
cor.reads <- cor.reads[,names(which(colSums(cor.reads) != 1))] # remove singletons
# one sample has no OTU reads left
cor.reads <- cor.reads[which(rowSums(cor.reads) > 0),]
cor.reads <- cor.reads[which(rowSums(cor.reads) > 1),] # remove singletons

# any samples without any read?
any(rowSums(cor.reads) <= 1)
# any OTU without any read?
any(colSums(cor.reads) <= 1) # we have singletons, problem with CSS = remove

# we have one sample, that has only one observation of an OTU but with many reads = remove
single.count<-apply(cor.reads, 1, function(row){
  zero <- length(row[row <= 1])
  if(length(row) - zero == 1){
    return(TRUE)
  } else {return(FALSE)}
}) # which sample has only one count?
cor.reads <- cor.reads[-which(single.count == T),]

# join back to phyloseq so that orders of samples and ASVs match across data frames
cor.all <- phyloseq(otu_table(cor.reads, taxa_are_rows = F),
                   sample_data(met.df),
                   tax_table(tax.df[row.names(tax.df) %in% colnames(cor.reads),]))

# save initial sumdf data frame to be used later as reference (sample ~ corrected sample name)
#saveRDS(sumdf, "./Objects/summary.meta.with.oldnames.rds")

# 7. CSS transformation ---------------------------------------------------------------
# We apply the cumuluative sum scaling transformation (Paulson et al. Nature Methods 2013)
# metagenomeSeq needs, samples in columns, "features" (= ASVs) in rows
# export phyloseq object into MRexperiment
all.ms <- physeq_to_metagenomeSeq_mod(cor.all) # encountered error with original function use customised function

# inspect data
head(pData(all.ms), 3) # meta data
head(fData(all.ms), 3) # taxonomy table
head(MRcounts(all.ms[, 1:2])) # count table

any(rowSums((MRcounts(all.ms))) <= 1)
any(colSums((MRcounts(all.ms))) <= 1)

# we've done all the subsetting in phyloseq, directly move to CSS transformation

# calculate normalisation factor
p <- round(cumNormStatFast(all.ms), digits = 2)

# apply CSS transformation
all.ms <- cumNorm(all.ms, p = p)

# check if normalisation was succesful
t(otu_table(cor.all))[1100:1110, 1:5]
MRcounts(all.ms, norm = T)[1100:1110, 1:5]

# export normalised count matrices
all.mat <- MRcounts(all.ms, norm = T)
#all.mat <- round(all.mat, digits = 0) # export as count data instead of decimals
exportMat(all.mat, file = paste0("./Output/201520162017_CSS_otu99_otutab_paper1_", Sys.Date(),".tsv"))

# all.mat <- read.csv("./Output/201520162017_CSS_otu99_otutab_paper1_2021-04-08.tsv",
#                     sep = "\t", stringsAsFactors = F)

# export sample statistics
# create folder for saving
dir.create(file.path("./Output/StatTables"))
exportStats(all.ms, file = paste0("./Output/StatTables/201520162017_CSS_otu99_paper1_transf_stats_",Sys.Date(),".tsv"))
# read with read.csv(file, sep = "\t")

# transform data into long format to merge with fin
# melt asv table into long format
# to do this easily we transpose and melt matrix
css <- melt.data.table(
  setDT(as.data.frame(t(all.mat)), keep.rownames = "Sample"),
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "css.reads"
)

fin<- merge(fin, css, by = c("Sample", "OTU"))

# 8. Calculate relative abundances ---------------------------------------------------------------
# !! not recommended to use !! #
# but calculate for completeness
# calculate relative abundances
fin[, library.size := sum(css.reads), .(Sample)]
fin[, rel.abund := css.reads / library.size]

#sanity check
san <- fin[, .(check = sum(rel.abund)), .(Sample)]
san$check # all 1

# 9. Identify who became absent after quality control ---------------------------------------------------------------
# Classify those as absent after correcting by bin
fin[, n := .N, by = .(ID)] # number of samples by factorial combination
fin[, n.obs := nrow(.SD[css.reads > 0]), by = .(ID)] # how many of those have an actual observation of ASV?

# add z-standardised css.reads to compare rare and abundant things later in regressions
fin[, z.css.reads := (css.reads - mean(css.reads, na.rm = T)) / sd(css.reads, na.rm = T)]

# for all samples, we export the whole file
out <- fin

# calculate mean reads for duplicates
dupl.mean <- fin[, .(reads = mean(reads, na.rm = T),
                  cor.reads = mean(cor.reads, na.rm = T),
                  css.reads = mean(css.reads, na.rm = T), 
                  rel.abund = mean(rel.abund, na.rm = T),
                  z.css.reads = mean(z.css.reads, na.rm = T)), by = .(dr_match_name, dna_type, OTU)]

out <- dupl.mean[fin, c("Sample","Season","year","sample.type.year","catchment.area",
              "distance.from.mouth", "n", "n.obs", "PA", "ID", "library.size") := 
            list(i.Sample,i.Season,i.year,i.sample.type.year,i.catchment.area,
                 i.distance.from.mouth, i.n, i.n.obs, i.PA, i.ID, i.library.size), on = .(dr_match_name, dna_type, OTU)]

# 11. Export ---------------------------------------------------------------
# rearrange data frame before saving
setcolorder(
  out,
  c(
    "ID",
    "Sample",
    "dr_match_name",
    "OTU",
    "year",
    "Season",
    "dna_type",
    "sample.type.year",
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
        paste0("./Objects/201520162017_css_otu99_paper1_",
               Sys.Date(),
               ".rds"))
write.table(out, paste0("./Output/201520162017_css_otu99_paper1_",
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
  paste0("./Output/201520162017_fin_css_otu99_table_paper1_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = T
)

#otu.tab <- read.csv("./Output/201520162017_fin_css_otu99_table_paper1_2021-04-08.csv", sep = ";", dec = ".", stringsAsFactors = F)

# export shrinked meta data
met.df

# only take those samples in OTU table, will remove duplicates and one lost sample
met.df <- met.df[met.df$seq_name %in% levels(factor(row.names(otu.tab))),]
setDT(met.df)
met.df <- met.df[out, c("dr_match_name", "sample.type.year") := 
                 list(i.dr_match_name, i.sample.type.year), on = .(seq_name == Sample)]
setDF(met.df)
# Assign rownames to be Sample ID's
rownames(met.df) <- met.df$seq_name
write.table(
  met.df,
  paste0("./Output/201520162017_meta_otu99_paper1_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = T
)

tax.df <- tax_mat(no.dups)
tax.df <- tax.df[row.names(tax.df) %in% levels(factor(colnames(otu.tab))),]
# orders need to match between tax.tab and otu.tab
tax.df <- tax.df[order(row.names(tax.df)),]

write.table(
  tax.df,
  paste0("./Output/201520162017_tax_otu99_table_paper1_", Sys.Date(), ".csv"),
  sep = ";",
  dec = ".",
  row.names = T
)

#rm(met.df, out, all, all.mat, all.ms, present, san, seqdf, sum.reads, sumdf, tax.df, tax.tab, test,notindna, p, absent,
#   otu.tab, cor.reads, css, dupl.mean, fin)

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