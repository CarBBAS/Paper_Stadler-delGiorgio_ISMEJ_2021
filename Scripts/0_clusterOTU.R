#-- Script for the publication:
#-- Title
#-- Responsible author: Masumi Stadler

# This script is the first of a series of scripts that were used to analyse the data
# used in the publication.

###--------------------------###
#-   OTU clustering of ASVs  - #
###--------------------------###

# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("DECIPHER", "Biostrings", # OTU clustering
              "tidyverse", "plyr", "data.table", # wrangling
              "doMC") # parallel computing

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

# 2. Read and prepare data ---------------------------------------------------------------

# Read in dada2 output
seqtab <- readRDS("./MotherData/nochim_seqtab_2018.rds") # ASV table with raw sequences
tax <-
  readRDS("./MotherData/taxtab_idtaxa_silva_v138_2018.rds") # assigned taxonomy

# extract taxonomy
tax <- as.data.frame(tax, stringsAsFactors = F)
# split into taxnomy group, which ASVs belong to the same taxonomical classification
align.group <- unite(tax, col = "align.group", sep = "$", na.rm = T)
align.group$Seq <- rownames(align.group)

# extract sequences into list bins, where bins are the deepest taxonomical classification
# bin names are taxonomical classifications
align.ls <- dlply(align.group, .(align.group), function(x){
  vec <- x$Seq
  return(vec)
}, .parallel = T)

# do not include un-classified ASVs, they will be removed downstream anyway
align.ls[[1]] <- NULL


# 3. OTU clustering ---------------------------------------------------------------
# set to number of cpus/processors to use for the clustering
nproc <- (detectCores() / 2) - 2

#align.ls <- align.ls[2:length(align.ls)]

# create empty list to fill
cl.out <- list()

# run loop to do sequence alignment and similarity collapse at 99% within each taxonomical classification
for(i in 1:length(align.ls)){
  if(length(align.ls[[i]]) >= 2){
    sub.seq <- seqtab[,colnames(seqtab) %in% align.ls[[i]]]
    asv_sequences <- colnames(sub.seq)
    sample_names <- rownames(sub.seq)
    dna <- Biostrings::DNAStringSet(asv_sequences)
    
    ## Find clusters of ASVs to form the new OTUs
    aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
    d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
    clusters <- DECIPHER::IdClusters(
      d, 
      method = "complete",
      cutoff = 0.01,
      # use `cutoff = 0.01` for a 99% OTU ; 0.003 for 99.7%
      # cutoff 0.02 for 98%
      # cutoff 0.03 for 97%
      processors = nproc)
    
    cl.out[[i]] <- clusters
    
  } else {
    cl.out[[i]] <- align.ls[[i]]
    names(cl.out[i]) <- names(align.ls[i])
  }
}

# at this point we have a list with
# each bin representing a taxonomic group (e.g. 'Bacteria') with individual ASVs
# the bin contains the OTU number each ASV was assigned to
# OTU numbers go from 1 to maximum OTU number within each taxonomic group
# Meaning that there are many OTU1, OTU2, OTU3, important is that they stay within the bin
# Actual OTU assignment will happen later
#saveRDS(cl.out, "./Objects/OTU_clusters_95.rds") # save intermediate object
saveRDS(cl.out, "./Objects/OTU_clusters.rds") # save intermediate object

for(j in 98:95){
  print(paste0("Working on...",j,"% similarity"))
  # 2. Assign OTU numbers and keep ASV sequences  ---------------------------------------------------------------
  cl.out <- readRDS(paste0("./Objects/OTU_clusters_",j,".rds"))
  
  # This loop re-assigns the exact sequences to each "OTU" 
  for(n in 1:length(cl.out)){
    if(class(cl.out[[n]]) == 'data.frame'){
      # if there are many "ASVs" in a taxonomic classification, output is a data frame
      cl.out[[n]]$sequences <- align.ls[[n]] # fill with sequence
    } else{
      # if there is only one ASV in taxonomic classification, make it a data frame for easy merging
      cl.out[[n]] <- data.frame(cluster = NA, sequences = align.ls[[n]])
    }
  }
  
  # create new empty list
  # loops runs bin by bin, summing all reads that are of the same cluster (= OTU)
  fin.ls <- list()
  for(i in 1:length(cl.out)){
    #transpose back into OTU table and keep sequences as rownames
    sub.seq <- as.data.frame(t(seqtab[,colnames(seqtab) %in% align.ls[[i]]]))
    setDT(sub.seq, keep.rownames = "sequences", key = "sequences")
    
    # merge abundance table with clustering output
    # combines both tables by exact sequence matches, and adds the cluster ID
    merg.seq <- sub.seq[setDT(cl.out[[i]], key = "sequences")]
    # extract sequence of most abundant ASV in OTU cluster
    max.abun.seq <- merg.seq[, `:=`(sum = rowSums(.SD, na.rm = T)), .SDcols = -c("sequences","cluster")]
    max.seq <- max.abun.seq[, .SD[which.max(sum)], by = .(cluster)]$sequences
    
    # sum abundances for all ASVs within one OTU
    merg.seq <- merg.seq[, lapply(.SD[,-1], sum, na.rm = T), by = .(cluster)] 
    
    # Add taxonomic bin info
    merg.seq$tax <- names(align.ls)[i]
    merg.seq$sequences <- max.seq
    fin.ls[[i]] <- merg.seq
  }
  
  # bind lists into one big data frame
  fin.df <- bind_rows(fin.ls)
  
  # order rows by max sum abundance
  fin.df <- fin.df[order(fin.df$sum, decreasing = T),]
  
  #create OTU numbers = ascending numbers (smallest number = OTU with highest abundance)
  fin.df$OTU <-
    paste("OTU", seq(length = nrow(fin.df)), sep = "_")
  setDF(fin.df)
  row.names(fin.df) <- fin.df$OTU
  
  # save OTU table with sequences separately (as reference)
  ref.tab <- fin.df[,c("OTU","sequences")]
  
  # save taxonomic annotation for each OTU
  tax <- data.frame(OTU = fin.df$OTU, all = fin.df$tax)
  tax <- tax %>% separate(col = "all", into = c("domain","phylum","class","order","family","genus","species"),
                          sep = "[$]")
  rownames(tax) <- tax$OTU ; tax$OTU <- NULL
  tax <- as.matrix(tax)
  
  # remove unneccessary columns, convert back to OTU table
  fin.df$cluster <- NULL; fin.df$sum <- NULL; fin.df$OTU <- NULL; fin.df$tax <- NULL; fin.df$sequences <- NULL
  seqtab <- as.matrix(t(fin.df))
  
  # save as R object
  saveRDS(tax, paste0("./Objects/OTU_",j,"_taxonomy.rds"))
  saveRDS(seqtab, paste0("./Objects/OTU_",j,"_table.rds"))
  
  # save as csv file for repository archiving
  write.table(seqtab, paste0("./Output/OTU_",j,"_table.csv"), sep = ",", row.names = T)
  write.table(tax, paste0("./Output/OTU_",j,"_taxonomy.csv"), sep = ",", row.names = T)
  write.table(ref.tab, paste0("./Output/OTU_",j,"_sequences.csv"), sep = ",", row.names = F)
  
}

# 2. Assign OTU numbers and keep ASV sequences  ---------------------------------------------------------------
cl.out <- readRDS("./Objects/OTU_clusters.rds")

# This loop re-assigns the exact sequences to each "OTU" 
for(i in 1:length(cl.out)){
  if(class(cl.out[[i]]) == 'data.frame'){
    # if there are many "ASVs" in a taxonomic classification, output is a data frame
    cl.out[[i]]$sequences <- align.ls[[i]] # fill with sequence
  } else{
    # if there is only one ASV in taxonomic classification, make it a data frame for easy merging
    cl.out[[i]] <- data.frame(cluster = NA, sequences = align.ls[[i]])
  }
}

# create new empty list
# loops runs bin by bin, summing all reads that are of the same cluster (= OTU)
fin.ls <- list()
for(i in 1:length(cl.out)){
  #transpose back into OTU table and keep sequences as rownames
  sub.seq <- as.data.frame(t(seqtab[,colnames(seqtab) %in% align.ls[[i]]]))
  setDT(sub.seq, keep.rownames = "sequences", key = "sequences")
  
  # merge abundance table with clustering output
  # combines both tables by exact sequence matches, and adds the cluster ID
  merg.seq <- sub.seq[setDT(cl.out[[i]], key = "sequences")]
  # extract sequence of most abundant ASV in OTU cluster
  max.abun.seq <- merg.seq[, `:=`(sum = rowSums(.SD, na.rm = T)), .SDcols = -c("sequences","cluster")]
  max.seq <- max.abun.seq[, .SD[which.max(sum)], by = .(cluster)]$sequences
  
  # sum abundances for all ASVs within one OTU
  merg.seq <- merg.seq[, lapply(.SD[,-1], sum, na.rm = T), by = .(cluster)] 
  
  # Add taxonomic bin info
  merg.seq$tax <- names(align.ls)[i]
  merg.seq$sequences <- max.seq
  fin.ls[[i]] <- merg.seq
}

# bind lists into one big data frame
fin.df <- bind_rows(fin.ls)

# order rows by max sum abundance
fin.df <- fin.df[order(fin.df$sum, decreasing = T),]

#create OTU numbers = ascending numbers (smallest number = OTU with highest abundance)
fin.df$OTU <-
  paste("OTU", seq(length = nrow(fin.df)), sep = "_")
setDF(fin.df)
row.names(fin.df) <- fin.df$OTU

# save OTU table with sequences separately (as reference)
ref.tab <- fin.df[,c("OTU","sequences")]

# save taxonomic annotation for each OTU
tax <- data.frame(OTU = fin.df$OTU, all = fin.df$tax)
tax <- tax %>% separate(col = "all", into = c("domain","phylum","class","order","family","genus","species"),
                sep = "[$]")
rownames(tax) <- tax$OTU ; tax$OTU <- NULL
tax <- as.matrix(tax)

# remove unneccessary columns, convert back to OTU table
fin.df$cluster <- NULL; fin.df$sum <- NULL; fin.df$OTU <- NULL; fin.df$tax <- NULL; fin.df$sequences <- NULL
seqtab <- as.matrix(t(fin.df))

# save as R object
saveRDS(tax, "./Objects/OTU_99_taxonomy.rds")
saveRDS(seqtab, "./Objects/OTU_99_table.rds")

# save as csv file for repository archiving
write.table(seqtab, "./Output/OTU_99_table.csv", sep = ",", row.names = T)
write.table(tax, "./Output/OTU_99_taxonomy.csv", sep = ",", row.names = T)
write.table(ref.tab, "./Output/OTU_99_sequences.csv", sep = ",", row.names = F)

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