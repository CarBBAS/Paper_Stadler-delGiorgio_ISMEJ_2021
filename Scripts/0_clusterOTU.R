
#----------#
# PACKAGES #
#----------#s
library(DECIPHER)
library(Biostrings)
library(tidyverse)
library(plyr)
library(doMC)

#-----------#
# FUNCTIONS #
#-----------#
source("./Functions/custom_fun.R")
#-----------------#
# PARALLEL SET-UP #
#-----------------#
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

#-----------#
# Read data #
#-----------#
# Read in dada2 output
seqtab <- readRDS("./MotherData/nochim_seqtab_2018.rds") # ASV table with raw sequences
tax <-
  readRDS("./MotherData/taxtab_idtaxa_silva_v138_2018.rds") # assigned taxonomy

# extract taxonomy
tax <- as.data.frame(tax, stringsAsFactors = F)
# split into taxnomy group, which ASVs belong to the same taxonomical classification
align.group <- unite(tax, col = "align.group", na.rm = T)
align.group$Seq <- rownames(align.group)

# extract sequences into list bins, where bins are the deepest taxonomical classification
align.ls <- dlply(align.group, .(align.group), function(x){
  vec <- x$Seq
  return(vec)
}, .parallel = T)

# do not include un-classified ASVs, they will be removed downstream anyway
align.ls[[1]] <- NULL

# set to number of cpus/processors to use for the clustering
nproc <- (detectCores() / 2) - 2 

#align.ls <- align.ls[2:length(align.ls)]

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
      cutoff = 0.01, # use `cutoff = 0.01` for a 99% OTU ; 0.003 for 99.7%
      processors = nproc)
    
    cl.out[[i]] <- clusters
    
  } else {
    cl.out[[i]] <- align.ls[[i]]
    names(cl.out[i]) <- names(align.ls[i])
  }
}

saveRDS(cl.out, "./Objects/OTU_99_clusters.rds")
cl.out <- readRDS("./Objects/OTU_clusters.rds")

for(i in 1:length(cl.out)){
  if(class(cl.out[[i]]) == 'data.frame'){
    cl.out[[i]]$sequences <- align.ls[[i]]
  } else{
    cl.out[[i]] <- data.frame(cluster = NA, sequences = align.ls[[i]])
  }
}

fin.ls <- list()
for(i in 1:length(cl.out)){
  sub.seq <- as.data.frame(t(seqtab[,colnames(seqtab) %in% align.ls[[i]]]))
  setDT(sub.seq, keep.rownames = "sequences", key = "sequences")
  
  merg.seq <- sub.seq[setDT(cl.out[[i]], key = "sequences")]
  merg.seq <- merg.seq[, lapply(.SD[,-1], sum, na.rm = T), by = .(cluster)]
  
  merg.seq$tax <- names(align.ls)[i]
  fin.ls[[i]] <- merg.seq
}

fin.df <- bind_rows(fin.ls) #create OTU numbers

fin.df$OTU <-
  paste("OTU", seq(length = nrow(fin.df)), sep = "_")
setDF(fin.df)
row.names(fin.df) <- fin.df$OTU

tax <- data.frame(OTU = fin.df$OTU, all = fin.df$tax)
tax <- tax %>% separate(col = "all", into = c("domain","phylum","class","order","family","genus","species"),
                sep = "[_]")
rownames(tax) <- tax$OTU ; tax$OTU <- NULL
tax <- as.matrix(tax)

fin.df$cluster <- NULL; fin.df$OTU <- NULL; fin.df$tax <- NULL
seqtab <- as.matrix(t(fin.df))



saveRDS(tax, "./Objects/otu_taxtab_99.rds")
saveRDS(seqtab, "./Objects/otu_seqtab_99.rds")
