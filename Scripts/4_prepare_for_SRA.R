#-- Script for the publication:
#-- Terrestrial connectivity, upstream aquatic history and seasonality shape bacterial
#-- community assembly within a large boreal aquatic network. The ISME Journal.
#-- Authors: Masumi Stadler & Paul A. del Giorgio
#-- Responsible code author: Masumi Stadler

# This script is the fifth of a series of scripts that were used to analyse the data
# used in the publication.

#---------------------------#
#- Clean meta data for SRA -#
#---------------------------#

# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
  library(tidyverse)
  library(data.table)
  library(readxl)
  
  ### Functions -----------------------------------------------------------------------------
  source("./Functions/custom_fun.R")
  
  fin.df <- select_newest("./Output", "201520162017_fin_css_otu99_phantomcor_paper1_")
  fin.df <- as.matrix(read.csv(
    paste0("./Output/", fin.df),
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ))
  
  splitdf <- readRDS("./Objects/splitdf_new.rds")
  
  # Extract splitdf rows, that are represented in this data set
  sub <- splitdf[splitdf$final_name %in% row.names(fin.df),]
  sub <- sub[sub$final_name != "THW1D", ]
  # Some downriver samples were wrongly assigned to "Downriver"
  sub[grep("SWLR", sub$final_name),]$sample.type.year <- "Soilwater"
  sub[grep("SLR", sub$final_name),]$sample.type.year <- "Soil"
  
  # Save a sub splitdf to upload on Zenodo
  write.table(sub, "./Output/dada2_filt_splitdf.csv", sep = ",", dec = ".", row.names = F)
  
  # Read-in SRA file format
  sra <- read.csv("./Output/Metagenome.environmental.1.0.tsv", sep = "\t", stringsAsFactors = F, skip = 10)
  # must are sample_name, organism, isolation source, collection date, geographic location,
  # latitude, longitude
  # Unique identifiers such as nucleic_acid_type, replicate, sampling_depth(m)
  
  setDT(sub)
  sub[is.na(lat),]
  sub[, lat_lon := paste(paste(lat, "N", sep = " "), paste(abs(long), "W", sep = " "), sep = " ")]
  sub[, collection_date := sampling.date]
  sub[, geo_loc_name := "Canada:Quebec"]
  sub[, ecosystem := sample.type.year]
  sub$ecosystem <- factor(sub$sample.type.year, levels = c("Soil","Sediment",
                                                                        "Soilwater","Hyporheic", 
                                                                        "Wellwater","Stream", "Tributary",
                                                                        "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                        "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                                        "Marine"),
                                    labels = c("Soil","Sediment",
                                               "Soilwater","Hyporheic", 
                                               "Groundwater","Stream", "Tributary",
                                               "Lake", "Pond", "Lake", "Lake",
                                               "River",# "RO3", "RO2", "RO1","Deep",
                                               "Reservoir","Reservoir", "Reservoir","Reservoir",
                                               "River",
                                               "Estuary"))
  sub[, organism := paste(ecosystem, "metagenome", sep = " ")]
  sub[, isolation_source := ecosystem]
  sub[ecosystem == "Stream" | ecosystem == "Tributary" |
        ecosystem == "Lake" | ecosystem == "Pond" | ecosystem == "River" |
        ecosystem == "Reservoir" | ecosystem == "Estuary", sampling_depth_m := 0.3]
  
 # get sampling depths of hypolimnion
  hypo <- read.csv("./MotherData/hypolim_depths.csv"); setDT(hypo)
  sub[hypo, sampling_depth_m := i.sampling_depth, on = c("final_name" = "sample_name")]
  
  # rename RNA to cDNA
  sub[, nucleic_acid_type := dna_type][nucleic_acid_type == "RNA", nucleic_acid_type := "cDNA"]
  
  # keep DNA-RNA match name as sample title
  sub[, sample_title := dr_match_name]
  sub[, sample_name := final_name]
  # correct replicate number of sediment cores
  sub[grep("L330.C1", sample_name), lat_lon := "51.18157 N 63.92587 W"]
  
  # correct replicate number of sediment cores
  sub[grep("L330CM", sample_name), replicate := str_sub(sample_title, start = -1)]
  # duplicate soilwater
  sub[grep("SWLR02", sample_name), replicate := str_sub(sample_title, start = -1)]
  sub[grep("SWLR03", sample_name), replicate := str_sub(sample_title, start = -1)]
  
  # Merge with SRA
  sra <- read.csv("./Output/Metagenome.environmental.1.0.tsv", sep = "\t", skip = 10, header = T)
  colnames(sra)
  # must are:
  # sample_name
  # organism riverine metagenome
  # isolation source = habitat type
  # collection date "YYYY-MM-DD"
  # geographic location
  # latitude and longitude #d.dd N\S d.dd W\E e.g. 38.98 N 77.11 W
  setDT(sra)
  colnames(sra)<-c("sample_name","sample_title","bioproject_accession","organism",
                 "host","isolation_source","collection_date","geo_loc_name",
                 "lat_lon","ref_biomaterial","rel_to_oxygen","samp_collect_device",
                 "samp_mat_process","samp_size","source_material_id","description")
  sra.out <- data.frame(sample_name = sub$sample_name,
             sample_title = sub$sample_title,
             bioproject_accession = "PRJNA693020",
             organism = sub$organism,
             host = NA,
             isolation_source = sub$isolation_source,
             collection_date = sub$collection_date,
             geo_loc_name = sub$geo_loc_name,
             lat_lon = sub$lat_lon,
             ref_biomaterial = NA,
             rel_to_oxygen = NA,
             samp_collect_device = NA,
             samp_mat_process = NA,
             samp_size = NA,
             source_material_id = NA,
             description = NA,
             nucleic_acid_type = sub$nucleic_acid_type,
             replicate = sub$replicate,
             sampling_depth_m = sub$sampling_depth_m)
  
  View(sra.out[sra.out$sample_name %in% non.uni,])
  
  write.table(sra.out, "./Output/sra_lrpaper1_metagenome.tsv",  quote=FALSE, sep='\t',
              row.names = F)
  
  # Meta data ------------------------------------------------------------------------------
  
  meta.ov <- read_excel("./Output/SRA_metadata.xlsx", sheet = 2)
  # 2 rows per sample = forward and reverse
  colnames(meta.ov)
  
  meta <- data.frame(sample_name = sra.out$sample_name,
                     library_ID = sra.out$sample_name,
                     type = sra.out$isolation_source,
                     title = "unique",
                     library_strategy = "AMPLICON",
                     library_source = "METAGENOMIC",
                     library_selection = "PCR",
                     library_layout = "paired",
                     platform = "ILLUMINA",
                     instrument_model = "Illumina MiSeq",
                     design_description = "Conducted by Genome Quebec",
                     filetype = "fastq",
                     filename = "for",
                     filename2 = "rev")
  
  setDT(meta)
  
  meta[, title := paste("16S rRNA (V4) amplicon sequencing of", type, "bacterioplankton community", sep = " ")]
  meta[, filename := paste0(sample_name, "_R1.fastq.gz")]
  meta[, filename2 := paste0(sample_name, "_R2.fastq.gz")]
  
  meta[, type := NULL]
  
  write.table(meta, "./Output/sra_lrpaper1_metadata.tsv",  quote=FALSE, sep='\t',
              row.names = F)
  
  # Update phyloseq object
  pb <- phyloseq(otu_table(fin.df, taxa_are_rows = F),
                 sample_data(met.df),
                 tax_table(tax.tab))
  
  sample.names <- data.frame(stringsAsFactors = F)
  
  for(i in 2015:2017){
    file.names <- dir(paste0("../Raw/",i), pattern = "*fastq.gz")
    t <- strsplit(file.names, "_", fixed = T)
    names <- sapply(t, "[[", 1)
    fr <- sapply(t, '[[', 2)
    #names <- unique(names)
    oneyear <- data.frame(SeqNames = names, ext = fr, Year = i, stringsAsFactors = F)
    sample.names <- rbind(sample.names, oneyear)
  }

# copied files
#copied <- list.files("./Raw/Upload")
#cp <- data.frame(sample = sapply(strsplit(copied, "_"), "[[",1),
#                 ext = sapply(strsplit(copied, "_"), "[[",2))
  

cp <- data.frame(stringsAsFactors = F)
for(i in 2015:2017){
  files <- list.files(path = paste0("./Raw/Upload/",i))
  tmp <- data.frame(sample = sapply(strsplit(files, "_"), "[[", 1),
                    ext = sapply(strsplit(files, "_"), "[[", 2))
  cp <- rbind(cp, tmp)
}

nrow(cp) == nrow(sample.names)

# read in sequence legend
leg <- read.csv("../Meta/seq_sample_legend.csv", sep = ",", dec = ".", stringsAsFactors = F)
nrow(leg) == nrow(cp)/2

sample.names[!sample.names$SeqNames %in% leg$SeqPlateNames,]
leg[!leg$SeqPlateNames %in% sample.names$SeqNames,]
# add D for DNA in 2017 in leg
# add R for RNA in 2017 in leg
# exchange underscore to dot for 2017 in leg
# exchange hyphen to dot for 2016 in sample.names
leg$RawSeqNames <- leg$SeqPlateNames
leg$SeqPlateNames[leg$Year == 2017 & leg$DnaType == "DNA"] <- paste0(leg$SeqPlateNames[leg$Year == 2017 & leg$DnaType == "DNA"], "D")
leg$SeqPlateNames[leg$Year == 2017 & leg$DnaType == "cDNA"] <- paste0(leg$SeqPlateNames[leg$Year == 2017 & leg$DnaType == "cDNA"], "R")
leg$SeqPlateNames[leg$Year == 2017] <- gsub("[_]", ".", leg$SeqPlateNames[leg$Year == 2017])
sample.names$AsDir <- sample.names$SeqNames
sample.names$SeqNames[sample.names$Year == 2016] <- gsub("[-]", ".", sample.names$SeqNames[sample.names$Year == 2016])

# do all sequence names match?
sample.names[!sample.names$SeqNames %in% leg$SeqPlateNames,]
leg[!leg$SeqPlateNames %in% sample.names$SeqNames,]

sample.names$fullname <- paste(sample.names$SeqNames, sample.names$ext, sep = "_")
#for merging
sample.names$fullseq <- paste(sample.names$AsDir, sample.names$ext, sep = "_")
cp$fullseq <- paste(cp$sample, cp$ext, sep = "_")

# no duplicated samples names
sample.names[duplicated(sample.names$fullname),]


t <- left_join(cp, sample.names %>% select(SeqNames, fullseq), by = c("fullseq"))

# Read in meta data
meta <-
  read.csv(
    "./MotherData/master_bacmeta.csv",
    sep = ";",
    dec = ".",
    stringsAsFactors = F
  )
nrow(meta); length(unique(t$SeqNames))
from<-which(meta$DadaNames == "S33R")+1
cor <- rbind(meta[1:which(meta$DadaNames == "S33R"),],
      meta[which(meta$DadaNames == "S33R"),] %>% mutate(MatchNames = "S33D",
                                                        GQNames = "S33D",
                                                        DadaNames = "S33D",
                                                        DnaType = "DNA"),
      meta[from:nrow(meta),])


tt <- data.frame(meta = cor$DadaNames[order(cor$DadaNames)], plate = unique(t$SeqNames))
tt[tt$meta != tt$plate,]

all <- left_join(t %>% select(SeqNames) %>% distinct(), cor, by = c("SeqNames" = "DadaNames"))

# duplicates?
(dups <- all[duplicated(all$MatchNames),]$MatchNames)

setDT(all)
all[MatchNames %in% dups,]$MatchNames
overwrite <- c("PR29.1D","PR29.2D","PR30.1D","PR30.2D","PR31.1D","PR31.2D",
               "RO106.1D","RO106.2D","RO230.1D","RO230.2D","TR49.1D","TR49.2D")
# make new column
all[, NewNames := MatchNames]
all[NewNames %in% dups, NewNames := overwrite]
# any still duplicated?
all[duplicated(all$NewNames),]$NewNames

new.names <- rep(all$NewNames, each = 2)
last.ch <- data.frame(new.names, t$ext, copied)

overwrite <- paste(last.ch$new.names, last.ch$t.ext, sep = "_")

filepath <- "./Raw/Upload/"

paste0(filepath, copied[100])
paste0(filepath, overwrite[100])

for(f in 1:length(copied)){
  file.rename(paste0(filepath, copied[f]), paste0(filepath, overwrite[f]))  
}

last.ch <- last.ch %>%
  mutate(new.names = paste(new.names, t.ext, sep = "_"),
         old.names = copied,
         t.ext = NULL,
         copied = NULL)
write.table(last.ch, "./Raw/old-new.filesnames.csv", sep = ";", row.names = F)

# number of samples?
files <- list.files("./Raw/Upload"); length(files)
# Move files that are not part of this paper
# mv *DSeq*.fastq.gz ../Deep
# mv *BIOD*.fastq.gz ../Bioassay

sra <- read.csv("./Raw/Metagenome.environmental.1.0.tsv", sep = "\t", skip = 10, header = T)
colnames(sra)
# must are:
# sample_name
# organism riverine metagenome
# isolation source = habitat type
# collection date "YYYY-MM-DD"
# geographic location
# latitude and longitude #d.dd N\S d.dd W\E e.g. 38.98 N 77.11 W

#copied <- list.files("./Raw/Upload", pattern = "fastq.gz")

View(all[is.na(all$lat),])

sra.out <- all %>% filter(SeqDepth == "Shallow") %>%
    select(sample_name = NewNames, old_sample = MatchNames,
           isolation_source = sample.type.year,
                 collection_date = sampling.date,
                 lat, long, sample_description = DnaType)
  
  library(data.table)
  setDT(sra.out)
  # add lat long
  sra.out[!is.na(lat), latitude_and_longitude := paste(lat, "N", abs(long), "W")]
  # clean isolation source
  sra.out[isolation_source == "Deep", isolation_source := "Reservoir_Hypolimnion"]
  sra.out[isolation_source == "Downriver" | isolation_source == "Upriver", isolation_source := "River"]
  sra.out[isolation_source == "IslandLake" | isolation_source == "HeadwaterLakes", isolation_source := "Lake"]
  sra.out[isolation_source == "PRLake", isolation_source := "Pond"]
  sra.out[isolation_source == "Marine", isolation_source := "Estuary"]
  sra.out[isolation_source == "RO1" |
            isolation_source == "RO2" |
            isolation_source == "RO3", isolation_source := "Reservoir"]
  
  sra.out <- sra.out[isolation_source != "Bioassay",]
  
  sra.out[, organism := paste(isolation_source, "metagenome", sep = " ")]
  sra.out[isolation_source == "Reservoir_Hypolimnion", organism := "Reservoir metagenome"]
  sra.out[, geographic_location := "Canada:Quebec"]
  
 
#write.table(sra.out, "./Raw/SRA.out.table.csv", sep = ",", row.names = F)

hand <- read.csv("./Raw/SRA.out.table_handcor.csv", sep = ",", stringsAsFactors = F)
hand <- hand %>% distinct ()
setDT(hand)

# remove THW1D
sra.out <- sra.out[sample_name != "THW1D",]

# rm *THW1D*.fastq.gz

# which don't have latitudes or dates?
no.date <- sra.out[is.na(collection_date),]$sample_name; length(no.date)
no.lat <- sra.out[is.na(latitude_and_longitude),]$sample_name; length(no.lat)
nrow(hand[sample_name %in% no.date,])
head(hand[sample_name %in% no.lat,]); nrow(hand[sample_name %in% no.lat,])

# just take the date and lat long from the manually corrected ones
sra.out <- sra.out[hand[sample_name %in% no.date,], collection_date := i.collection_date,
                on = c(old_sample = "sample_name")]
sra.out <- sra.out[hand[sample_name %in% no.lat,], latitude_and_longitude := i.latitude_and_longitude,
                   on = c(old_sample = "sample_name")]
sra.out[sample_name %in% no.lat, ]


#order by file appearance
files <- unique(sapply(strsplit(list.files("./Raw/Upload"), "_"), "[[",1))
out <- sra.out[match(files, sample_name),]
# any duplicated?
out[duplicated(out$sample_name),]

out <- out %>% select(sample_name, organism, isolation_source, collection_date, 
                              geographic_location, latitude_and_longitude, nucleci_acid_type = sample_description)
# add attributes
setDT(out)
out[, replicate := 1]
out[grep(".2D", out$sample_name, fixed = T), replicate := 2]

# Depth
depths <- read.csv("hypolim_depths.csv", sep = ",", stringsAsFactors = F)
out[grep("90m", out$sample_name), isolation_source := "Reservoir_Hypolimnion"]
out[grep("60m", out$sample_name), isolation_source := "Reservoir_Hypolimnion"]
out[sample_name == "RO211160mD", sample_name := "RO211190mD"]
out[, Sampling_Depth_m := 0.5]
out[depths, Sampling_Depth_m := i.sampling_depth, on = .(sample_name)]


#write.table(out, "./Raw/SRA.out.table_final.csv", sep = ",", row.names = F)
out<- read.csv("./Raw/SRA.out.table_final.csv", sep = ",", stringsAsFactors = F)
# Meta data
library(readxl)
sra.meta <- read_xlsx("./Raw/SRA_metadata.xlsx", sheet = 2)

# must are:
# "sample_name"        "library_ID"         "title"             
# "library_strategy"   "library_source"     "library_selection" 
# "library_layout"     "platform"           "instrument_model"  
# "design_description" "filetype"           "filename"

fwd <- list.files("./Raw/Upload/", pattern = "*R1.fastq.gz")
rvs <- list.files("./Raw/Upload/", pattern = "*R2.fastq.gz")

View(data.frame(c(fwd,1), rvs))

sra.meta <- data.frame(sample_name = out$sample_name,
           library_ID = out$sample_name,
           title = "MiSeq sequencing of 16S rRNA (V4 region) of bacterial community",
           library_strategy = "AMPLICON",
           library_source = "METAGENOMIC",
           library_selection = "PCR",
           library_layout = "PAIRED-END",
           platform = "ILLUMINA",
           instrument_model = "Illumina MiSeq",
           design_description = "Sequencing conducted by Genome Quebec",
           filetype = "fastq",
           filename = fwd,
           filename2 = rvs)

sra.meta[, sample_name := out$sample_name]
