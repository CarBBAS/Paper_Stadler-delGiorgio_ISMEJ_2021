library(tidyverse)
library(data.table)

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
copied <- list.files("./Raw/Upload")
cp <- data.frame(sample = sapply(strsplit(copied, "_"), "[[",1),
                 ext = sapply(strsplit(copied, "_"), "[[",2))

# read in sequence legend
leg <- read.csv("../Meta/seq_sample_legend.csv", sep = ",", dec = ".", stringsAsFactors = F)

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
out[grep("90m", out$sample_name), isolation_source := "Reservoir_Hypolimnion"]
out[grep("60m", out$sample_name), isolation_source := "Reservoir_Hypolimnion"]
out[sample_name == "RO211160mD", sample_name := "RO211190mD"]
out[, Depth_m := 0.5]


write.table(out, "./Raw/SRA.out.table_final.csv", sep = ",", row.names = F)

out <- read.csv("./Raw/SRA.out.table_final.csv", sep = ",")

hand <- left_join(hand, sra.out %>% select(sample_name, sample_description))
View(t[t$sample_description == "DNA",])
  