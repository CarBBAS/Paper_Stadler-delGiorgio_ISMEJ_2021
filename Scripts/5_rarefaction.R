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
if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
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

# # Raw pre-CSS
# pb <- readRDS("./Objects/physeq_paper1_pre-css.rds")
# #pb <- prune_taxa(!taxa_sums(pb) == 0, pb)
# pb <- prune_samples(!sample_sums(pb) == 0, pb)
# 
# # remove THW1D - no meta data
# pb <- subset_samples(pb, seq_name != "THW1D")
# met.df <- sample_df(pb)
# 
# # Phantom taxa ---------------------------------------------------------------------------------------
# otu.tab <- setDT(as.data.frame(otu_table(pb), strings.As.Factors = F), keep.rownames = "Sample")
# # melt OTU table into long format
# otu.tab <- melt.data.table(
#   otu.tab,
#   id.vars = "Sample",
#   measure.vars = patterns("^OTU_"),
#   variable.name = "OTU",
#   value.name = "reads"
# )
# 
# # Join OTU and meta table
# sumdf <- left_join(
#   otu.tab,
#   met.df[met.df$seq_name %in% otu.tab$Sample,] %>%
#     dplyr::select(
#       seq_name,
#       dr_match_name,
#       year,
#       Season,
#       dna_type,
#       replicate,
#       sample.type,
#       sample.type.year,
#     ),
#   by = c("Sample" = "seq_name")
# )
# 
# # set back to data.table, order data.table by catchment.area
# setDT(sumdf)
# 
# # check for duplicates
# any(duplicated(sumdf[dna_type == "DNA" & OTU == "OTU_10",]$dr_match_name) == T) # yes
# any(duplicated(sumdf[dna_type == "RNA" & OTU == "OTU_10",]$dr_match_name) == T) # no
# 
# # calculate mean for duplicates
# # calculate mean reads for duplicates
# dupl.mean <- sumdf[, .(mean.reads = mean(reads, na.rm = T)), by = .(dr_match_name, dna_type, OTU)]
# 
# # check for duplicates
# any(duplicated(dupl.mean[dna_type == "DNA" & OTU == "OTU_10",]$dr_match_name) == T) # no
# 
# # Rename dna type factor
# dupl.mean[, dna_type := factor(dna_type, levels = c("DNA","cDNA"), labels = c("DNA","RNA"))]
# 
# # cast into wide format
# castdf <- dcast(dupl.mean, dr_match_name + OTU ~ dna_type, value.var = "mean.reads")
# 
# # fill NAs with 0, those are taxa not found in RNA or DNA
# castdf[is.na(DNA), DNA := 0]
# castdf[is.na(RNA), RNA := 0]
# 
# # how many samples have an OTU with RNA > 0 and DNA == 0? (phantom taxa)
# castdf[DNA > 0 | RNA > 0, n.all := .N, .(dr_match_name)] # number of any observations above 0 by sample
# temp <- castdf[RNA > 0 & DNA == 0,]
# temp <- temp[, .(n = .N, n.all = unique(n.all)), .(dr_match_name)][, prop := n * 100 / n.all]
# 
# setDT(met.df)
# # which samples have more than 50% phantom taxa
# View(met.df[dr_match_name %in% temp[prop != 100 & prop > 50,]$dr_match_name,] %>% select(dr_match_name,dna_type,LibrarySize))
# 
# # overwrite all DNA observations == 0, where RNA > 0 with 1
# castdf[RNA > 0 & DNA == 0, DNA := 1]
# castdf[RNA > 0 & DNA == 0,] # check if overwrite was successful
# 
# # format back to long format
# temp <- melt.data.table(castdf,
#                 id.vars = c("dr_match_name","OTU"),
#                 measure.vars = c("DNA","RNA"),
#                 variable.name = "dna_type",
#                 value.name = "reads")
# 
# # original and casted data frame have not the same length
# # possible we added rows for samples that didn't have any RNA, but with RNA = 0
# nrow(temp) == nrow(dupl.mean)
# 
# # do row merge/overwrite
# # we have dealt with the replicates already in the earlier script
# mer <- dupl.mean[temp, cor.reads := i.reads, on = .(dr_match_name, OTU, dna_type)]
# # Rename dna type factor
# met.df[, dna_type := factor(dna_type, levels = c("DNA","cDNA"), labels = c("DNA","RNA"))]
# mer <- mer[met.df, Sample := i.seq_name, on = .(dr_match_name, dna_type)]
# 
# # cast into wide format
# fin.df <- dcast(mer, Sample ~ OTU, value.var = "cor.reads")
# # make matrix and give row.names
# fin.df <- as.matrix(setDF(fin.df, rownames = fin.df$Sample)[,-1]) # remove first row with "Sample" move into rownames
# 
# # same dimensions as original df?
# dim(fin.df)
# dim(otu_table(pb)) # yes
# 
# # Read in the rest and update physeq object
# # Taxonomy table
# tax.tab <- select_newest("./Output", "201520162017_tax_otu99_table_paper1_")
# tax.tab <-
#   as.matrix(read.csv(
#     paste0("./Output/", tax.tab),
#     sep = ";",
#     dec = ".",
#     row.names = 1,
#     stringsAsFactors = F
#   ))
# tax.tab[1:4,]
# # orders need to match between tax.tab and otu.tab
# tax.tab <- tax.tab[order(row.names(tax.tab)),]
# #tax.tab <- tax.tab[,1:1000]
# 
# # Meta data
# met.df <-
#   select_newest(path = "./Output", file.pattern = "201520162017_meta_otu99_paper1_")
# met.df <-
#   read.csv(
#     paste0("./Output/", met.df),
#     sep = ";",
#     dec = ".",
#     row.names = 1,
#     stringsAsFactors = F
#   )
# 
# met.df[1:4,]
# # Checked by PCoA THW1D is a soilwater sample (actually Hyporheic)
# met.df[met.df$seq_name == "THW1D", ]$sample.type.year <- "Hyporheic"
# # Some downriver samples were wrongly assigned to "Downriver"
# met.df[grep("SWLR", met.df$seq_name),]$sample.type.year <- "Soilwater"
# met.df[grep("SLR", met.df$seq_name),]$sample.type.year <- "Soil"
# 
# 
# # merge some sample types
# met.df$sample.type.year <- factor(met.df$sample.type.year, levels = c("Soil","Sediment",
#                                                                       "Soilwater","Hyporheic", 
#                                                                       "Wellwater","Stream", "Tributary",
#                                                                       "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
#                                                                       "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
#                                                                       "Marine"),
#                                   labels = c("Soil","Sediment",
#                                              "Soilwater","Soilwater", 
#                                              "Groundwater","Stream", "Tributary",
#                                              "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
#                                              "Upriver",# "RO3", "RO2", "RO1","Deep",
#                                              "Reservoirs","Reservoirs", "Reservoirs","Reservoirs",
#                                              "Downriver",
#                                              "Estuary"))
# 
# met.df$Season <- factor(met.df$Season, levels = c("Spring","Summer","Autumn"),
#                         labels = c("Spring","Summer","Autumn"))
# 
# met.df$dna_type <- factor(met.df$dna_type, levels = c("DNA","cDNA"),
#                           labels = c("DNA","RNA"))
# 
# 
# # Update phyloseq object
# pb <- phyloseq(otu_table(fin.df, taxa_are_rows = F),
#                sample_data(met.df),
#                tax_table(tax.tab))

# save physeq object
#saveRDS(pb, "./Objects/physeq_paper1_pre-css_phantomcor.rds")

pb <- readRDS("./Objects/physeq_paper1_pre-css_phantomcor.rds")

#' # do we have several files per object? -> take newest version
#' # read in
#' fin.df <- select_newest("./Output", "201520162017_fin_css_otu99_phantomcor_paper1_")
#' fin.df <- as.matrix(read.csv(
#'   paste0("./Output/", fin.df),
#'   sep = ";",
#'   dec = ".",
#'   row.names = 1,
#'   stringsAsFactors = F
#' ))
#' 
#' # Taxonomy table
#' tax.tab <- select_newest("./Output", "201520162017_tax_otu99_table_paper1_")
#' tax.tab <-
#'   as.matrix(read.csv(
#'     paste0("./Output/", tax.tab),
#'     sep = ";",
#'     dec = ".",
#'     row.names = 1,
#'     stringsAsFactors = F
#'   ))
#' tax.tab[1:4,]
#' # orders need to match between tax.tab and otu.tab
#' tax.tab <- tax.tab[order(row.names(tax.tab)),]
#' #tax.tab <- tax.tab[,1:1000]
#' 
#' # Meta data
#' met.df <-
#'   select_newest(path = "./Output", file.pattern = "201520162017_meta_otu99_paper1_")
#' met.df <-
#'   read.csv(
#'     paste0("./Output/", met.df),
#'     sep = ";",
#'     dec = ".",
#'     row.names = 1,
#'     stringsAsFactors = F
#'   )
#' 
#' met.df[1:4,]
#' # Checked by PCoA THW1D is a soilwater sample (actually Hyporheic)
#' met.df[met.df$seq_name == "THW1D", ]$sample.type.year <- "Hyporheic"
#' # Some downriver samples were wrongly assigned to "Downriver"
#' met.df[grep("SWLR", met.df$seq_name),]$sample.type.year <- "Soilwater"
#' met.df[grep("SLR", met.df$seq_name),]$sample.type.year <- "Soil"
#' 
#' 
#' # merge some sample types
#' met.df$sample.type.year <- factor(met.df$sample.type.year, levels = c("Soil","Sediment",
#'                                                                       "Soilwater","Hyporheic", 
#'                                                                       "Wellwater","Stream", "Tributary",
#'                                                                       "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
#'                                                                       "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
#'                                                                       "Marine"),
#'                                   labels = c("Soil","Sediment",
#'                                              "Soilwater","Soilwater", 
#'                                              "Groundwater","Stream", "Tributary",
#'                                              "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
#'                                              "Upriver",# "RO3", "RO2", "RO1","Deep",
#'                                              "Reservoirs","Reservoirs", "Reservoirs","Reservoirs",
#'                                              "Downriver",
#'                                              "Estuary"))
#' 
#' met.df$Season <- factor(met.df$Season, levels = c("Spring","Summer","Autumn"),
#'                         labels = c("Spring","Summer","Autumn"))
#' 
#' met.df$dna_type <- factor(met.df$dna_type, levels = c("DNA","cDNA"),
#'                           labels = c("DNA","RNA"))
#' 
#' # Update phyloseq object
#' pb <- phyloseq(otu_table(fin.df, taxa_are_rows = F),
#'                sample_data(met.df),
#'                tax_table(tax.tab))
#' 
#' # remove THW1D - no meta data
#' pb <- subset_samples(pb, seq_name != "THW1D")
#' pb <- prune_taxa(taxa_sums(pb) != 0, pb)
#' pb <- prune_samples(sample_sums(pb) != 0, pb) # remove samples with no reads
#' 
#' # export factors for colouring
#' sample.factors <- levels(met.df$sample.type.year)
#' 
#' # more colour blind friendly
#' colvec <- c("#FCFDBFFF", #"#FEC589FF", #Soil
#'             "#FEC589FF", #"#FDEBACFF", #Sediment
#'             "#F9795DFF", #"#F9795DFF", #Soilwater
#'             "#DE4968FF", #"#DE4968FF", #Groundwater,
#'             "skyblue", #Stream
#'             "#AD347CFF",# Tributary, 
#'             "palegreen", #Riverine Lakes, 
#'             "#7AD151FF", #Headwater Ponds,
#'             "#FDE725FF",# Lake, 
#'             "#1F9F88FF", # Upriver, 
#'             "orchid", #"#471063FF", #Reservoir,
#'             #'salmon',
#'             #"pink",
#'             #'navy',
#'             "#375A8CFF", #Downriver,
#'             "gray40") #Estuary)
#' 
#' names(colvec) <- as.character(sample.factors) # assign names for later easy selection

# set theme for plotting
theme_set(theme_bw())

# export factors for colouring
sample.factors <- levels(sample_df(pb)$sample.type.year)

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

# Re-do all analyses with rarefaction tresholds --------------------------------------------------------------
# Do not run this code if you have a slow computer, or if you want to go quickly through the script
# Read in intermediate data frame in line 214 onwards
# df <- sample_df(pb)
# hist(df$LibrarySize, breaks = 20)
# 
# #choose rarefaction thresholds
# min_lib <- c(min(df$LibrarySize), 5000, 10000, 20000, 30000)
# set.seed(3)
# rar.ls <- list()
# for(i in 1:length(min_lib)){
#   out <- rarefy_even_depth(pb, sample.size = min_lib[i], rngseed = 3)
#   rar.ls[[i]] <- out
#   names(rar.ls)[i] <- paste("rar",min_lib[i], sep = "_")
# }
# 
# # amend original matrix
# rar.ls[[length(rar.ls)+1]] <- pb
# names(rar.ls)[length(rar.ls)] <- "orig"
# 
# # save rarefaction list
# saveRDS(rar.ls, "./Objects/rarefaction_list_phyloseq_precss.rds")

rar.ls <- readRDS("./Objects/rarefaction_list_phyloseq_precss.rds")
# Re-run all analayses and save plots as comparisons
# DNA PCoA ----------------------------------------------------------------------------------------------------------

dna.pcoa <- llply(rar.ls, function(x){
  # subset only DNA samples
  dna <- subset_samples(x, dna_type == "DNA")
  
  # remove ASVs that do not appear in this dataset
  dna <- prune_taxa(taxa_sums(dna) != 0, dna)
  dna <- prune_samples(sample_sums(dna) != 0, dna) # remove samples with no reads
  
  # extract ASV matrix
  pb.mat <- otu_mat(dna) # convert phyloseq obj to matrix
  pb.mat <- decostand(pb.mat, "hellinger")
  
  # PCoA with Bray-Curtis
  pb.bray <- vegdist(pb.mat, method = "bray")
  is.euclid(pb.bray) # FALSE
  pb.bray <- sqrt(pb.bray) # make Euclidean
  is.euclid(pb.bray) # TRUE
  # we need euclidean distance for the PERMANOVA later on
  
  ncol(pb.mat) # 16322 OTUs
  nrow(pb.mat) # 389 Samples
  
  # make PCoA
  pb.bray.pcoa <- ape::pcoa(pb.bray)
  # extract scores and variance explained
  dna.pcoa <- plot_pcoa(pb.bray.pcoa, physeq = dna,  colours = colvec, output = T)
  
  # main PCoA
  p <- dna.pcoa$plot + theme_cust(base_theme = "pubr",
                                  border = T) + guides(colour = "none")
  return(p)
})

dna.comp <- ggarrange(dna.pcoa[[6]] + labs(title = "CSS"),
          dna.pcoa[[1]] + labs(title = paste("Rarefaction: ", 
                                             sapply(strsplit(names(dna.pcoa[1]), split = "_"), "[",2))),
          dna.pcoa[[2]] + labs(title = paste("Rarefaction: ", 
                                             sapply(strsplit(names(dna.pcoa[2]), split = "_"), "[",2))),
          dna.pcoa[[3]] + labs(title = paste("Rarefaction: ", 
                                             sapply(strsplit(names(dna.pcoa[3]), split = "_"), "[",2))),
          #dna.pcoa[[4]] + labs(title = paste("Rarefaction: ", 
          #                                  sapply(strsplit(names(dna.pcoa[4]), split = "_"), "[",2))),
          #dna.pcoa[[5]] + labs(title = paste("Rarefaction: ", 
          #                                  sapply(strsplit(names(dna.pcoa[5]), split = "_"), "[",2))),
          ncol = 2, nrow = 2, common.legend = T, align = "hv")

# save
ggsave(paste0("./Figures/Final/Figure2_rarcomp.png"), dna.comp,
       width = 20, height = 20, unit = "cm", dpi = 300)

# DNA-RNA PCoA -----------------------------------------------------------------------------------------

dnarna.pcoa <- llply(rar.ls, function(x){
  # subset only DNA samples
  dna <- subset_samples(x, dna_type == "DNA")
  
  # remove ASVs that do not appear in this dataset
  dna <- prune_taxa(taxa_sums(dna) != 0, dna)
  dna <- prune_samples(sample_sums(dna) != 0, dna) # remove samples with no reads
  # remove OTUs that only appear in RNA
  x <- prune_taxa(taxa_names(x) %in% taxa_names(dna), x)
  
  # removing of OTUs that only appear in DNA did not have a difference on results
  
  # extract species table with species in columns
  pb.mat <- otu_mat(x)
  
  # PCoA with Bray-Curtis
  pb.bray <- vegdist(pb.mat, method = "bray")
  is.euclid(pb.bray) # FALSE
  pb.bray <- sqrt(pb.bray) # make Euclidean
  is.euclid(pb.bray)
  # make PCoA
  pb.bray.pcoa <- ape::pcoa(pb.bray)
  
  # plot with custom function
  all.pcoa <- plot_pcoa(pb.bray.pcoa, physeq = x, axes = c(1:3), colours = colvec, output = T)
  # plot second and third axes
  pcoa.23 <- plot_pcoa(pb.bray.pcoa, physeq = x, plot.axes = c(3,2), colours = colvec, output = T)
  
  find_hull <- function(x, axes){x[chull(x[,paste0("Axis.",axes[1])], x[,paste0("Axis.",axes[2])]),]}
  hulls <- ddply(all.pcoa$df, "dna_type", find_hull, axes = c(1,2))
  
  pcoa.plot <- ggplot(all.pcoa$df, aes(x =Axis.1, 
                                       y = Axis.2)) +
    theme_cust(base_theme = "pubr",
               border = T) +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    scale_fill_manual(values = colvec[names(colvec) %in% as.character(levels(all.pcoa$df$sample.type.year))],
                      name = "Habitat Type") +
    geom_point(aes(fill = sample.type.year, shape = Season), size=2.5) + #colour = dna_type,
    geom_polygon(data = hulls, alpha = 0, aes(linetype = dna_type), fill = "white", colour = "black") +
    scale_linetype_manual(values = c("solid","dotted"), name = "Nucleic Acid Type") +
    coord_fixed(1) + # ensure aspect ratio
    scale_shape_manual(values = c(21,23,25)) +
    #scale_size_manual(values = c(2.5, 2.6), name = "Nucleic Acid \nType") +
    labs(x = paste0("PC1 (", round(all.pcoa$var[1,"var"],
                                   digits = 1),"%)"), 
         y = paste0("PC2 (", round(all.pcoa$var[2,"var"],
                                   digits = 1),"%)")) +
    guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
           fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)),
           size = FALSE)
  
  # 2nd dimension
  hulls <- ddply(pcoa.23$df, "dna_type", find_hull, axes = c(3,2))
  
  pcoa.plot2 <- ggplot(pcoa.23$df, aes(x =Axis.3, 
                                       y = Axis.2)) +
    theme_cust(base_theme = "pubr",
               border = T) +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
    scale_fill_manual(values = colvec[names(colvec) %in% as.character(levels(pcoa.23$df$sample.type.year))],
                      name = "Habitat Type") +
    geom_point(aes(fill = sample.type.year, shape = Season), size=2.5) + #colour = dna_type,
    scale_linetype_manual(values = c("solid","dotted"), name = "Nucleic Acid Type") +
    geom_polygon(data = hulls, alpha = 0, aes(linetype = dna_type), fill = "white", colour = "black") +
    coord_fixed(1) + # ensure aspect ratio
    scale_shape_manual(values = c(21,23,25)) +
    
    #scale_size_manual(values = c(2.5, 2.6), name = "Nucleic Acid \nType") +
    labs(x = paste0("PC3 (", round(pcoa.23$var[1,"var"],
                                   digits = 1),"%)"), 
         y = paste0("PC2 (", round(pcoa.23$var[2,"var"],
                                   digits = 1),"%)")) +
    guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
           fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)),
           size = FALSE)
  p<-list(pcoa.plot, pcoa.plot2)
  
  # p <- ggarrange(pcoa.plot, pcoa.plot2, ncol = 2,
  #                 align = "hv", labels = "auto", common.legend = T, # legend.grob = legs,
  #                 legend = "right")
  return(p)
})

dnarna.comp <- ggarrange(dnarna.pcoa[[6]][[1]] + labs(title = "CSS"),
                      dnarna.pcoa[[6]][[2]],
                      dnarna.pcoa[[1]][[1]] + labs(title = paste("Rarefaction: ", 
                                                         sapply(strsplit(names(dnarna.pcoa[1]), split = "_"), "[",2))),
                      dnarna.pcoa[[1]][[2]],
                      dnarna.pcoa[[2]][[1]] + labs(title = paste("Rarefaction: ", 
                                                         sapply(strsplit(names(dnarna.pcoa[2]), split = "_"), "[",2))),
                      dnarna.pcoa[[2]][[2]],
                      dnarna.pcoa[[3]][[1]] + labs(title = paste("Rarefaction: ", 
                                                         sapply(strsplit(names(dnarna.pcoa[3]), split = "_"), "[",2))),
                      dnarna.pcoa[[3]][[2]],
                      #dnarna.pcoa[[4]] + labs(title = paste("Rarefaction: ", 
                      #                                  sapply(strsplit(names(dnarna.pcoa[4]), split = "_"), "[",2))),
                      #dnarna.pcoa[[5]] + labs(title = paste("Rarefaction: ", 
                      #                                  sapply(strsplit(names(dnarna.pcoa[5]), split = "_"), "[",2))),
                      ncol = 2, nrow = 4, common.legend = T, align = "hv", legend = "right")


# save
ggsave(paste0("./Figures/Final/Figure3_rarcomp.png"), dnarna.comp,
       width = 20, height = 40, unit = "cm", dpi = 300)


# Delta distance -------------------------------------------------------------------------------------------------------

dist.comp <- llply(rar.ls, function(x){
  
  pb.mat <- otu_mat(x)
  
  # PCoA with Bray-Curtis
  pb.bray <- vegdist(pb.mat, method = "bray")
  pb.bray <- sqrt(pb.bray) # make Euclidean
  # make PCoA
  pb.bray.pcoa <- ape::pcoa(pb.bray)
  # Calculate incidence based dissimilarity
  # PCoA with Sorensen, incidence based equivalent of Bray Curtis
  pb.soren <- vegdist(pb.mat, method = "bray", binary = T)
  pb.soren <- sqrt(pb.soren) # make Euclidean
  # make PCoA
  pb.soren.pcoa <- ape::pcoa(pb.soren)
  
  # Extract scores and meta data
  pb.scores <- rbind(data.frame(Sample = as.character(row.names(pb.bray.pcoa$vectors)), Metric = "Bray",
                                pb.bray.pcoa$vectors, stringsAsFactors = F),
                     data.frame(Sample = as.character(row.names(pb.soren.pcoa$vectors)), Metric = "Sorensen",
                                pb.soren.pcoa$vectors, stringsAsFactors = F))# get first three axes
  # merge with a selection of meta data
  meta <- data.frame(Sample = as.character(row.names(sample_df(x))),
                     sample_df(x) %>% dplyr::select(sample.type.year, Season, year, 
                                                     dna_type, distance.from.mouth, dr_match_name), 
                     stringsAsFactors = F)
  data <- merge(pb.scores, meta, by = "Sample")
  data$Sample <- as.character(data$Sample)
  
  setDT(data)
  setcolorder(data, c("Sample","Metric",
                      "sample.type.year","Season","year", "dna_type","distance.from.mouth","dr_match_name",
                      colnames(data)[!(colnames(data) %in% c("Sample","Metric",
                                                             "sample.type.year","Season","year", "dna_type","distance.from.mouth","dr_match_name"))]))
  
  
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
  
  # merge with a selection of meta data
  meta <- data.frame(Sample = as.character(row.names(sample_df(x))),
                     sample_df(x) %>% dplyr::select(sample.type.year, Season, year, 
                                                     dna_type, distance.from.mouth, dr_match_name), 
                     stringsAsFactors = F)
  bray.df <- setDT(merge(meta, bray.df, by = "Sample"))
  soren.df <- setDT(merge(meta, soren.df, by = "Sample"))
  
  # melt and combine bray-curtis and sorensen results into one data frame
  temp.75 <- rbind(melt.data.table(bray.df, id.vars = c("dr_match_name","dna_type", "Metric"), measure.vars = patterns("^Axis."),
                                   variable.name = "Axis", value.name = "Coordinates"),
                   melt.data.table(soren.df, id.vars = c("dr_match_name","dna_type", "Metric"), measure.vars = patterns("^Axis."),
                                   variable.name = "Axis", value.name = "Coordinates"))
  setDT(temp.75)
  
  temp.75 <- dcast(temp.75, dr_match_name + Axis + Metric ~ dna_type, value.var = c("Coordinates"))
  # remove NAs
  temp.75 <- na.omit(temp.75)
  temp.75[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
  temp.75 <- temp.75[, .(sum.dist = sum(pnt.dist)), by = .(Metric, dr_match_name)] # sum temp axes
  temp.75 <- temp.75[, dist := sqrt(sum.dist)] # take sqrt
  
  # 75% variance
  dist.75 <- temp.75[data, c("sample.type.year",
                             "year", "Season") := list(i.sample.type.year,
                                                       i.year, i.Season), on = .(dr_match_name)]
  
  # Calculate delta of the two metrics
  diff.df <- dcast(dist.75, dr_match_name + sample.type.year + Season ~ Metric, value.var = "dist")
  diff.df <- diff.df[, delta := Bray - Sorensen]
  
  diff.df <- melt(diff.df, id.vars = c("dr_match_name","sample.type.year","Season"),
                  measure.vars = c("delta","Sorensen","Bray"),
                  variable.name = "Metric", value.name = "dist")
  sum.delta <- diff.df[, .(mean =  mean(dist, na.rm = T),
                           conf.int = conf.int(dist),
                           stdev = sd(dist, na.rm = T),
                           n = .N),
                       by = .(sample.type.year, Season, Metric)]
  
  # add new column to split plot into main and side panel
  sum.ls <- lapply(list(sum.delta), function(x){
    x[, panels := "main"]
    x[sample.type.year == "Tributary" |
        sample.type.year == "Lake" |
        sample.type.year == "Riverine \nLakes" |
        sample.type.year == "Sediment", panels := "side"]
    x[, sample.type.year := factor(sample.type.year, levels = c("Soil","Sediment",
                                                                "Soilwater",
                                                                "Stream", "Tributary",
                                                                "Riverine \nLakes", "Headwater \nPonds", "Lake",
                                                                "Upriver",
                                                                "Reservoirs","Downriver",
                                                                "Estuary"))]
  })
  
  ##########################
  # Lollipop plots of distance
  temp <- sum.ls[[1]][!is.nan(sum.ls[[1]]$mean),]
  # Rename level for plotting
  levels(temp$sample.type.year)[levels(temp$sample.type.year) == "Riverine \nLakes"] <- "Riv. Lakes"
  names(colvec)[7] <- "Riv. Lakes"
  
  plot.df <- dcast(temp, sample.type.year + Season + panels ~ Metric, value.var = "mean")
  
  temp$Metric <- factor(temp$Metric, levels = c('delta','Bray','Sorensen'), labels = c("Delta","Abundance",
                                                                                       "Incidence"))
  
  
  # Make a fake plot with all the legend units, habitat type, season and metric
  legs <- get_legend(ggplot(temp[!temp$Metric == "Delta",], aes(x = sample.type.year)) +
                       theme_cust("pubr") +
                       geom_point(aes(y=mean, colour = sample.type.year),
                                  size = 3) +
                       scale_colour_manual(values =  colvec[names(colvec) %in% temp$sample.type.year], name = "Habitat type") +
                       geom_point(aes(y=mean, fill = Metric, shape = Season), size = 3, colour = 'black') +
                       scale_fill_manual(values = c("gray40","white"), name = "Metric") +
                       scale_shape_manual(values = c(21, 23, 25)) +
                       guides(shape = guide_legend(order = 3, override.aes=list(size = 2)),
                              colour = guide_legend(order = 2,
                                                    override.aes=list(size = 2, shape = 21,
                                                                      fill = colvec[names(colvec) %in% temp$sample.type.year], 
                                                                      colour = "black")),
                              fill = guide_legend(order = 1, override.aes=list(shape = 21, size = 2))))
  
  
  m <- ggplot(plot.df[panels == "main",], aes(x = sample.type.year)) +
    theme_cust(base_theme = "pubr") +
    facet_grid(.~Season, scales = "free_x", space = "free") +
    geom_linerange(aes(ymin = Sorensen, ymax = Bray,
                       colour = sample.type.year, group = Season),
                   colour = "gray30", linetype = "solid", position = position_dodge(0.7)) +
    geom_jitter(aes(y = Sorensen, colour = sample.type.year,
                    shape = Season), colour = "black",
                fill = "white", position = position_dodge(0.7), size = 3) +
    scale_fill_manual(values =  colvec[names(colvec) %in% temp$sample.type.year], name = "Habitat type") +
    scale_shape_manual(values = c(21, 23, 25)) +
    geom_jitter(aes(y = Bray, fill = sample.type.year,
                    shape = Season), position = position_dodge(0.7), size = 3) +
    scale_colour_manual(values = colvec) +
    labs(y= "Mean DNA-RNA distance", x = "Habitat Type") +
    lims(y = c(min(plot.df$Sorensen, na.rm = T),
               max(plot.df$Bray, na.rm = T))) +
    theme(
      axis.text.x = element_blank(),
      axis.text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title = element_text(size = 10),
      strip.background = element_rect(fill = "gray20"),
      strip.text = element_text(colour = "white", size = 10)) +
    guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
           fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
  
  
  s <- ggplot(plot.df[panels == "side",], aes(x = sample.type.year)) +
    theme_cust(base_theme = "pubr") +
    facet_grid(.~"Meta-community") +
    geom_linerange(aes(ymin = Sorensen, ymax = Bray,
                       colour = sample.type.year, group = Season),
                   colour = "gray30", linetype = "solid", position = position_dodge(0.7)) +
    geom_jitter(aes(y = Sorensen, colour = sample.type.year,
                    shape = Season), colour = "gray20",
                fill = "white", position = position_dodge(0.7), size = 3) +
    scale_fill_manual(values =  colvec[names(colvec) %in% temp$sample.type.year]) +
    scale_shape_manual(values = c(21, 23, 25)) +
    geom_jitter(aes(y = Bray, fill = sample.type.year,
                    shape = Season), position = position_dodge(0.7), size = 3) +
    scale_colour_manual(values = colvec) +
    labs(y= "DNA-RNA distance", x = "Habitat Type") +
    lims(y = c(min(plot.df$Sorensen, na.rm = T),
               max(plot.df$Bray, na.rm = T))) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text = element_text(size = 8),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_rect(fill = "gray20"),
      strip.text = element_text(colour = "white", size = 10)
    ) + guides(shape = "none", fill = "none")
  
  # make a numeric sample type column to add smooth line
  plot.df[sample.type.year == "Soil", num.hab := 1]
  plot.df[sample.type.year == "Soilwater", num.hab := 2]
  plot.df[sample.type.year == "Stream", num.hab := 3]
  plot.df[sample.type.year == "Upriver", num.hab := 4]
  plot.df[sample.type.year == "Reservoirs", num.hab := 5]
  plot.df[sample.type.year == "Downriver", num.hab := 6]
  plot.df[sample.type.year == "Estuary", num.hab := 7]
  
  m.down <- ggplot(plot.df[panels == "main",], ) +
    theme_cust(base_theme = "pubr") +
    facet_grid(.~Season, scales = "free_x", space = "free") +
    geom_line(aes(x = num.hab,
                  y = delta),
              stat = "smooth", method = "loess", span = 0.7, se = F,
              alpha = 0.5, colour= "gray10", size = 0.8) +
    geom_jitter(aes(x = num.hab, y = delta, fill = sample.type.year,
                    shape = Season), 
                position = position_dodge(0.7), size = 3) +
    scale_fill_manual(values =  colvec[names(colvec) %in% temp$sample.type.year], name = "Habitat type") +
    scale_linetype_manual(values = c("solid","dashed","dotted")) +
    scale_shape_manual(values = c(21, 23, 25)) +
    labs(y= expression(paste(Delta, " Distances")), x = "Habitat Type") +
    scale_x_continuous(breaks = 1:length(unique(as.character(plot.df[panels == "main",]$sample.type.year))),
                       labels = unique(as.character(plot.df[panels == "main",]$sample.type.year))) +
    theme(axis.text.x = element_text(angle = 45, hjust =1),
          axis.text = element_text(size = 8),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    guides(fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2))) +
    scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(min(plot.df$delta, na.rm = T)-0.03,
                                                         max(plot.df$delta, na.rm = T)+0.01))
  
  s.down <-ggplot(plot.df[panels == "side",], ) +
    theme_cust(base_theme = "pubr") +
    facet_grid(.~"Meta-community") +
    geom_jitter(aes(x = sample.type.year, y = delta, fill = sample.type.year,
                    shape = Season), 
                position = position_dodge(0.7), size = 3) +
    scale_fill_manual(values = colvec[names(colvec) %in% temp$sample.type.year], name = "Habitat type") +
    scale_shape_manual(values = c(21, 23, 25)) +
    labs(y= expression(paste(Delta, " Distances")), x = "Habitat Type") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 8),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    scale_y_continuous(breaks = c(0,0.2,0.4), limits = c(min(plot.df$delta, na.rm = T)-0.03,
                                                         max(plot.df$delta, na.rm = T)+0.01)) +
    guides(fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
  
  down <- ggarrange(m.down,s.down, ncol = 2, widths = c(3,1), legend = "none")
  
  up <- ggarrange(m,s, ncol = 2, nrow = 1,  widths = c(3,1), legend = "none")
  
  out <- list(up, down, legs)
  return(out)
  
})

delta.p <- ggarrange(annotate_figure(dist.comp[[6]][[1]], top = "CSS"),
                     annotate_figure(dist.comp[[1]][[1]], top = paste("Rarefaction: ", 
                                                                    sapply(strsplit(names(dist.comp[1]), split = "_"), "[",2))),
                     dist.comp[[6]][[2]], 
                     dist.comp[[1]][[2]], 
                     annotate_figure(dist.comp[[2]][[1]], top = paste("Rarefaction: ", 
                                                                    sapply(strsplit(names(dist.comp[2]), split = "_"), "[",2))),
                    annotate_figure(dist.comp[[3]][[1]], top = paste("Rarefaction: ", 
                                                                    sapply(strsplit(names(dist.comp[3]), split = "_"), "[",2))),
                         dist.comp[[2]][[2]],
                         dist.comp[[3]][[2]],
                         #dnarna.pcoa[[4]] + labs(title = paste("Rarefaction: ", 
                         #                                  sapply(strsplit(names(dnarna.pcoa[4]), split = "_"), "[",2))),
                         #dnarna.pcoa[[5]] + labs(title = paste("Rarefaction: ", 
                         #                                  sapply(strsplit(names(dnarna.pcoa[5]), split = "_"), "[",2))),
                         ncol = 2, nrow = 4, common.legend = T, legend.grob = dist.comp[[6]][[3]],
                         align = "hv", legend = "right", heights = c(2,1,2,1))

# save
ggsave(paste0("./Figures/Final/Figure4_rarcomp.png"), delta.p,
       width = 33, height = 30, unit = "cm", dpi = 300)

# Proportion reactive --------------------------------------------------------------------------------------------------------------------


reac.unreactive <- llply(rar.ls, function(x){
  
  # melt community matrix for tidy data set
  commat <- melt.data.table(
    setDT(as.data.frame(otu_mat(x)), keep.rownames = "Sample"),
    id.vars = "Sample",
    measure.vars = patterns("^OTU_"),
    variable.name = "OTU",
    value.name = "reads"
  )
  
  # add meta variables
  commat[setDT(sample_df(x), keep.rownames = "Sample"), 
         c("dr_match_name", "dna_type","sample.type.year","Season") := list(i.dr_match_name, i.dna_type,
                                                                            i.sample.type.year, i.Season),
         on = .(Sample)]
  
  # cast so that we can calculate difference between DNA and RNA of each OTU
  temp <- dcast(commat, dr_match_name + sample.type.year + OTU ~ dna_type, value.var = c("reads"))
  temp[, diff := DNA - RNA]
  temp <- temp[!is.na(diff),]
  
  # Calculate percentage of activity of taxa first found in...
  # Soil > Soilwater > Stream > Upriver > Reservoir > Downriver > Estuary
  # Focus only on the direct continuum
  # where were the OTUs first observed (based on DNA)
  temp <- temp[(sample.type.year == "Soil" |
                  sample.type.year == "Soilwater" |
                  sample.type.year == "Stream" |
                  sample.type.year == "Upriver" |
                  sample.type.year == "Reservoirs" |
                  sample.type.year == "Downriver" |
                  sample.type.year == "Estuary"),]
  
  # all OTUs in...
  soil <- unique(as.character(temp[DNA > 0 & sample.type.year == "Soil",]$OTU))
  soilwater <- unique(as.character(temp[DNA > 0 & sample.type.year == "Soilwater",]$OTU))
  stream <- unique(as.character(temp[DNA > 0 & sample.type.year == "Stream",]$OTU))
  upriver <- unique(as.character(temp[DNA > 0 & sample.type.year == "Upriver",]$OTU))
  reservoir <- unique(as.character(temp[DNA > 0 & sample.type.year == "Reservoirs",]$OTU))
  downriver <- unique(as.character(temp[DNA > 0 & sample.type.year == "Downriver",]$OTU))
  estuary <- unique(as.character(temp[DNA > 0 & sample.type.year == "Estuary",]$OTU))
  
  # now, remove OTUs from vectors that were found before
  soilwater <- soilwater[!soilwater %in% soil]
  stream <- stream[!stream %in% unique(c(soil,soilwater))]
  upriver <- upriver[!upriver %in% unique(c(soil,soilwater,stream))]
  reservoir <- reservoir[!reservoir %in% unique(c(soil,soilwater,stream,upriver))]
  downriver <- downriver[!downriver %in% unique(c(soil,soilwater,stream, upriver, reservoir))]
  estuary <- estuary[!estuary %in% unique(c(soil,soilwater,stream, upriver, reservoir, downriver))]
  
  # sanity check
  all.otus <- unique(as.character(temp[DNA > 0 & (sample.type.year == "Soil" |
                                                    sample.type.year == "Soilwater" |
                                                    sample.type.year == "Stream" |
                                                    sample.type.year == "Upriver" |
                                                    sample.type.year == "Reservoirs" |
                                                    sample.type.year == "Downriver" |
                                                    sample.type.year == "Estuary"),]$OTU))
  length(c(soil, soilwater, stream, upriver, reservoir, downriver, estuary)) == length(all.otus) # need to be TRUE
  
  # 12795 of original 16368
  # some OTUs are unique to ecosystems outside of the continuum
  
  # make OTU character
  temp[, OTU := as.character(OTU)]
  # create new column with first observed category for each OTU
  temp[OTU %in% soil, first.obs := "Soil"]
  temp[OTU %in% soilwater, first.obs := "Soilwater"]
  temp[OTU %in% stream, first.obs := "Stream"]
  temp[OTU %in% upriver, first.obs := "Upriver"]
  temp[OTU %in% reservoir, first.obs := "Reservoirs"]
  temp[OTU %in% downriver, first.obs := "Downriver"]
  temp[OTU %in% estuary, first.obs := "Estuary"]
  
  # make factor
  temp[, first.obs := factor(first.obs, levels = c("Soil", "Soilwater", "Stream",
                                                   "Upriver", "Reservoirs", "Downriver",
                                                   "Estuary"))]
  # remove all OTUs that are found outside the direct continuum
  first.df <- temp[!is.na(first.obs),]
  # all OTUs are characterized
  any(!(first.df$OTU %in% c(soil,soilwater,stream, upriver,reservoir,downriver,estuary)) == T) # need to be FALSE
  
  # sanity check
  first.df[is.na(DNA) & RNA >0,]
  first.df[DNA == 0 & RNA >0,]
  
  # Calculate ratio as a comarison
  first.df[, ratio := DNA / RNA]
  
  # merge with some meta data
  first.df[setDT(sample_df(x), keep.rownames = "Sample"),
           c("sample.type.year","Season") := list(i.sample.type.year,
                                                          i.Season), on = .(dr_match_name)]
  
  # remove any NAs in the dataset
  first.df[is.na(DNA),]
  first.df[is.na(RNA),]
  
  setDT(commat)
  # calculate means by sample type
  means <- commat[, .(mean.reads = mean(reads, na.rm = T),
                      sd.css = sd(reads, na.rm = T)), by = .(sample.type.year, dna_type, OTU)]
  # order the abundances to make ranks
  means <- means[mean.reads != 0,] # remove all 0 observations
  means <- means[order(mean.reads, decreasing = T)]
  means[, rank.abun := 1:.N, by = .(dna_type, sample.type.year)]
  means[, log.mean := log1p(mean.reads), by = .(dna_type, sample.type.year)]
  
  # smooth and get derivative
  classif.thres <- ddply(means, .(dna_type, sample.type.year), function(x){
    spl <- smooth.spline(x$rank.abun, x$log.mean, spar = 0.7)
    #pred <- predict(spl)
    #first <- predict(spl, deriv = 1) # first derivative
    sec <- predict(spl, deriv = 2) # second derivative
    setDT(x)
    x[rank.abun <= x$rank.abun[localMaxima(sec$y)[1]], ab.group := "Abundant", by = .(OTU)]
    x[rank.abun > x$rank.abun[localMaxima(sec$y)[1]] &
        rank.abun <= x$rank.abun[localMaxima(sec$y)[2]], ab.group := "Medium", by = .(OTU)]
    x[rank.abun > x$rank.abun[localMaxima(sec$y)[2]], ab.group := "Rare", by = .(OTU)]
    
    out <- x[, .(max = max(mean.reads),
                 min = min(mean.reads)), by = .(ab.group)]
    return(out)
  }, .parallel = T)
  
  setDT(classif.thres)
  thres.df <- classif.thres[, .(mean.max = mean(max),
                    sd.max = sd(max),
                    mean.min = mean(min),
                    sd.min = sd(min)), by = .(dna_type, ab.group)]
  a.thres <- round(thres.df[ab.group == "Abundant" & dna_type == "DNA",]$mean.max, digits = 0)
  m.thres <- round(thres.df[ab.group == "Medium" & dna_type == "DNA",]$mean.max, digits = 0)
  r.thres <- round(thres.df[ab.group == "Rare" & dna_type == "DNA",]$mean.max, digits = 0)
  
  # Add abundance group to each observation
  first.df[DNA >= a.thres, abg.dna := "Abundant"]
  first.df[DNA < a.thres & DNA >= m.thres, abg.dna := "Moderate"]
  first.df[DNA < m.thres, abg.dna := "Rare"]

  # Extract reactive fraction: RNA > 0
  rna.df <- first.df[RNA > 0, 
                     .(sum.rna.byotu = sum(RNA, na.rm = T),
                       mean.ratio = mean(ratio, na.rm = T)),
                     by = .(OTU, sample.type.year, abg.dna, first.obs, Season)]
  
  # Extract all (DNA > 0 and RNA >0 + DNA > 0 and RNA = 0)
  dna.df <- first.df[DNA > 0, 
                     .(sum.dna.byotu = sum(DNA, na.rm = T)),
                     by = .(OTU, sample.type.year, abg.dna, first.obs, Season)]
  mean.df <- rna.df[dna.df, , on = .(OTU, sample.type.year, abg.dna, first.obs, Season)]
  
  # Calculate the number of reads only in reactive fraction
  ac.all <- first.df[RNA > 0,
                     .(otu.n.all.rna = .N,
                       sum.reads.all.rna = sum(RNA, na.rm = T)), by = .(sample.type.year, Season)]  # dr_match_name
  mean.df[ac.all, "sum.reads.all.rna" := i.sum.reads.all.rna, on = .(sample.type.year,Season)]
  
  # Sum of all reads of all observations of DNA (not only DNA >0 and RNA = 0, like before)
  all <- first.df[DNA > 0,
                  .(otu.n.all.dna = .N,
                    sum.reads.all.dna = sum(DNA, na.rm = T)), by = .(sample.type.year, Season)] #dr_match_name
  mean.df[all, "sum.reads.all.dna" := i.sum.reads.all.dna, on = .(sample.type.year,Season)]
  
  # Calculate percentage
  mean.df[, perc.sum.rna.byotu := sum.rna.byotu * 100 / sum.reads.all.rna]
  mean.df[, perc.sum.dna.byotu := sum.dna.byotu * 100 / sum.reads.all.dna]
  
  
  # sanity check
  mean.df[, .(sum.all = sum(perc.sum.rna.byotu, na.rm = T)),
          by = .(sample.type.year, Season)]
  mean.df[, .(sum.all = sum(perc.sum.dna.byotu, na.rm = T)),
          by = .(sample.type.year, Season)] # all 100%
  
  # make quartile bins to represent variation better
  # overall mean simplifies pattern too much. Tried 4 bins first by:
  # 1. < 25th quartile, 2. 25th quartile - Median, 3. Median - 75th quartile, 4. > 75th quartile
  # 1+2 and 3+4 were relatively similar, and hence, decide on only two bins.
  # mean.df[perc.sum.rna.byotu < exp(-5), quan.bin := 1] #  < 50%
  # mean.df[perc.sum.rna.byotu >= exp(-5), quan.bin := 2] # >= 50%
  mean.df[perc.sum.rna.byotu < summary(perc.sum.rna.byotu)[3], quan.bin := 1] #  < 50%
  mean.df[perc.sum.rna.byotu >= summary(perc.sum.rna.byotu)[3], quan.bin := 2] # >= 50%

  # Any OTUs not assigned a quartile bin?
  any(is.na(mean.df$quan.bin))
  # only OTUs that don't have RNA, that's fine
  
  # save by OTU data frame for later
  byotu.df <- mean.df
  
  # summarise by bins
  mean.df <- mean.df[,.(perc.sum.rna.otu = mean(perc.sum.rna.byotu, na.rm = T),
                        perc.sum.dna.otu = mean(perc.sum.dna.byotu, na.rm = T),
                        mean.ratio = mean(mean.ratio, na.rm = T)),
                     by = .(sample.type.year, abg.dna, first.obs, Season, quan.bin)]
  
  mean.mean <- mean.df[quan.bin == 2, .(mean = mean(perc.sum.rna.otu),
                                        sd = sd(perc.sum.rna.otu)), by = .(first.obs)]
  
  lins <- dlply(mean.df[quan.bin == 2,], .(sample.type.year, Season), function(x){
    lm(log(perc.sum.dna.otu) ~ log(perc.sum.rna.otu), data = x)
  })
  
  (slopes.df <- ldply(lins, coef))
  # Prepare for plotting
  # overwrite soil colour, original beige is too hard to see
  new.col <- colvec
  new.col[names(new.col) == "Soil"] <- "gold"
  new.col <- new.col[names(new.col) %in% mean.df$sample.type.year]
  # make df for annotation of Median
  ann_df <- data.frame(perc.sum.dna.otu = -5,
                       perc.sum.rna.otu = as.numeric(summary(log(mean.df$perc.sum.rna.otu))[3]),
                       Season = factor("Summer",levels = c("Spring","Summer","Autumn")),
                       sample.type.year = factor("Soil"))
  # leave panels empty where we don't have any data
  line_df <- data.frame(Season = c("Summer","Spring", "Summer", "Spring", "Summer",
                                   "Spring", "Summer", "Autumn","Spring", "Summer", "Autumn",
                                   "Spring", "Summer", "Autumn",
                                   "Summer"),
                        sample.type.year = c("Soil","Soilwater","Soilwater", "Stream","Stream",
                                             "Upriver","Upriver","Upriver",
                                             "Downriver","Downriver","Downriver",
                                             "Reservoirs","Reservoirs","Reservoirs",
                                             "Estuary"))
  line_df$x <- 0 ; line_df$y <- summary(log(mean.df$perc.sum.rna.otu))[3]
  #line_df$x <- 0 ; line_df$y <- -5
  line_df$Season <- factor(line_df$Season, levels = c("Spring","Summer","Autumn"))
  setDT(line_df)
  line_df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                    "Upriver", "Reservoirs", "Downriver",
                                                                    "Estuary"))]
  # Make factor
  mean.df[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                    "Upriver", "Reservoirs", "Downriver",
                                                                    "Estuary"))]
  
  # Categorize median bins
  mean.df[, quan.bin := factor(quan.bin, levels = c(2,1),
                               labels = c("> Median","< Median"))]
  
  
  perc.cont.dna.rna <- ggplot(mean.df[!is.na(quan.bin),]) +
      theme_cust("pubr")+
      facet_grid(sample.type.year~Season) +
      geom_hline(data = line_df,
                 aes(yintercept = y),
                 linetype = "solid", colour = "gray70") +
      geom_smooth(aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), group = as.character(quan.bin)),
                  method = "lm", se = F, colour = "gray30", size = 0.5, linetype = "solid",
                  data = mean.df[!is.na(quan.bin) & quan.bin == "> 0.0067",]) +
      #ggforce::geom_mark_ellipse(data = subset(mean.df[!is.na(quan.bin) & Season == "Summer" 
      #                                        & sample.type.year == "Soil",], (quan.bin == 2)),
      #                  aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), group = as.character(quan.bin)),
      #                  linetype = "dashed", colour = "gray30", expand = unit(1.5, "mm")) +
      geom_point(aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), colour = first.obs,
                     shape = abg.dna, fill = quan.bin), alpha = 0.9, size =2) +
      geom_text(data = ann_df,
                aes(y = perc.sum.rna.otu +1, x = perc.sum.dna.otu+1.2), label = "Threshold", size =2.5,
                colour = "gray50") +
      scale_colour_manual(values = new.col,
                          name = "First \ndetected in") +
      scale_shape_manual(values = c(21,23,25),
                         name = "Local DNA\nabundance") +
      scale_fill_manual(values = c("white","gray80"),
                        name = "Reactivity\nthreshold (%)") +
      scale_y_continuous(limits = c(min(log(mean.df$perc.sum.rna.otu), na.rm = T) - 1, 
                                    max(log(mean.df$perc.sum.rna.otu), na.rm = T) + 1)) +
      theme(panel.grid = element_blank(),
            legend.position = "right",
            strip.background = element_rect(fill = "gray20"),
            strip.text = element_text(colour = "white", size = 9)) +
      labs(x = "Mean OTU contribution to\nlocal DNA pool (%, log-scale)",
           y = "Mean OTU contribution to\nlocal RNA pool (%, log-scale)") +
      guides(shape = guide_legend(order = 2, override.aes=list(size = 3)),
             fill = guide_legend(order = 1, override.aes=list(shape=21, size = 3)),
             colour = guide_legend(order = 3, override.aes=list(size = 3)))
  
  #
  #stat_ellipse(data = subset(mean.df[!is.na(quan.bin),], (quan.bin == 2)),
  #             aes(x = log(perc.sum.dna.otu), y = log(perc.sum.rna.otu), group = as.character(quan.bin)),
  #             linetype = "dashed", colour = "grey50") +
  
  # removed observations are OTUs with no RNA
  # Annotate second x-axis (panels)
  ann.perc.cont.dna.rna <- annotate_figure(perc.cont.dna.rna + theme(legend.position = "none"),
                                            right = text_grob("Habitat Type", rot = 270, hjust = 0.7))
  # merge with earlier percentage first detected bar plots
  #(first.obs.dna.rna <- ggarrange(frac.first.obs, perc.cont.dna.rna, nrow = 1, ncol = 2, labels = c("","c"), hjust = -10, widths = c(1,2)))
  # Save
  #ggsave("./Figures/Final/first.observed.contribution.dna.rna.png", 
  #       first.obs.dna.rna, width = 25, height = 15, units = "cm", dpi = 300)
  #ggsave("./Figures/Final/first.observed.contribution.dna.rna.tiff", 
  #       first.obs.dna.rna, width = 25, height = 15, units = "cm", dpi = 300)
  
  # Calculate the percentage they occupy in DNA ---------------------------------------------------------------
  
  # Calculate perc of first observed of < Median, > Median, and no RNA OTUs
  subdf <- byotu.df[, .(sum.first.medlow = sum(sum.dna.byotu, na.rm = T)), .(first.obs, Season, sample.type.year, quan.bin)]
  all.sum <- byotu.df[,.(sum.dna = sum(sum.dna.byotu, na.rm = T)), .(Season, sample.type.year, quan.bin)]
  
  subdf <- subdf[all.sum, , on = .(Season, sample.type.year, quan.bin)]
  subdf[, perc.dna := sum.first.medlow * 100 / sum.dna]
  
  # Categorize median bins
  subdf[, quan.bin := factor(quan.bin, levels = c(2,1),
                             labels = c("> Threshold\n(reactive)","< Threshold\n(unreactive)"))]
  # Add a new bin for the non RNA OTUs
  subdf[is.na(quan.bin), quan.bin := "RNA = 0\n(unreactive)"]
  
  # Level sample types
  subdf[, sample.type.year := factor(sample.type.year, levels = c("Soil", "Soilwater", "Stream",
                                                                  "Upriver", "Reservoirs", "Downriver",
                                                                  "Estuary"))]
  un.perc.dna <- ggplot(subdf[quan.bin == "RNA = 0\n(unreactive)",], aes(x = sample.type.year)) +
      theme_cust("pubr") +
      geom_col(aes(y = perc.dna, fill = first.obs), colour = "gray20", size = 0.3) +
      facet_grid(quan.bin~Season, scale = "free_x", space = "free") +
      theme(axis.text.x = element_blank(), #axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
            axis.title.x = element_blank(),
            strip.background = element_rect(fill = "gray20"),
            strip.text = element_text(colour = "white", size = 9)) +
      labs(x = "Habitat type", y = "Reads in local\nDNA pool (%)") +
      scale_fill_manual(values = colvec[names(colvec) %in% subdf$sample.type.year],
                        name = "First observed in") +
      guides(fill = "none")
  
  perc.dna <- ggplot(subdf[!quan.bin == "RNA = 0\n(unreactive)",], aes(x = sample.type.year)) +
      theme_cust("pubr") +
      geom_col(aes(y = perc.dna, fill = first.obs), colour = "gray20", size = 0.3) +
      facet_grid(quan.bin~Season, scale = "free_x", space = "free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
            #axis.title.y = element_blank(),
            strip.background.x = element_blank(),
            strip.background = element_rect(fill = "gray20"),
            strip.text = element_text(colour = "white", size = 9)) +
      labs(x = "Habitat type", y = "Reads in local\nDNA pool (%)") +
      scale_fill_manual(values = colvec[names(colvec) %in% subdf$sample.type.year],
                        name = "First observed in") +
      guides(fill = "none")
  
  p.list <- list(un.perc.dna, perc.dna, ann.perc.cont.dna.rna, get_legend(perc.cont.dna.rna))
  return(p.list)
  
  
})

perc.rarcomp <- 
    ggarrange(
      annotate_figure(
        ggarrange(ggarrange(reac.unreactive[[6]][[1]], reac.unreactive[[6]][[2]], nrow = 2, heights = c(1,2.1), labels = c("a","c")),
                reac.unreactive[[6]][[3]], ncol = 2, labels = c("","b"),
                widths = c(1.4,2), legend = "none"), top = "CSS"),
      annotate_figure(
        ggarrange(ggarrange(reac.unreactive[[1]][[1]], reac.unreactive[[1]][[2]], nrow = 2, heights = c(1,2.1), labels = c("a","c")),
              reac.unreactive[[1]][[3]], ncol = 2, labels = c("","b"),
              widths = c(1.4,2), legend = "none"),
        top = paste("Rarefaction: ", 
                    sapply(strsplit(names(reac.unreactive[1]), split = "_"), "[",2))),
      annotate_figure(
        ggarrange(ggarrange(reac.unreactive[[2]][[1]], reac.unreactive[[2]][[2]], nrow = 2, heights = c(1,2.1), labels = c("a","c")),
                reac.unreactive[[2]][[3]], ncol = 2, labels = c("","b"),
                widths = c(1.4,2), legend = "none"),
        top = paste("Rarefaction: ", 
                    sapply(strsplit(names(reac.unreactive[2]), split = "_"), "[",2))),
      annotate_figure(
        ggarrange(ggarrange(reac.unreactive[[3]][[1]], reac.unreactive[[3]][[2]], nrow = 2, heights = c(1,2.1), labels = c("a","c")),
                reac.unreactive[[3]][[3]], ncol = 2, labels = c("","b"),
                widths = c(1.4,2), legend = "none"),
        top = paste("Rarefaction: ", 
                    sapply(strsplit(names(reac.unreactive[3]), split = "_"), "[",2))),
      nrow = 2, ncol = 2, legend.grob = reac.unreactive[[1]][[4]], legend = "right")



ggsave("./Figures/Final/Figure5_rarcomp.png", 
       perc.rarcomp, width =38, height = 25, units = "cm", dpi = 300)
