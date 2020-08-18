# Apply beta partitioning framework of Schmera et al 2020
# Partitions are:
# Intersection (I) of nestedness and beta = importance of nestedness
# Relative complement (RC) of nestedness in beta

# Weiher Boylen = interpreted as raw measures, gives number of species froming pairwise pattern components (PPCs)
# make calculation independent of total number of species present at the two sites: 
# Jaccard, interpreted as relativised measures expressed as percentage of the number of species in both sites
# Sorensen: gives more weight to species common between compared samples
# interpreted as relativised measures expressed as percentage of the average number of species in both sites

# Input to functions are:
# presence-absence matrix, where rows are sites, columns are species

# extract species table with species in columns
pb.mat <- t(otu_mat(pb))
# convert to presence-absence
pb.mat.pa <- decostand(pb.mat, "pa")

# calculate weiher boylen partitions
wb.dists <- setpart.wb(pb.mat.pa)

library(simba)

# combine all into a melted dataframe
parts.wb <- Reduce(function(x,y, ...) merge(x, y, by = c("Sample.x", "Sample.y"), ...),
                   list(melt.dist(wb.dists[[1]], var.name = "BD"),
                        melt.dist(wb.dists[[2]], var.name = "I"),
                        melt.dist(wb.dists[[3]], var.name = "RC"),
                        melt.dist(sim(pb.mat.pa, method = "weiher"), var.name = "S")))

# calculate jaccard partitions
j.dists <- setpart.j(pb.mat.pa)

parts.j <- Reduce(function(x,y, ...) merge(x, y, by = c("Sample.x", "Sample.y"), ...),
                  list(melt.dist(j.dists[[1]], var.name = "BD"),
                       melt.dist(j.dists[[2]], var.name = "I"),
                       melt.dist(j.dists[[3]], var.name = "RC"),
                       melt.dist(sim(pb.mat.pa, method = "jaccard"), var.name = "S")))

# merge both distances
bd.df <- rbind(parts.wb %>% mutate(Metric = "Weiher Boylen"), 
               parts.j %>% mutate(Metric = "Jaccard"))

# Get meta data to rename DNA and RNA data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(DnaType, Year, Season, sample.type.year), 
                   stringsAsFactors = F)

# add meta data for both sample x and sample y
bd.df <- merge(bd.df, meta, by.x =  "Sample.x", by.y = "Sample")
bd.df <- merge(bd.df, meta %>% select(DnaType, Sample, Year), by.x =  "Sample.y", by.y = "Sample")

# omit all samples of 2015 (no RNA was taken, and sample name strategy changed -> creates duplicates)
bd.df <- bd.df[bd.df$Year.x != 2015,]
bd.df <- bd.df[bd.df$Year.y != 2015,]

# correct a few wrong sample names for matching DNA and RNA
bd.df <- as.data.frame(cbind(apply(bd.df[1:2], 2, function(x){
  x <- str_replace(x, "RO2R52R","RO2.52R")
  x <- str_replace(x, "SWR34R","SW34R")
  x <- str_replace(x, "RO2.36pD","RO2.36D")
  x <- str_replace(x, "RO2.36pR","RO2.36R")
  x <- str_replace(x, "RO2111.60mD", "RO2111.90mD")
  x <- str_replace(x, "RO2.30DPR", "RO2.30R") # two DNA
  x <- str_replace(x, "RO301.HypoR", "RO31.HypoR")
  x <- str_replace(x, "RO301R", "RO31R")
  x <- str_replace(x, "RO304R", "RO34R")
  x <- str_replace(x, "RO307R", "RO37R")
  x <- str_replace(x, "L230R", "L330R")
}), bd.df[3:ncol(bd.df)]))

# remove Ds and Rs to match counterpart samples
bd.df$ID.x[bd.df$DnaType.x == "DNA"] <- str_replace(bd.df$Sample.x[bd.df$DnaType.x == "DNA"], "D$", "")
bd.df$ID.x[bd.df$DnaType.x == "cDNA"] <- str_replace(bd.df$Sample.x[bd.df$DnaType.x == "cDNA"], "R$", "")
bd.df$ID.y[bd.df$DnaType.y == "DNA"] <- str_replace(bd.df$Sample.y[bd.df$DnaType.y == "DNA"], "D$", "")
bd.df$ID.y[bd.df$DnaType.y == "cDNA"] <- str_replace(bd.df$Sample.y[bd.df$DnaType.y == "cDNA"], "R$", "")

# keep all rows where Sample.x and Sample.y are the same
bd.df <- bd.df[bd.df$ID.x == bd.df$ID.y,]
bd.df <- bd.df[!(bd.df$DnaType.x == bd.df$DnaType.y),] # omit all distances between same DnaType

# check duplicates
#bd.df[bd.df$ID.x %in% bd.df[duplicated(bd.df$ID.x),]$ID.x,]

bd.df <- bd.df %>% select(Metric, ID = ID.x, Year = Year.x, Season, sample.type.year,
                          BD, I, RC, S)

# calculate mean dissimilarities for duplicates
sum <- bd.df %>% 
  dplyr::group_by(ID, Metric) %>%
  dplyr::summarise(BD = mean(BD),
                   I = mean(I),
                   RC = mean(RC),
                   S = mean(S),
                   n = n()) %>%
  ungroup()

# combine back with categories
setDT(sum)
bd.dt <- sum[bd.df, c("sample.type.year",
                      "Year", "Season" ) := list(i.sample.type.year,
                                                 i.Year, i.Season), on = .(ID)]

bd.dt$sample.type.year <- factor(bd.dt$sample.type.year, levels = c("Soil","Sediment",
                                                                    "Soilwater","Hyporheicwater", 
                                                                    "Wellwater","Stream", "Tributary",
                                                                    "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                    "Upriver","RO3", "RO2", "RO1","Deep",
                                                                    "Downriver",
                                                                    "Marine"),
                                 labels = c("Soil","Sediment",
                                            "Soilwater","Hyporheicwater", 
                                            "Groundwater","Stream", "Tributary",
                                            "Headwater \nLakes", "Upstream \nPonds", "Lake", "Lake",
                                            "Upriver","RO3","RO2", "RO1","Hypolimnion","Downriver",
                                            "Estuary"))
bd.dt$Season <- factor(bd.dt$Season, levels = c("spring", "summer", "autumn"), 
                       labels = c("Spring", "Summer","Autumn"))


ggtern(bd.dt[Metric == "Jaccard",], aes(x = I, y = RC, z = S)) +
  theme_bw() +
  geom_point(aes(fill = sample.type.year), shape = 21) +
  theme_showarrows() +
  facet_grid(.~Season) +
  scale_fill_manual(values = colvec)

# add new column to split plot into main and side panel
bd.dt[, panels := "main"]
bd.dt[sample.type.year == "Tributary" |
        sample.type.year == "Lake" |
        sample.type.year == "HeadwaterLakes" |
        sample.type.year == "Sediment", panels := "side"]

# calculate confidence interval and means of sample type and season combinations
plot.df <- bd.dt[, .(BD =  mean(BD, na.rm = T),
                     BD.conf = conf.int(BD),
                     I = mean(I, na.rm = T),
                     I.conf = conf.int(I),
                     RC = mean(RC, na.rm = T),
                     RC.conf = conf.int(RC)), 
                 by = .(Metric, sample.type.year, Season, panels)]




wb <- plot.df[Metric == "Weiher Boylen",]
wb <- plot.df[Metric == "Jaccard",]


wb.melt <- merge(melt(wb, id.vars = c("sample.type.year", "Season", "panels"),
                      measure.vars = c("BD", "I", "RC"),
                      variable.name = "PPC", value.name = "mean"),
                 melt(wb, id.vars = c("sample.type.year", "Season", "panels"),
                      measure.vars = c("BD.conf", "I.conf", "RC.conf"),
                      variable.name = "PPC", value.name = "conf.int")[, PPC := str_replace(.SD$PPC, ".conf", "")],
                 by = c("sample.type.year", "Season", "panels", "PPC"))

# Plot Bray
(
  main.b <-
    ggplot(wb.melt[panels == "main", ], 
           aes(x = sample.type.year, y = mean, fill = PPC
           )) +
    theme_cust(base_theme = "pubr")
  +
    geom_errorbar(aes(ymin = mean - conf.int, ymax = mean + conf.int, colour = Season),
                  position = position_dodge(0.7), width = 0) +
    geom_jitter(aes(fill = PPC), shape = 21, 
                position = position_dodge(0.7), size = 2.5) +
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
    labs(x = "Habitat type", y = "BD Partition") +
    facet_grid(Season~.) +
    #lims(y = c(0.25, 1)) +
    #facet_grid(.~Year, scales = "free") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title = element_text(size = 10)
    )
)

# side panel
(
  side.b <-
    ggplot(plot.df[panels == "side"& Metric == "Bray", ], 
           aes(x = sample.type.year, y = mean, fill = Season
           )) +
    theme_cust(base_theme = "pubr") +
    geom_errorbar(aes(ymin = mean - conf.int, ymax = mean + conf.int, colour = Season),
                  position = position_dodge(0.7), width = 0) +
    geom_jitter(aes(fill = Season), shape = 21, 
                position = position_dodge(0.7), size = 2.5) +
    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
    scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
    labs(x = "Sample type", y = "Bray-Curtis dissimilarity") +
    #lims(y = c(0.25, 1)) +
    #facet_grid(.~Year, scales = "free") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 8),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
)