# Customly written functions to be used

#####################
## Data management ##
#####################

select_newest <- function(path, file.pattern, by = NULL) {
  # some error dealing
  if (length(path) > 1) {
    stop("'path' needs to be a single string.")
  } else if(length(path) == 0L){
    stop("'path' is empty.")
  }
  if (length(file.pattern) > 1) {
    stop("'file.pattern' needs to be a single string.")
  } else if(length(file.pattern) == 0L){
    stop("'file.pattern' is empty.")
  }
  
  # read in files in directory
  files <- list.files(path, pattern = file.pattern)
  
  if(length(files) == 0L){
    stop("No file that matches pattern exists.")
  }
  
  if(!is.null(by)){
    # do we have several files per object? -> take newest version
    for (h in 1:length(by)) {
      versions <- files[grepl(by[h], files)]
      versions.sans <-tools::file_path_sans_ext(files)
      if (length(versions) > 1) {
        find.date <- dplyr::bind_rows(lapply(stringr::str_split(versions.sans,"_"), function(x){
          col.names <- paste0("V", seq(length = length(x)))
          as.data.frame.list(x, col.names = col.names, stringsAsFactors = F)
        }))
        dates <- as.Date(unlist(Filter(Negate(is.null), apply(find.date, 2, function(x){
          is.date <- try(as.Date(x), silent = T)
          if(inherits(is.date, "try-error")){
            x <- NULL
          }
          return(x)
        }))))
        
        newest <- max(dates)
        
        files <-
          c(versions[grepl(newest, versions)], files[!grepl(by[h], files)])
      }
    }
  } else {
    versions <- files
    versions.sans <-tools::file_path_sans_ext(files)
    find.date <- dplyr::bind_rows(lapply(stringr::str_split(versions.sans,"_"), function(x){
      col.names <- paste0("V", seq(length = length(x)))
      as.data.frame.list(x, col.names = col.names, stringsAsFactors = F)
    }))
    dates <- as.Date(unlist(Filter(Negate(is.null), apply(find.date, 2, function(x){
      is.date <- try(as.Date(x), silent = T)
      if(inherits(is.date, "try-error")){
        x <- NULL
      }
      return(x)
    }))))
    
    newest <- max(dates)
    files <- versions[grepl(newest, versions)]
  }
  return(files)
}

########################
## phyloseq wrangling ##
########################
## Functions to index in phloseq objects

otu_mat <- function(ps) as(otu_table(ps), "matrix")
tax_mat <- function(ps) as(tax_table(ps), "matrix")
sample_df <- function(ps) as(sample_data(ps), "data.frame")

## Function to import physeq into metagenomeSeq, modified
## (original gave errors, matrix conversion issues)
physeq_to_metagenomeSeq_mod <- function (physeq, ...) 
{
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  countData = round(as(otu_table(physeq), "matrix"), digits = 0)
  if (!is.null(sample_data(physeq, FALSE))) {
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
  }
  else {
    ADF = NULL
  }
  if (!is.null(tax_table(physeq, FALSE))) {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                        data.frame(as(tax_table(physeq), "matrix")), 
                                        row.names = taxa_names(physeq)))
  }
  else {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                        row.names = taxa_names(physeq)))
  }
  if (requireNamespace("metagenomeSeq")) {
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, 
                                           phenoData = ADF, featureData = TDF, ...)
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
      p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
      p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}


########################
## Big data functions ##
########################

minhead <- function(x){x[1:5,1:5]}

#####################
## Curve functions ##
#####################
# Identify positions of local maxima
# this function is very sensitive some quality control has to be done afterwards
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(-Inf, x)) > 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples the every second location (i.e. switch to TRUE)
  
  y
  # return locations of local maxima
}

localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(Inf, x)) < 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples the every second location (i.e. switch to TRUE)
  
  y
  # return locations of local maxima
}


############
# Plotting #
############
# Custom theme with smaller legends
theme_cust <- function(base_theme = "bw", base_size = 11, half_line = base_size/2){
  if(base_theme == "bw"){
    t <- theme_bw(base_size = base_size)
  }
  
  if(base_theme == "pubr"){
    t <- theme_pubr(base_size = base_size)
  }
  
  t %+replace% theme(legend.position = "right", 
          legend.spacing = unit(half_line / 2, "pt"), # spacing between legends
          legend.key.size = unit(0.7, "lines"), # size of legend symbol box
          legend.box.spacing = unit(1.5 * half_line, "pt"), # spacing between legend and plot
          legend.text = element_text(size = unit(base_size - 4, "pt")), 
          legend.title = element_text(size = unit(base_size - 3, "pt"), hjust = 0)
          )
  
}

# originals from theme_grey
#legend.key.size = unit(1.5, "lines") # size of legend symbol box
#legend.spacing = unit(2 * half_line, "pt") # spacing between legends
#legend.box.spacing = unit(2 * half_line, "pt") # spacing between legend and plot
#plot.title = element_text(size = rel(1.2)


####################################
## Plotting multivariate analysis ##
####################################
# From phyloseq, they have it from NeatMap-package -> cite both if results are reported

## Make into function to save scripting space
## (as we're repeating the same steps for different datasets)

plot_pcoa <- function(pcoa, physeq, plot.axes = c(1,2), axes = plot.axes,
                      colours = NULL, output = F){
  # extract scores and variance explained
  pdataframe <- data.frame(Sample = as.character(row.names(pcoa$vectors)),
                          pcoa$vectors[,axes[1:length(axes)]],
                          stringsAsFactors = F)  # extract site scores
  
  pb.var <- data.frame(Axes = axes,
                       var = round(100 * pcoa$values$Eigenvalues[axes] / sum(pcoa$values$Eigenvalues), 2),
                       stringsAsFactors = F)
                 
  # merge with a selection of meta data
  meta <- data.frame(Sample = as.character(row.names(sample_df(physeq))),
                     sample_df(physeq) %>% dplyr::select(sample.type.year, Season, Year, 
                                                     DnaType, distance.from.mouth, DR.names), 
                     stringsAsFactors = F)
  pdataframe <- merge(pdataframe, meta, by = "Sample")
  pdataframe$Sample <- as.character(pdataframe$Sample)
  
  # extract only factors represented in df for colvec
  colvec <- colvec[names(colvec) %in% as.character(levels(pdataframe$sample.type.year))]
  
  # get legend of plot separately
  
    # main plot
    p <- ggplot(pdataframe, aes_string(x = paste0("Axis.", plot.axes[1]), 
                                                y = paste0("Axis.", plot.axes[2]))) +
        theme_cust() +
        geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
        geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
        geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
                   size = 2.5) +
        coord_fixed(1) + # ensure aspect ratio
        scale_shape_manual(values = c(21,23,25)) +
        scale_alpha_manual(values = c(1,0.7), name = "Nucleic Acid \nType") +
        labs(x = paste0("PC",  plot.axes[1]," [ ", pb.var[pb.var$Axes == plot.axes[1],"var"]," %]"), 
             y = paste0("PC",  plot.axes[2]," [ ", pb.var[pb.var$Axes == plot.axes[2],"var"]," %]")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
               alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
               fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
    
    if(!is.null(colours)){
      p <- p + scale_fill_manual(values = colvec, name = "Habitat Type")
    }
  
  if(output == T){
    return(list(df = pdataframe, var = pb.var, plot = p))
  } else {
    p
  }
}

# Function to plot Bray-Curtis and Jaccard PCoA plots
plot_bray_n_jacc <- function(bray, jacc, colours = colvec, plot.name = NULL, output = F){
  if (is.null(plot.name) == T) {
    stop("'plot.name' needs to be specified.")
  }
  
  # extract scores and variance explained
  pb.scores <- rbind(data.frame(Sample = as.character(row.names(bray$vectors)),
                                Metric = "Bray", bray$vectors[,1:3], stringsAsFactors = F),
                     data.frame(Sample = as.character(row.names(jacc$vectors)),
                                Metric = "Jaccard", jacc$vectors[,1:3], stringsAsFactors = F))  # get first three axes

  
  pb.var <- rbind(data.frame(x = round(100 * bray$values$Eigenvalues[1] / sum(bray$values$Eigenvalues), 2),
                             y = round(100 * bray$values$Eigenvalues[2] / sum(bray$values$Eigenvalues), 2),
                             z = round(100 * bray$values$Eigenvalues[3] / sum(bray$values$Eigenvalues), 2),
                             Metric = "Bray", stringsAsFactors = F),
                  data.frame(x = round(100 * jacc$values$Eigenvalues[1] / sum(jacc$values$Eigenvalues), 2),
                             y = round(100 * jacc$values$Eigenvalues[2] / sum(jacc$values$Eigenvalues), 2),
                             z = round(100 * jacc$values$Eigenvalues[3] / sum(jacc$values$Eigenvalues), 2),
                             Metric = "Jaccard", stringsAsFactors = F))
  pdataframe <- merge(pb.scores, pb.var, by = "Metric")
  
  # calculate centroids and confidence interval of factorial clouds (season-sample type)
  
  # merge with a selection of meta data
  meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                     sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, DnaType), 
                     stringsAsFactors = F)
  pdataframe <- merge(pdataframe, meta, by = "Sample")
  pdataframe$Sample <- as.character(pdataframe$Sample)
  
  # overwrite factor levels
  pdataframe$sample.type.year <- factor(pdataframe$sample.type.year, levels = c("Soil","Sediment",
                                                                                "Soilwater","Hyporheicwater", 
                                                                                "Wellwater","Stream", "Tributary",
                                                                                "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                                "Upriver", "Downriver","RO3", "RO2", "RO1","Deep",
                                                                                "Marine"),
                                        labels = c("Soil","Sediment",
                                                   "Soilwater","Soilwater", 
                                                   "Groundwater","Stream", "Tributary",
                                                   "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
                                                   "Upriver","Downriver",
                                                   "Reservoirs","Reservoirs", "Reservoirs","Reservoirs",
                                                   "Estuary"))
  pdataframe$Season <- factor(pdataframe$Season, levels = c("spring", "summer", "autumn"), 
                              labels = c("Spring", "Summer","Autumn"))
  pdataframe$DnaType <- factor(pdataframe$DnaType, levels = c("DNA", "cDNA"), labels = c("DNA", "RNA"))
  
  # extract only colours that are in data frame
  colvec <- colvec[names(colvec) %in% as.character(levels(pdataframe$sample.type.year))]
  
  # get legend of plot separately
  (
    bot.leg <-
      get_legend(
        ggplot(pdataframe %>% filter(Metric == "Bray"), aes(x = Axis.1, y = Axis.2)) +
          theme_bw() +
          geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
                     size = 2.5) +
          scale_fill_manual(values = colvec, name = "Sample Type") +
          theme(legend.position = "bottom") +
          scale_shape_manual(values = c(21, 23, 25)) +
          scale_alpha_manual(values = c(1, 0.5), name = "Nucleic Acid Type") +
          guides(shape = guide_legend(order = 1),
                 alpha = guide_legend(order = 2), fill = "none")
      )
  )
  
  # main plot
  (pcoa.bray <- ggplot(pdataframe %>% filter(Metric == "Bray"), aes(x = Axis.1, y = Axis.2)) +
      theme_bw() +
      geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
      geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
      geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
                 size = 2.5) +
      scale_fill_manual(values = colvec, name = "Sample Type") +
      scale_shape_manual(values = c(21,23,25)) +
      scale_alpha_manual(values = c(1,0.7), name = "Nucleic Acid \nType") +
      coord_fixed(1) + # ensure aspect ratio
      labs(x = paste("PC1 [", pdataframe[pdataframe$Metric == "Bray",]$x,"%]"), 
           y = paste("PC2 [", pdataframe[pdataframe$Metric == "Bray",]$y,"%]"),
           title = "Bray-Curtis (sqrt) (is euclidean = T)") +
      theme(legend.key.size = unit(1.5, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(fill = guide_legend(override.aes=list(shape=21)),
             alpha = "none", shape = "none"))
  
  (pcoa.jac <- ggplot(pdataframe %>% filter(Metric == "Jaccard"), aes(x = Axis.1, y = Axis.2)) +
      theme_bw() +
      geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
      geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
      geom_point(aes(fill = sample.type.year, shape = Season, alpha = DnaType),
                 size = 2.5) +
      scale_fill_manual(values = colvec, name = "Sample Type") +
      scale_shape_manual(values = c(21,23,25)) +
      scale_alpha_manual(values = c(1,0.7), name = "Nucleic Acid \nType") +
      coord_fixed(1) + # ensure aspect ratio
      labs(x = paste("PC1 [", pdataframe[pdataframe$Metric == "Jaccard",]$x,"%]"), 
           y = paste("PC2 [", pdataframe[pdataframe$Metric == "Jaccard",]$y,"%]"),
           title = "Jaccard (is euclidean = T)") +
      theme(legend.key.size = unit(1.5, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(fill = guide_legend(override.aes=list(shape=21)),
             alpha = "none", shape = "none"))
  
  # for centroids
  #geom_errorbar(aes(ymax = upper.ci.pc2, ymin = lower.ci.pc2, colour = sample.type.year),
  #              width = 0.02, size = 0.3) +
  #geom_errorbarh(aes(xmax = upper.ci.pc1, xmin = lower.ci.pc1, colour = sample.type.year),
  #               height = 0.02, size = 0.3) +
  #scale_colour_manual(values = colvec, name = "Sample Type") +
  
  # arrange main plot and legend
  (pcoa.arr <- ggarrange(ggarrange(pcoa.bray, pcoa.jac,
                                   ncol = 2,
                                   common.legend = T,
                                   labels = c("a", "b"),
                                   legend = "right"
  ), bot.leg, nrow = 2, heights = c(3, 0.2)))
  
  ggsave(paste0("./Figures/General/PCoA_", plot.name,"_SampleType.tiff"), pcoa.arr,
         width = 30, height = 15, unit = "cm")
  ggsave(paste0("./Figures/General/PCoA_", plot.name,"_SampleType.png"), pcoa.arr,
         width = 30, height = 15, unit = "cm")
  
  if(output == T){
    return(list(df = pdataframe, bray = pcoa.bray, jacc = pcoa.jac))
  }
  
}

# define function for confidence interval
conf.int <- function(x){
  mean <- mean(x, na.rm = T)
  SD <- sd(x, na.rm = T)
  N <- length(x)
  
  SE <- SD / sqrt(N)
  ci <- qt((0.975/2 + 0.5), N - 1) * SE
  return(ci)
}

scatter_panels <- function(data, labs = c(x,y)){
  (
    main <-
      ggplot(data[data$panels == "main",], aes(
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = labs[1], 
           y = labs[2]) +
      lims(y = c(min(data$mean - data$stdev, na.rm = T),
                 max(data$mean + data$stdev, na.rm = T))) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 10)
      )
  )
  
  # side panel
  (
    side <-
      ggplot(data[data$panels == "side",], aes(
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = labs[1], 
           y = labs[2]) +
      lims(y = c(min(data$mean - data$stdev, na.rm = T),
                 max(data$mean + data$stdev, na.rm = T))) +
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
  
  # only dist
  
  p <- ggarrange(main, 
                 side,
                 widths = c(3,1),
                 ncol = 2, nrow = 1, 
                 common.legend = T,
                 align = "h",
                 font.label = list(size = 10))
  
  out <- list(plot = p, main = main, side = side)
  return(out)
}

## Function to calculate and plot DNA-RNA distance within PCoA space

# two metrics version
dist.dnarna.metrics <- function(data, save.name = NULL, output = F, dimensions = 2){
  if (is.null(save.name) == T) {
    stop("'save.name' needs to be specified.")
  }
  if(dimensions <= 1 | dimensions > 3){
    stop("'dimensions' needs to be 2 or 3 to calculate distance in ordination space.")
  }
  
  # correct a few wrong sample names for matching DNA and RNA
  data[data$Sample == "RO2R52R", "Sample"] <- "RO2.52R"
  data[data$Sample == "SWR34R", "Sample"] <- "SW34R"
  data[data$Sample == "RO2.36pD", "Sample"] <- "RO2.36D"
  data[data$Sample == "RO2.36pR", "Sample"] <- "RO2.36R"
  data[data$Sample == "RO2111.60mD", "Sample"] <- "RO2111.90mD"
  data[data$Sample == "RO2.30DPR", "Sample"] <- "RO2.30R" # two DNA
  data[data$Sample == "RO301.HypoR", "Sample"] <- "RO31.HypoR"
  data[data$Sample == "RO301R", "Sample"] <- "RO31R" 
  data[data$Sample == "RO304R", "Sample"] <- "RO34R" 
  data[data$Sample == "RO307R", "Sample"] <- "RO37R" 
  data[data$Sample == "L230R", "Sample"] <- "L330R" # L230 does not exist
  
  # remove Ds and Rs to match counterpart samples
  data$ID[data$DnaType == "DNA"] <- str_replace(data$Sample[data$DnaType == "DNA"], "D$", "")
  data$ID[data$DnaType == "RNA"] <- str_replace(data$Sample[data$DnaType == "RNA"], "R$", "")
  
  # export table to look at point positions in GIS
  #write.table(pb.scores, "./Output/BrayCurtis_scores_withmeta.csv", sep = ",", dec = ".", row.names = F)
  
  if(dimensions == 2){
    # calculate mean coordinates for duplicates
    sum <- data %>% 
      filter(!Year == 2015) %>% 
      dplyr::group_by(ID, DnaType, Metric) %>%
      dplyr::summarise(x = mean(Axis.1), y = mean(Axis.2),
                       n = n()) %>%
      ungroup()
    
    setDT(sum)
    # Calculate distance
    temp <- dcast(sum, ID + Metric ~ DnaType, value.var = c("x","y"))
    temp[, distance := sqrt((x_DNA - x_RNA)^2 + (y_DNA - y_RNA)^2)]
    temp <- temp[!is.na(temp$distance),]
    
  } else if(dimensions == 3){
    # calculate mean coordinates for duplicates
    sum <- data %>% 
      filter(!Year == 2015) %>% 
      dplyr::group_by(ID, DnaType, Metric) %>%
      dplyr::summarise(x = mean(Axis.1), y = mean(Axis.2), z = mean(Axis.3),
                       n = n()) %>%
      ungroup()
    
    setDT(sum)
    # Calculate distance
    temp <- dcast(sum, ID + Metric ~ DnaType, value.var = c("x","y","z"))
    temp[, distance := sqrt((x_DNA - x_RNA)^2 + (y_DNA - y_RNA)^2 + (z_DNA - z_RNA)^2)]
    temp <- temp[!is.na(temp$distance),]
  }
  
  
  # combine back with categories
  dist.dr <- temp[data, c("sample.type.year",
                                "Year", "Season" ) := list(i.sample.type.year,
                                                           i.Year, i.Season), on = .(ID)]
  # add new column to split plot into main and side panel
  dist.dr[, panels := "main"]
  dist.dr[sample.type.year == "Tributary" |
            sample.type.year == "Lake" |
            sample.type.year == "Riverine \nLakes" |
            sample.type.year == "Sediment", panels := "side"]
  
  write.table(dist.dr, paste0("./Output/", save.name ,"_PCoA_distance.csv"),
              sep = ";", dec = ".", row.names = F)
  
  # calculate confidence interval and means of sample type and season combinations
  dist.dr <- dist.dr[, .(mean =  mean(distance, na.rm = T),
                         conf.int = conf.int(distance)), by = .(Metric, sample.type.year, Season, panels)]
  
  
  # Plot
  # Bray Curtis
  # plot main plot
  (
    main.b <-
      ggplot(dist.dr[panels == "main" & Metric == "Bray", ], aes(
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - conf.int, ymax = mean + conf.int, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = "Sample type", 
           y = paste0("Distance between DNA and RNA \nin ordination space"),
           title = "Bray-Curtis (sqrt)") +
      lims(y = c(min(dist.dr[Metric == "Bray",]$mean - dist.dr[Metric == "Bray",]$conf.int, na.rm = T),
                 max(dist.dr[Metric == "Bray",]$mean + dist.dr[Metric == "Bray",]$conf.int, na.rm = T))) +
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
      ggplot(dist.dr[panels == "side" & Metric == "Bray", ], aes(
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - conf.int, ymax = mean + conf.int, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = "Sample type", 
           y = paste0("Distance between DNA and RNA \nin ordination space")) +
      lims(y = c(min(dist.dr[Metric == "Bray",]$mean - dist.dr[Metric == "Bray",]$conf.int, na.rm = T),
                 max(dist.dr[Metric == "Bray",]$mean + dist.dr[Metric == "Bray",]$conf.int, na.rm = T))) +
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
  
  # Jaccard
  # plot main plot
  (
    main.j <-
      ggplot(dist.dr[panels == "main" & Metric == "Jaccard", ], aes(
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - conf.int, ymax = mean + conf.int, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = "Sample type", 
           y = paste0("Distance between DNA and RNA \nin ordination space"),
           title = "Jaccard") +
      lims(y = c(min(dist.dr[Metric == "Jaccard",]$mean - dist.dr[Metric == "Jaccard",]$conf.int, na.rm = T),
                 max(dist.dr[Metric == "Jaccard",]$mean + dist.dr[Metric == "Jaccard",]$conf.int, na.rm = T))) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 10)
      )
  )
  
  # side panel
  (
    side.j <-
      ggplot(dist.dr[panels == "side" & Metric == "Jaccard", ], aes(
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - conf.int, ymax = mean + conf.int, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = "Sample type", 
           y = paste0("Distance between DNA and RNA \nin ordination space")) +
      lims(y = c(min(dist.dr[Metric == "Jaccard",]$mean - dist.dr[Metric == "Jaccard",]$conf.int, na.rm = T),
                 max(dist.dr[Metric == "Jaccard",]$mean + dist.dr[Metric == "Jaccard",]$conf.int, na.rm = T))) +
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
  
  
  # combine both plots
  (p <- ggarrange(
    ggarrange(main.b, 
              side.b,
              widths = c(3,1),
              ncol = 2, nrow = 1, 
              common.legend = T,
              legend = "top",
              align = "h",
              font.label = list(size = 10)),
    ggarrange(main.j, 
              side.j,
              widths = c(3,1),
              ncol = 2, nrow = 1, 
              common.legend = F,
              legend = "none",
              align = "h",
              font.label = list(size = 10)),
    nrow = 2))
  # add x axis title to be in the middle of two panels
  (p <- annotate_figure(p, bottom = text_grob("Sample Type")))
  
  # save
  ggsave(paste0("./Figures/General/DNARNA_withinPCoA_distance_pnts_", save.name,".tiff"), 
         p, width = 18, height = 13, unit = "cm")
  ggsave(paste0("./Figures/General/DNARNA_withinPCoA_distance_pnts_", save.name,".png"), 
         p, width = 18, height = 13, unit = "cm")
  
  if(output == T){
    return(list(df = dist.dr, plot = ggarrange(main.b, 
                                               side.b,
                                               widths = c(3,1),
                                               ncol = 2, nrow = 1, 
                                               common.legend = T,
                                               legend = "top",
                                               align = "h",
                                               font.label = list(size = 10))))  
  }
}


## Function to calculate and plot DNA-RNA distance within PCoA space
dist.dnarna <- function(bray, save.name = NULL, output = F){
  if (is.null(save.name) == T) {
    stop("'save.name' needs to be specified.")
  }
  
  pb.scores <- data.frame(Sample = as.character(row.names(bray$vectors)),
                          bray$vectors, stringsAsFactors = F)  # get first three axes
  # merge with a selection of meta data
  meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                     sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, 
                                                     DnaType, distance.from.mouth, DR.names), 
                     stringsAsFactors = F)
  data <- merge(pb.scores, meta, by = "Sample")
  data$Sample <- as.character(data$Sample)
  
  # export table to look at point positions in GIS
  #write.table(pb.scores, "./Output/BrayCurtis_scores_withmeta.csv", sep = ",", dec = ".", row.names = F)
  
    setDT(data)
    # melt datatable
    temp <- melt(data, id.vars = c("DR.names","DnaType"), measure.vars = patterns("^Axis."),
                 variable.name = "Axis", value.name = "Coordinates")
    temp <- dcast(temp, DR.names + Axis ~ DnaType, value.var = c("Coordinates"))
    # remove NAs
    temp <- na.omit(temp)
    
    # Calculate distance of all axes
    temp[, pnt.dist := (abs(DNA - RNA))^2] # calculate point distance for each axis and square root
    temp <- temp[, .(sum.dist = sum(pnt.dist)), by = .(DR.names)] # sum all axes
    temp <- temp[, dist := sqrt(sum.dist)]
    
    # Calculate distance
    #temp[, distance.1D := abs(DNA - RNA)]
    
    #temp.2d <- dcast(temp, DR.names ~ Axis, value.var = c("DNA","RNA"))
    #temp.2d[, distance.12 := sqrt(abs((DNA_Axis.1 - RNA_Axis.1))^2 + abs((DNA_Axis.2 - RNA_Axis.2))^2)]
    #temp.2d[, distance.13 := sqrt(abs((DNA_Axis.1 - RNA_Axis.1))^2 + abs((DNA_Axis.3 - RNA_Axis.3))^2)]
    
  # combine back with categories
  dist.dr <- temp[data, c("sample.type.year",
                          "Year", "Season") := list(i.sample.type.year,
                                                     i.Year, i.Season), on = .(DR.names)]
  #dist.2d <- temp.2d[data, c("sample.type.year",
  #                           "Year", "Season") := list(i.sample.type.year,
  #                                                     i.Year, i.Season), on = .(DR.names)]
  indiv.df <- dist.dr
  
  # calculate confidence interval and means of sample type and season combinations
  dist.dr <- dist.dr[, .(mean =  mean(dist, na.rm = T),
                         conf.int = conf.int(dist),
                         stdev = sd(dist, na.rm = T)),
                     by = .(sample.type.year, Season)]
  
  # add new column to split plot into main and side panel
  dist.dr[, panels := "main"]
  dist.dr[sample.type.year == "Tributary" |
            sample.type.year == "Lake" |
            sample.type.year == "Riverine \nLakes" |
            sample.type.year == "Sediment", panels := "side"]
  
  
  
  # Plot
  # Bray Curtis
  # plot main plot
  (
    main.b <-
      ggplot(dist.dr[panels == "main", ], aes( #& Axis == "Axis.2"
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = "Habitat type", 
           y = paste0("Distance in \nordination space (PC2)")) +
      #lims(y = c(min(dist.dr[Axis == "Axis.2",]$mean - dist.dr[Axis == "Axis.2",]$stdev, na.rm = T),
      #           max(dist.dr[Axis == "Axis.2",]$mean + dist.dr[Axis == "Axis.2",]$stdev, na.rm = T))) +
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
      ggplot(dist.dr[panels == "side" , ], aes( # & Axis == "Axis.2"
        x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = "Habitat type", 
           y = paste0("Distance in ordination space (PC2)")) +
      #lims(y = c(min(dist.dr[Axis == "Axis.2",]$mean - dist.dr[Axis == "Axis.2",]$stdev, na.rm = T),
      #           max(dist.dr[Axis == "Axis.2",]$mean + dist.dr[Axis == "Axis.2",]$stdev, na.rm = T))) +
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
  
  #(
  #  main.s <-
  #    ggplot(dist.dr[panels == "main", ], aes(
  #      x = sample.type.year, y = mean, fill = Season
  #    )) +
  #    theme_cust(base_theme = "pubr") +
  #    geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
  #                  position = position_dodge(0.7), width = 0) +
  #    geom_jitter(aes(fill = Season), shape = 21, 
  #                position = position_dodge(0.7), size = 2.5) +
  #    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
  #    scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
  #    labs(x = "Habitat type", 
  #         y = paste0("Distance in \nordination space (PC3)")) +
  #    lims(y = c(min(dist.dr[Axis == "Axis.3",]$mean - dist.dr[Axis == "Axis.3",]$stdev, na.rm = T),
  #               max(dist.dr[Axis == "Axis.3",]$mean + dist.dr[Axis == "Axis.3",]$stdev, na.rm = T))) +
  #    theme(
  #      axis.text.x = element_text(angle = 45, hjust = 1),
  #      axis.text = element_text(size = 8),
  #      axis.title.x = element_blank(),
  #      axis.title = element_text(size = 10)
  #    )
  #)
  
  # side panel
  #(
  #  side.s <-
  #    ggplot(dist.dr[panels == "side" & Axis == "Axis.3", ], aes(
  #      x = sample.type.year, y = mean, fill = Season
  #    )) +
  #    theme_cust(base_theme = "pubr") +
  #    geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
  #                  position = position_dodge(0.7), width = 0) +
  #    geom_jitter(aes(fill = Season), shape = 21, 
  #                position = position_dodge(0.7), size = 2.5) +
  #    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
  #    scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
  #    labs(x = "Habitat type", 
  #         y = paste0("Distance in ordination space (PC3)")) +
  #    lims(y = c(min(dist.dr[Axis == "Axis.3",]$mean - dist.dr[Axis == "Axis.3",]$stdev, na.rm = T),
  #               max(dist.dr[Axis == "Axis.3",]$mean + dist.dr[Axis == "Axis.3",]$stdev, na.rm = T))) +
  #    theme(
  #      axis.title.x = element_blank(),
  #      axis.text.x = element_text(angle = 45, hjust = 1),
  #      axis.text = element_text(size = 8),
  #      axis.title.y = element_blank(),
  #      axis.text.y = element_blank(),
  #      axis.line.y = element_blank(),
  #      axis.ticks.y = element_blank()
  #    )
  #)
  
  # combine both plots
  (p <- ggarrange(
    ggarrange(main.b, 
              side.b,
              widths = c(3,1),
              ncol = 2, nrow = 1, 
              common.legend = T,
              legend = "top",
              align = "h",
              font.label = list(size = 10)),
    #ggarrange(main.s, 
    #          side.s,
    #          widths = c(3,1),
    #          ncol = 2, nrow = 1, 
    #          common.legend = T,
    #          legend = "none",
    #          align = "h",
    #          font.label = list(size = 10)), 
    nrow = 2, align = "hv", common.legend = T)
  )
    
  # add x axis title to be in the middle of two panels
  (p <- annotate_figure(p, bottom = text_grob("Habitat Type")))
  
  if(output == T){
    return(list(df = dist.dr, indiv.df = indiv.df, raw.df = data, #dist.2d = dist.2d,
                plot.main = main.b, plot.side = side.b))  
  }
}

## Function to calculate and plot DNA-RNA dissimilarity (not PCoA space)

melt.dist <- function(x, var.name = "dist"){
  # convert distance matrix into long format
  dist.mat <- as(x, "matrix")
  xy <- t(combn(colnames(dist.mat), 2))
  dist.df <- data.frame(xy, dist=dist.mat[xy], stringsAsFactors = F)
  colnames(dist.df)[1:3] <- c("Sample.x", "Sample.y", var.name)
  return(dist.df)
}

subset_asv <- function (physeq, names){
    OTU <- access(physeq, "otu_table", TRUE)
    if (!taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    OTU <- as(OTU, "matrix")
    new <- subset(OTU, rownames(OTU) %in% names)
    
    if (inherits(physeq, "otu_table")) {
      return(otu_table(new))
    } else {
      otu_table(physeq) <- otu_table(new, taxa_are_rows = T)
      return(physeq)
    }
}


dissim.dnarna <- function(physeq, save.name = NULL, output = F){
  if (is.null(save.name) == T) {
    stop("'save.name' needs to be specified.")
  }
  pb.mat <- otu_mat(physeq)
  #pb.mat <- log2(pb.mat + 1)
  bray <- vegdist(pb.mat, method = "bray")
  bray <- sqrt(bray) # make euclidean
  
  # convert distance matrix into long format
  dist.df <- melt.dist(bray)
  
  # Get meta data to rename DNA and RNA data
  meta <- data.frame(Sample = as.character(row.names(sample_df(physeq))),
                     sample_df(physeq) %>% dplyr::select(DR.names, DnaType, Year, Season, sample.type.year), 
                     stringsAsFactors = F)
  
  # add meta data for both sample x and sample y
  dist.df <- merge(dist.df, meta, by.x =  "Sample.x", by.y = "Sample")
  dist.df <- merge(dist.df, meta %>% select(DnaType, Sample, Year, DR.names), by.x =  "Sample.y", by.y = "Sample")
  
  # omit all samples of 2015 (no RNA was taken, and sample name strategy changed -> creates duplicates)
  dist.df <- dist.df[dist.df$Year.x != 2015,]
  dist.df <- dist.df[dist.df$Year.y != 2015,]
  
  # keep all rows where Sample.x and Sample.y are the same
  dist.df <- dist.df[dist.df$DR.names.x == dist.df$DR.names.y,]
  dist.df <- dist.df[!(dist.df$DnaType.x == dist.df$DnaType.y),] # omit all distances between same DnaType
  
  # check duplicates
  #dist.df[dist.df$ID.x %in% dist.df[duplicated(dist.df$ID.x),]$ID.x,]
  
  dist.dr <- dist.df %>% select(ID = DR.names.x, Year = Year.x, Season, sample.type.year, dist)
  setDT(dist.dr)
  
    # add new column to split plot into main and side panel
  dist.dr[, panels := "main"]
  dist.dr[sample.type.year == "Tributary" |
            sample.type.year == "Lake" |
            sample.type.year == "HeadwaterLakes" |
            sample.type.year == "Sediment", panels := "side"]
  
  # calculate confidence interval and means of sample type and season combinations
  plot.df <- dist.dr[, .(mean =  mean(dist, na.rm = T),
                         conf.int = conf.int(dist),
                         stdev = sd(dist, na.rm = T)), by = .(sample.type.year, Season, panels)]
  
  # Plot Bray
  (
    main.b <-
      ggplot(plot.df[panels == "main", ], 
             aes(x = sample.type.year, y = mean, fill = Season
      )) +
      theme_cust(base_theme = "pubr")
    +
      geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
                    position = position_dodge(0.7), width = 0) +
      geom_jitter(aes(fill = Season), shape = 21, 
                  position = position_dodge(0.7), size = 2.5) +
      scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
      scale_colour_manual(values = c("#009E73", "#FFAA1D", "#D55E00")) +
      labs(x = "Sample type", y = "Pair-wise \nBray-Curtis dissimilarity") +
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
      ggplot(plot.df[panels == "side", ], 
             aes(x = sample.type.year, y = mean, fill = Season
             )) +
      theme_cust(base_theme = "pubr") +
      geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev, colour = Season),
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
  
  # combine both plots
  (p <- ggarrange(main.b, 
              side.b,
              widths = c(3,1),
              ncol = 2, nrow = 1, 
              common.legend = T,
              legend = "top",
              align = "h",
              font.label = list(size = 10)))
  # add x axis title to be in the middle of two panels
  (p <- annotate_figure(p, bottom = text_grob("Sample Type")))
  
  # save
  ggsave(paste0("./Figures/General/DNARNA_dissimilarity_", save.name,".tiff"), 
         p, width = 18, height = 13, unit = "cm")
  ggsave(paste0("./Figures/General/DNARNA_dissimilarity_", save.name,".png"), 
         p, width = 18, height = 13, unit = "cm")
  
  
  #wide.format <- dcast(dist.dr, ID + sample.type.year + Season ~ , value.var = c("dist"))
  
  if(output == T){
    out <- list(original.df = dist.dr, #wide = wide.format,
                sum.df = plot.df, plot.main = main.b, plot.side = side.b)
    
    return(out)  
  }
}



################################################################################
# Adapted from NeatMap-package and phyloseq.
# Vectorized for speed and simplicity, also only calculates theta and not r.
RadialTheta <- function(pos){
  pos = as(pos, "matrix")
  xc  = mean(pos[, 1])
  yc  = mean(pos[, 2])
  theta = atan2((pos[, 2] - yc), (pos[, 1] - xc))
  names(theta) <- rownames(pos)
  return(theta)
}

# Translate gneiss functions into R
#mean_nich_estimator <- function(abundances, gradient){
  #Estimates the mean niche along a gradient of an organism.
  
  #Calculates the mean niche of an organism along a gradient.
  #This is done by calculating the mean gradient values that
  #an organism is observed in.
  
  #Specifically, this module calculates the following
  #.. math::
  #  f(g , x) =
  #  \sum\limits_{i=1}^N g_i \frac{x_i}{\sum\limits_{j=1}^N x_j}
  
  
  #Where :math:`N` is the number of samples, :math:`x_i` is the proportion of
  #species :math:`x` in sample :math:`i`, :math:`g_i` is the gradient value
  #at sample `i`.
  
  #Parameters
  #abundances : pd.DataFrame or pd.Series, np.float
  #Vector of fraction abundances of an organism over a list of samples.
  #gradient : pd.Series, np.float
  #Vector of numerical gradient values.
 # len_abundances <- length(abundances)
#  len_gradient <- length(gradient)
#  if(len_abundances != len_gradient){
#    stop("Length of 'abundances' does not match the length of 'gradient'")
#  }
#  if(any(is.na(gradient)) == T){
#    stop("'gradient' cannot have NAs")
#  }
# v <- abundances / sum(abundances)
#  m <- 
#  
#}







# Original boxplots for PCoA distance
#(
#  main.b <-
#    ggplot(dist.dr[panels == "main" & Metric == "Bray", ], aes(
#      x = sample.type.year, y = distance, fill = Season
#    )) +
#    theme_pubr() +
#    geom_boxplot(width = 0.7, outlier.size = 0.5, size = 0.3) + # outlier.alpha = 0
#    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
#    new_scale_fill() +
    #geom_point(aes(fill = as.character(Year)), position = position_jitterdodge(), colour = "black", shape = 21) + #alpha = 0.5,
    #scale_fill_manual(name = "Year", values = c("white","gray20")) +
#    labs(x = "Sample type", 
#         y = paste0("Distance between DNA and RNA \nin ordination space (", dimensions, "D)")) +
    #facet_grid(.~Year, scales = "free") +
#    theme(
#      axis.text.x = element_text(angle = 45, hjust = 1),
#      axis.text = element_text(size = 8),
#      axis.title.x = element_blank(),
#      axis.title = element_text(size = 10)
#    )
#)

# side panel
#(
#  side.b <-
#    ggplot(dist.dr[panels == "side" & Metric == "Bray", ], aes(
#      x = sample.type.year, y = distance, fill = Season
#    )) +
#    theme_pubr() +
#    geom_boxplot(width = 0.7, outlier.size = 0.5, size = .3) + # outlier.alpha = 0
#    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
#   new_scale_fill() +
    #geom_point(
    #  aes(fill = as.character(Year)),
    #  position = position_jitterdodge(),
    #  colour = "black",
    #  shape = 21
    #) + #alpha = 0.5,
    #scale_fill_manual(name = "Year", values = c("white", "gray20")) +
#    labs(x = "Sample type", 
#         y = paste0("Distance between DNA and RNA \nin ordination space (", dimensions, "D)")) +
#    lims(y = c(0, 0.6)) +
    #facet_grid(.~Year, scales = "free") +
#    theme(
#      axis.title.x = element_blank(),
#      axis.text.x = element_text(angle = 45, hjust = 1),
#      axis.text = element_text(size = 8),
#      axis.title.y = element_blank(),
#      axis.text.y = element_blank(),
#      axis.line.y = element_blank(),
#      axis.ticks.y = element_blank()
#    )
#)

# Jaccard
# plot main plot
#(
#  main.j <-
#    ggplot(dist.dr[panels == "main" & Metric == "Jaccard", ], aes(
##      x = sample.type.year, y = distance, fill = Season
#    )) +
#    theme_pubr() +
#    geom_boxplot(width = 0.7, outlier.size = 0.5, size = 0.3) + # outlier.alpha = 0
#    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
##    new_scale_fill() +
#    #geom_point(aes(fill = as.character(Year)), position = position_jitterdodge(), colour = "black", shape = 21) + #alpha = 0.5,
#    #scale_fill_manual(name = "Year", values = c("white","gray20")) +
#    labs(x = "Sample type", 
#         y = paste0("Distance between DNA and RNA \nin ordination space (", dimensions, "D)")) +
    #facet_grid(.~Year, scales = "free") +
 #   theme(
#      axis.text.x = element_text(angle = 45, hjust = 1),
#      axis.text = element_text(size = 8),
#      axis.title.x = element_blank(),
#      axis.title = element_text(size = 10)
#    )
#)

# side panel
#(
#  side.j <-
#    ggplot(dist.dr[panels == "side" & Metric == "Jaccard", ], aes(
#      x = sample.type.year, y = distance, fill = Season
#    )) +
#    theme_pubr() +
#    geom_boxplot(width = 0.7, outlier.size = 0.5, size = .3) + # outlier.alpha = 0
#    scale_fill_manual(values = c("#009E73", "#F0E442", "#D55E00")) + # colour-blind friendly
#    new_scale_fill() +
    #geom_point(
    #  aes(fill = as.character(Year)),
    #  position = position_jitterdodge(),
    #  colour = "black",
    #  shape = 21
    #) + #alpha = 0.5,
#    #scale_fill_manual(name = "Year", values = c("white", "gray20")) +
#   labs(x = "Sample type", 
#         y = paste0("Distance between DNA and RNA \nin ordination space (", dimensions, "D)")) +
#    lims(y = c(0, 0.6)) +
#    #facet_grid(.~Year, scales = "free") +
#    theme(
#      axis.title.x = element_blank(),
#      axis.text.x = element_text(angle = 45, hjust = 1),
#      axis.text = element_text(size = 8),
#      axis.title.y = element_blank(),
#      axis.text.y = element_blank(),
#      axis.line.y = element_blank(),
#      axis.ticks.y = element_blank()
#    )
#)