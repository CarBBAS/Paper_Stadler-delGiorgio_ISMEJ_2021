
# ecologically ordered heat map

tax.order <- llply(ls.asvs, function(x){
  # extract abundance group specific OTUs
  df <- rel.df[OTU %in% x,]
  com.mat <- dcast(df, ID ~ OTU, value.var = "css.reads")
  setDF(com.mat)
  rownames(com.mat) <- com.mat$ID; com.mat$ID <- NULL # remove ID col
  # make PCoA to get OTU order in heatmap
  pb.bray <- vegdist(com.mat, method = "bray")
  pb.bray <- sqrt(pb.bray) # make Euclidean
  bray.pcoa <- ape::pcoa(pb.bray)
  # extract scores of first three axes
  sitesDF <- data.frame(bray.pcoa$vectors[,1:3])
  speciesDF <- data.frame(wascores(sitesDF, w = com.mat), stringsAsFactors = F)
  speciesDF <- na.omit(speciesDF)
  taxa.order <- row.names(speciesDF)[order(RadialTheta(speciesDF))]
  return(taxa.order)
}, .parallel = T)

rel.df$ID <- factor(rel.df$ID, levels = c("Soil_DNA", "Soil_cDNA",
                                          "Sediment_DNA", "Sediment_cDNA",
                                          "Soilwater_DNA","Soilwater_cDNA",
                                          "Hyporheicwater_DNA", "Hyporheicwater_cDNA",
                                          "Stream_DNA", "Stream_cDNA",
                                          "Tributary_DNA", "Tributary_cDNA",
                                          "HeadwaterLakes_DNA", "HeadwaterLakes_cDNA",
                                          "Lake_DNA", "Lake_cDNA",
                                          "Upriver_DNA", "Upriver_cDNA",
                                          "RO3_DNA", "RO3_cDNA",
                                          "RO2_DNA", "RO2_cDNA",
                                          "RO1_DNA", "RO1_cDNA",
                                          "Deep_DNA", "Deep_cDNA",
                                          "Downriver_DNA", "Downriver_cDNA",
                                          "Marine_DNA", "Marine_cDNA"),
                    labels = c("Soil DNA", "Soil RNA",
                               "Sediment DNA", "Sediment RNA",
                               "Soilwater DNA","Soilwater RNA",
                               "Hyporheicwater DNA", "Hyporheicwater RNA",
                               "Stream DNA", "Stream RNA",
                               "Tributary DNA", "Tributary RNA",
                               "Headwater Lakes DNA", "Headwater Lakes RNA",
                               "Lake DNA", "Lake RNA",
                               "Upriver DNA", "Upriver RNA",
                               "RO3 DNA", "RO3 RNA",
                               "RO2 DNA", "RO2 RNA",
                               "RO1 DNA", "RO1 RNA",
                               "Hypolimnion DNA", "Hypolimnion RNA",
                               "Downriver DNA", "Downriver RNA",
                               "Estuary DNA", "Estuary RNA"))

x.labs<-c("D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R","D", "R",
          "D", "R")

heat.plots <- llply(tax.order, function(x){
  #if(length(x) > 5000){
  #  sub.x <- sample(x, 5000)
  #  df <- sum[sum$ASV %in% sub.x,]
  #} else {
  #  
  #}
  df <- rel.df[rel.df$OTU %in% x,]
  df$OTU <- factor(df$OTU, levels = x)
  
  # get max value of all
  max.lim <- round(max(log2(rel.df$css.reads+1)),0)
  
  #low <- "#000033"; high <- "#66CCFF"
  p <- ggplot(df, aes(x = ID, y = OTU, fill = log2(css.reads + 1))) +
    theme_bw() +
    geom_raster() +
    labs(x = "Habitat Type", y = "OTUs") +
    scale_fill_viridis_c(option = "magma", name = "Mean\nReads",
                         limits = c(min(log2(df$css.reads + 1)),
                                    max.lim),
                         breaks = c(min(log2(df$css.reads + 1)),
                                    max.lim/ 2,
                                    max.lim),
                         labels = c(as.character(min(log2(df$css.reads + 1))),
                                    as.character(max.lim / 2),
                                    as.character(max.lim))) +
    scale_x_discrete(labels = x.labs) +
    coord_cartesian(ylim = c(1,
                             length(levels(df$OTU))), clip = "off") +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.title.x = element_blank(),
          legend.position = "right", plot.title = element_text(size = 12)) +
    guides(fill = guide_colourbar(barwidth = unit(10, "pt"),
                                  barheight = unit(50, "pt")))
  #scale_fill_gradient(low = low, high = high)
  return(p)
}, .parallel = T)

(univ.rare <- heat.plots$universal.rare + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = "Universally rare"))

(spec <- heat.plots$specialist + 
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank()
    ) + labs(title = "Specialist"))

(rare.shif <- heat.plots$rare.shifter +
    theme(plot.margin=unit(c(2,2,18,2), "mm")) +
    labs(title = "Rare shifter") +
    annotate("segment", x = 0.75, xend = 2.25, y = -30, yen = -30) +
    annotate("text", x = 1.6, y = - 40, label = "Soil", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 2.75, xend = 4.25, y = -30, yen = -30) +
    annotate("text", x = 3.6, y = - 40, label = "Sediment", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 4.75, xend = 6.25, y = -30, yen = -30) +
    annotate("text", x = 5.6, y = - 40, label = "Soilwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 6.75, xend = 8.25, y = -30, yen = -30) +
    annotate("text", x = 7.6, y = - 40, label = "Hyporheicwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 8.75, xend = 10.25, y = -30, yen = -30) +
    annotate("text", x = 9.6, y = - 40, label = "Stream", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 10.75, xend = 12.25, y = -30, yen = -30) +
    annotate("text", x = 11.6, y = - 40, label = "Tributary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 12.75, xend = 14.25, y = -30, yen = -30) +
    annotate("text", x = 13.6, y = - 40, label = "Headwater Lakes", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 14.75, xend = 16.25, y = -30, yen = -30) +
    annotate("text", x = 15.6, y = -40, label = "Lake", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 16.75, xend = 18.25, y = -30, yen = -30) +
    annotate("text", x = 17.6, y = -40, label = "Upriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 18.75, xend = 20.25, y = -30, yen = -30) +
    annotate("text", x = 19.6, y = -40, label = "RO3", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 20.75, xend = 22.25, y = -30, yen = -30) +
    annotate("text", x = 21.6, y = -40, label = "RO2", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 22.75, xend = 24.25, y = -30, yen = -30) +
    annotate("text", x = 23.6, y = -40, label = "RO1", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 24.75, xend = 26.25, y = -30, yen = -30) +
    annotate("text", x = 25.6, y = -40, label = "Hypolimnion", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 26.75, xend = 28.25, y = -30, yen = -30) +
    annotate("text", x = 27.6, y = -40, label = "Downriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 28.75, xend = 30.25, y = -30, yen = -30) +
    annotate("text", x = 29.6, y = -40, label = "Estuary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5))


(ab.shif <- heat.plots$abundant.shifter + 
    theme(axis.title.y = element_blank(),
          plot.margin=unit(c(2,2,18,2), "mm")) + 
    labs(title = "Abundant shifter") +
    annotate("segment", x = 0.75, xend = 2.25, y = -70, yen = -70) +
    annotate("text", x = 1.6, y = - 90, label = "Soil", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 2.75, xend = 4.25, y = -70, yen = -70) +
    annotate("text", x = 3.6, y = - 90, label = "Sediment", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 4.75, xend = 6.25, y = -70, yen = -70) +
    annotate("text", x = 5.6, y = - 90, label = "Soilwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 6.75, xend = 8.25, y = -70, yen = -70) +
    annotate("text", x = 7.6, y = - 90, label = "Hyporheicwater", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 8.75, xend = 10.25, y = -70, yen = -70) +
    annotate("text", x = 9.6, y = - 90, label = "Stream", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 10.75, xend = 12.25, y = -70, yen = -70) +
    annotate("text", x = 11.6, y = - 90, label = "Tributary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 12.75, xend = 14.25, y = -70, yen = -70) +
    annotate("text", x = 13.6, y = - 90, label = "Headwater Lakes", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 14.75, xend = 16.25, y = -70, yen = -70) +
    annotate("text", x = 15.6, y = - 90, label = "Lake", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 16.75, xend = 18.25, y = -70, yen = -70) +
    annotate("text", x = 17.6, y = - 90, label = "Upriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 18.75, xend = 20.25, y = -70, yen = -70) +
    annotate("text", x = 19.6, y = - 90, label = "RO3", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 20.75, xend = 22.25, y = -70, yen = -70) +
    annotate("text", x = 21.6, y = - 90, label = "RO2", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 22.75, xend = 24.25, y = -70, yen = -70) +
    annotate("text", x = 23.6, y = - 90, label = "RO1", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 24.75, xend = 26.25, y = -70, yen = -70) +
    annotate("text", x = 25.6, y = - 90, label = "Hypolimnion", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 26.75, xend = 28.25, y = -70, yen = -70) +
    annotate("text", x = 27.6, y = - 90, label = "Downriver", angle = 45, vjust = 0.25, hjust = 1, size = 2.5) +
    annotate("segment", x = 28.75, xend = 30.25, y = -70, yen = -70) +
    annotate("text", x = 29.6, y = - 90, label = "Estuary", angle = 45, vjust = 0.25, hjust = 1, size = 2.5))

p <- ggarrange(univ.rare, spec, rare.shif, ab.shif,
               ncol = 2, nrow = 2, align = "v", common.legend = T, heights = c(0.45, 0.55),
               labels = "auto", legend = "right")
p


# save
ggsave(paste0("./Figures/Final/Abundance_profiles_heat.tiff"), p,
       width = 25, height = 15, unit = "cm")
ggsave(paste0("./Figures/Final/Abundance_profiles_heat.png"),  p,
       width = 25, height = 15, unit = "cm")

####################################################################################

# Plotting species points onto PCoA

(p <- ggplot() +
   theme_cust() +
   geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
   geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
   geom_point(data = all.pcoa$df,
              aes(x = Axis.1, y = Axis.2, 
                  shape = Season, alpha = DnaType),
              size = 2.5, fill = "grey40") +
   scale_shape_manual(values = c(21,23,25)) +
   scale_alpha_manual(values = c(1,0.5), name = "Nucleic Acid \nType") +
   geom_point(data = pb.species[ab.names == "universal.rare",], 
              aes(x = Axis.1, y = Axis.2, fill = ab.names), shape = 21) +
   geom_point(data = pb.species[ab.names != "universal.rare",], 
              aes(x = Axis.1, y = Axis.2, fill = ab.names), shape = 21) +
   #scale_fill_manual(values = colvec, name = "Habitat Type") +
   coord_fixed(1) + # ensure aspect ratio
   labs(x = paste("PC1 [", unique(all.pcoa$df$x),"%]"), 
        y = paste("PC2 [", unique(all.pcoa$df$y),"%]")) +
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
   guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
          alpha = guide_legend(order = 3, override.aes=list(size = 2)), 
          fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))
)

ggsave("./Figures/General/PCoA_all_ab.groups_withunivrare_pc12.png", p,
       height = 10, width = 15, units = "cm")

##############################################################################


# Ecologically ordered heatmap, colours correspond to species correlation with axes
cor.heat <- llply(tax.order, function(x){
  df <- cor.df[cor.df$OTU %in% x,]
  df$OTU <- factor(df$OTU, levels = x)
  
  # get max value of all
  max.lim <- round(max(cor.df$pear.cor),0)
  
  p <- ggplot(df, aes(x = Axes, y = OTU, fill = pear.cor)) +
    theme_bw() +
    geom_raster() +
    labs(x = "PCoA Axes", y = "ASVs") +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "red",
                         name = "Pearson's r") +
    scale_x_discrete(labels = c("PC1", "PC2", "PC3")) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.title.x = element_blank(),
          legend.position = "right", plot.title = element_text(size = 12)) +
    guides(fill = guide_colourbar(barwidth = unit(10, "pt"),
                                  barheight = unit(50, "pt")))
  
  return(p)
}, .parallel = T)

p <- ggarrange(cor.heat$universal.rare +labs(title="Universally rare"), 
               cor.heat$specialist +labs(title="Specialist"),
               cor.heat$rare.shifter +labs(title="Rare shifter"), 
               cor.heat$abundant.shifter+labs(title="Abundant shifter"),
               ncol = 2, nrow = 2, align = "v", common.legend = T,
               labels = "auto", legend = "right")
p

ggsave("./Figures/General/heatplot_ab.groups_PC_pearcor.png", p,
       height = 12, width = 12)
# PC1: negative = soil, positive = water
# PC2: negative = RNA, positive = DNA
# PC3: negative = summer/autumn, positive = spring