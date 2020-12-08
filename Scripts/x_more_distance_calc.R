# Calculate other distance metrics, depnding on need

pb.scores <- rbind(data.frame(Sample = as.character(row.names(pb.bray.pcoa$vectors)), Metric = "Bray",
                              pb.bray.pcoa$vectors, stringsAsFactors = F),
                   data.frame(Sample = as.character(row.names(pb.soren.pcoa$vectors)), Metric = "Sorensen",
                              pb.soren.pcoa$vectors, stringsAsFactors = F))# get first three axes
# merge with a selection of meta data
meta <- data.frame(Sample = as.character(row.names(sample_df(pb))),
                   sample_df(pb) %>% dplyr::select(sample.type.year, Season, Year, 
                                                   DnaType, distance.from.mouth, DR.names), 
                   stringsAsFactors = F)
data <- merge(pb.scores, meta, by = "Sample")
data$Sample <- as.character(data$Sample)

setDT(data)
setcolorder(data, c("Sample","Metric",
                    "sample.type.year","Season","Year", "DnaType","distance.from.mouth","DR.names",
                    colnames(data)[!(colnames(data) %in% c("Sample","Metric",
                                                           "sample.type.year","Season","Year", "DnaType","distance.from.mouth","DR.names"))]))



# Calculate 1D distance for a selection of axes
temp.1d <- melt(data, id.vars = c("DR.names","DnaType","Metric"), measure.vars = patterns("^Axis."),
                variable.name = "Axis", value.name = "Coordinates")

temp.1d <- dcast(temp.1d[Axis == "Axis.3" | Axis == "Axis.2" | Axis == "Axis.1",], 
                 DR.names + Axis + Metric ~ DnaType, value.var = "Coordinates")
temp.1d <- temp.1d[, .(dist = abs(DNA - RNA)), by = .(Metric,DR.names, Axis)]

# Calculate distance within 3D
temp.3d <- melt(data, id.vars = c("DR.names","DnaType","Metric"), measure.vars = patterns("^Axis."),
                variable.name = "Axis", value.name = "Coordinates")
temp.3d <- dcast(temp.3d[Axis == "Axis.3" | Axis == "Axis.2" | Axis == "Axis.1",], 
                 DR.names + Metric ~ DnaType + Axis, value.var = "Coordinates")

temp.3d[, dist.3d := sqrt(abs((DNA_Axis.1 - RNA_Axis.1))^2 + abs((DNA_Axis.2 - RNA_Axis.2))^2 + abs((DNA_Axis.3 - RNA_Axis.3))^2)]
#temp.2d[, distance.13 := sqrt(abs((DNA_Axis.1 - RNA_Axis.1))^2 + abs((DNA_Axis.3 - RNA_Axis.3))^2)]

# combine back with categories
dist.1d <- temp.1d[data, c("sample.type.year",
                           "Year", "Season") := list(i.sample.type.year,
                                                     i.Year, i.Season), on = .(DR.names)]

dist.3d <- temp.3d[data, c("sample.type.year",
                           "Year", "Season") := list(i.sample.type.year,
                                                     i.Year, i.Season), on = .(DR.names)]

# calculate confidence interval and means of sample type and season combinations
sum.dist1d <- dist.1d[, .(mean =  mean(dist, na.rm = T),
                          conf.int = conf.int(dist),
                          stdev = sd(dist, na.rm = T),
                          n = .N),
                      by = .(Metric, Axis, sample.type.year, Season)]

sum.dist3d <- dist.3d[, .(mean =  mean(dist.3d, na.rm = T),
                          conf.int = conf.int(dist.3d),
                          stdev = sd(dist.3d, na.rm = T),
                          n = .N),
                      by = .(Metric, sample.type.year, Season)]


# Correlate OTU to individual axes --------------------------------------------------------------

cor.mat <- cor(all.pcoa$df[,2:4], pb.mat, method = "pearson")
cor.mat <- setDT(as.data.frame(cor.mat), keep.rownames = "Axes")
cor.df <- melt.data.table(
  cor.mat,
  id.vars = "Axes",
  measure.vars = patterns("OTU_"),
  variable.name = "OTU",
  value.name = "pear.cor"
)

# does not work on private computer
#sim <- with(sample_df(pb), vegan::simper(t(otu_mat(pb)), DnaType, permutations = 1, parallel = cl))
#sim <- big.simper(otu_mat(pb), sample_df(pb)$DnaType, permutations = 100, parallel = 2)
#av.sim <- readRDS("./Objects/simper_orig_avg_mat.rds")

cor.df <- dcast(cor.df, OTU ~ Axes, value.var = "pear.cor")
colnames(cor.df)[2:4] <- c("cor.pc1","cor.pc2","cor.pc3")

cor.df <- cor.df[pb.species, c("cordi.pc1","cordi.pc2","cordi.pc3",
                               "ab.groups") := list(i.Axis.1, i.Axis.2, i.Axis.3,
                                                    i.ab.names), on = .(OTU)]

# remove OTUs only in RNA
cor.df <- cor.df[!is.na(cor.df$ab.groups),]


ggplot(cor.df, aes(x = cor.pc1, y = cor.pc2, fill = ab.groups)) +
  theme_bw() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
  geom_point(data = subset(cor.df, ab.groups == "universal.rare"), shape = 21, alpha = 0.5) +
  geom_point(data = subset(cor.df, ab.groups != "universal.rare"), shape = 21) +
  #geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"), shape = 21) +
  #geom_point(data = subset(cor.df, ab.groups == "rare.shifter"), shape = 21) +
  #geom_point(data = subset(cor.df, ab.groups == "specialist"), shape = 21) +
  scale_fill_viridis_d() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#-----------------------------------------------------------------------------------------------
# compare slopes

cor.df$ab.groups <- factor(cor.df$ab.groups, levels = c('abundant.shifter', 'present.as', 
                                                        'rare.shifter', 'present.rs', 
                                                        'specialist', 'universal.rare'),
                           labels = c('Abundant shifter', 'Universal \nabundant shifter',
                                      'Rare shifter', 'Universal \nrare shifter',
                                      'Specialist', 'Universal rare'))

okabe.ito.sel <- c("#D55E00", "#CC79A7",  "#0072B2", "#56B4E9", "#009E73", "#E69F00")

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("","Terrestrial"
                   ,"","Water"),
  hjustvar = c(0,-0.08,1,1.2) ,
  vjustvar = c(0,1.2,0,1.2))

#hjust: lower values right, larger values left
#vjust lower values up, higher values down

#Bottom Left (h0,v0)","Top Left (h0,v1)"
#,"Bottom Right h1,v0","Top Right h1,v1"

(p1 <- ggplot(shape = 19) +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "Universal rare"), 
               aes(x = cordi.pc1, y = abs(cor.pc1), colour = ab.groups),
               shape = 21, alpha = 0.2, size = 0.5, colour = "grey80") +
    geom_point(data = subset(cor.df, ab.groups != "Universal rare"),
               aes(x = cordi.pc1, y = abs(cor.pc1), colour = ab.groups), shape = 21, alpha = 0.2, size = 0.5) +
    #geom_smooth(data = cor.df[ab.groups != "present.rs",], 
    #              aes(x = cordi.pc1, y = abs(cor.pc1), group = ab.groups,
    #                   colour = ab.groups), size = 1.5,
    #              method = 'gam', formula = y ~ s(x, bs = "cs"), se = F) +
    geom_smooth(data = cor.df[cordi.pc1 > 0,], 
                aes(x = cordi.pc1, y = abs(cor.pc1),
                    colour = ab.groups), size = 1.2,
                method = 'lm', formula = y ~ x, se = F) +
    geom_smooth(data = cor.df[cordi.pc1 < 0,], 
                aes(x = cordi.pc1, y = abs(cor.pc1),
                    colour = ab.groups), size = 1.2,
                method = 'lm', formula = y ~ x, se = F) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText), size = 3) +
    scale_y_continuous(limits = c(0, 0.75)) +
    scale_colour_manual(values = okabe.ito.sel,
                        name = "Abundance \ngroup",
                        breaks = c('Abundant shifter', 'Universal \nabundant shifter',
                                   'Rare shifter', 'Universal \nrare shifter',
                                   'Specialist', 'Universal rare')) +
    labs(x = "Species coordinates [PC1]", 
         y = expression(paste(" | Pearsons's ", italic("r")," | of OTUs with PC1"))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title.align=0.5, legend.text.align = 0.5)
)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("","RNA"
                   ,"","DNA"),
  hjustvar = c(0,-0.08,1,1.2) ,
  vjustvar = c(0,1.2,0,1.2))

(p2 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "Universal rare"), 
               aes(x = cordi.pc2, y = abs(cor.pc2), colour = ab.groups),
               shape = 21, alpha = 0.2, size = 0.5, colour = "grey80") +
    geom_point(data = subset(cor.df, ab.groups != "Universal rare"),
               aes(x = cordi.pc2, y = abs(cor.pc2), colour = ab.groups), shape = 21, alpha = 0.2, size = 0.5) +
    #geom_smooth(data = cor.df[ab.groups != "present.rs",], 
    #              aes(x = cordi.pc2, y = abs(cor.pc2), group = ab.groups,
    #                   colour = ab.groups), size = 1.5,
    #              method = 'gam', formula = y ~ s(x, bs = "cs"), se = F) +
    geom_smooth(data = cor.df[cordi.pc2 > 0,], 
                aes(x = cordi.pc2, y = abs(cor.pc2),
                    colour = ab.groups), size = 1.2,
                method = 'lm', formula = y ~ x, se = F) +
    geom_smooth(data = cor.df[cordi.pc2 < 0,], 
                aes(x = cordi.pc2, y = abs(cor.pc2),
                    colour = ab.groups), size = 1.2,
                method = 'lm', formula = y ~ x, se = F) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText), size = 3) +
    scale_y_continuous(limits = c(0, 0.75)) +
    scale_colour_manual(values = okabe.ito.sel,
                        name = "Abundance \ngroup",
                        breaks = c('Abundant shifter', 'Universal \nabundant shifter',
                                   'Rare shifter', 'Universal \nrare shifter',
                                   'Specialist', 'Universal rare')) +
    labs(x = "Species coordinates [PC2]", 
         y = expression(paste(" | Pearsons's ", italic("r")," | of OTUs with PC2"))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title.align=0.5, legend.text.align = 0.5)
)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("","Summer/Autumn"
                   ,"","Spring"),
  hjustvar = c(0,-0.05,1,1.2) ,
  vjustvar = c(0,1.2,0,1.2))

(p3 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "Universal rare"), 
               aes(x = cordi.pc3, y = abs(cor.pc3), colour = ab.groups),
               shape = 21, alpha = 0.2, size = 0.5, colour = "grey80") +
    geom_point(data = subset(cor.df, ab.groups != "Universal rare"),
               aes(x = cordi.pc3, y = abs(cor.pc3), colour = ab.groups), shape = 21, alpha = 0.2, size = 0.5) +
    #geom_smooth(data = cor.df[ab.groups != "present.rs",], 
    #              aes(x = cordi.pc3, y = abs(cor.pc3), group = ab.groups,
    #                   colour = ab.groups), size = 1.5,
    #              method = 'gam', formula = y ~ s(x, bs = "cs"), se = F) +
    geom_smooth(data = cor.df[cordi.pc3 < 0,], 
                aes(x = cordi.pc3, y = abs(cor.pc3), group = ab.groups,
                    colour = ab.groups), size = 1.2,
                method = 'lm', formula = y ~ x, se = F) +
    geom_smooth(data = cor.df[cordi.pc3 > 0,], 
                aes(x = cordi.pc3, y = abs(cor.pc3),
                    colour = ab.groups), size = 1.2,
                method = 'lm', formula = y ~ x, se = F) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText), size = 3) +
    scale_y_continuous(limits = c(0, 0.75)) +
    scale_colour_manual(values = okabe.ito.sel,
                        name = "Abundance \ngroup",
                        breaks = c('Abundant shifter', 'Universal \nabundant shifter',
                                   'Rare shifter', 'Universal \nrare shifter',
                                   'Specialist', 'Universal rare')) +
    labs(x = "Species coordinates [PC3]", 
         y = expression(paste(" | Pearsons's ", italic("r")," | of OTUs with PC3"))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title.align=0.5, legend.text.align = 0.5)
)

# combine
(p.trio <- ggarrange(p1, p2, p3, ncol = 3, common.legend = T))

ggsave("./Figures/Final/abgroups_cor.png", p.trio, width = 20, height = 9, units = "cm")

#-----------------------------------------------------------------------------------------------

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Terrestrial - RNA","Terrestrial - DNA"
                   ,"Water - RNA","Water - DNA"),
  hjustvar = c(-0.05, -0.05, 
               1.05, 1.05) ,
  vjustvar = c(-0.2, 1.2, 
               -0.2, 1.2))

(pp1 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.2) +
    geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"), 
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "rare.shifter"),
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "specialist"), 
               aes(x = cor.pc1, y = cor.pc2, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    labs(x = expression(paste("Pearsons's ", italic("r")," of ASVs with PC1")), 
         y = expression(paste("Pearsons's ", italic("r")," of ASVs with PC2"))) +
    scale_fill_viridis_d() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Terrestrial - Summer","Terrestrial - Spring"
                   ,"Water - Summer","Water - Spring"),
  hjustvar = c(-0.05, -0.05, 
               1.05, 1.05) ,
  vjustvar = c(-0.2, 1.2, 
               -0.2, 1.2))

(pp2 <- ggplot() +
    theme_bw() +
    geom_hline(yintercept =  0, colour = "grey80", size = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", size = 0.6) +
    geom_point(data = subset(cor.df, ab.groups == "universal.rare"), 
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.2) +
    geom_point(data = subset(cor.df, ab.groups == "abundant.shifter"), 
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "rare.shifter"),
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.8) +
    geom_point(data = subset(cor.df, ab.groups == "specialist"), 
               aes(x = cor.pc1, y = cor.pc3, fill = ab.groups), shape = 21, alpha = 0.8) +
    labs(x = expression(paste("Pearsons's ", italic("r")," of ASVs with PC1")), 
         y = expression(paste("Pearsons's ", italic("r")," of ASVs with PC3"))) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,
                                   hjust=hjustvar,vjust=vjustvar,
                                   label=annotateText)) +
    scale_fill_viridis_d() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)

p.cordcor <- ggarrange(p1, p2, p3, ncol = 3, common.legend = T, align = "hv")

p.corcor <- ggarrange(pp1, pp2, ncol = 2, common.legend = T, align = "hv")

ggsave("./Figures/Final/abundance_grps_pear.cor_cord.png", p.cordcor,
       height = 12, width = 30, units = "cm")

ggsave("./Figures/Final/abundance_grps_pear.cors_pcs.png", p.corcor,
       height = 10, width = 20, units = "cm")
