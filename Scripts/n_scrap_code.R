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
pb.mat <- otu_mat(pb)
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
bd.df$ID.x[bd.df$DnaType.x == "RNA"] <- str_replace(bd.df$Sample.x[bd.df$DnaType.x == "RNA"], "R$", "")
bd.df$ID.y[bd.df$DnaType.y == "DNA"] <- str_replace(bd.df$Sample.y[bd.df$DnaType.y == "DNA"], "D$", "")
bd.df$ID.y[bd.df$DnaType.y == "RNA"] <- str_replace(bd.df$Sample.y[bd.df$DnaType.y == "RNA"], "R$", "")

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
                                             "Soilwater", "Stream", "Tributary",
                                             "Riverine \nLakes", "Lake",
                                             "Upriver","Downriver",
                                             "Reservoirs",
                                             "Estuary"))

library(ggtern)
(terti <- ggtern(bd.dt[Metric == "Jaccard",], aes(x = I, y = RC, z = S)) +
  theme_bw() +
  geom_point(aes(fill = sample.type.year), shape = 21, size = 2, alpha =.8) +
  theme_showarrows() +
  facet_grid(.~Season) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  labs(x = "", y = "", z = "",
       xarrow = "I", yarrow = "RC", zarrow = "S"))

ggsave("./Figures/Final/indicandum_tertiary_plot.png",
       height = 10, width = 20, units = "cm")

# I = intersection of nestedness and beta div = richness difference
# RC = Relative complement of nestedness in beta diversity = replacement
# S = Sorensen

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


###
# DNA
###

# Execute regressions first on DNA diversity
df <- melt(dna.alpha, id.vars = c("DR.names","Data"),
           measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
           variable.name = "Index",
           value.name = "Diversity")

# Only focus on CSS, non rarefied data and two diversity indices
df <- df[Data == "css" & (Index == "Shannon" | Index == "Pielou"),]

# merge with distance
reg.df <- df[dist.75[Metric == "Sorensen",], c("distance", "sample.type.year", "Season") := 
               list(i.dist, i.sample.type.year, i.Season), on = .(DR.names)]

# remove any NAs, e.g. samples from 2015
reg.df <- reg.df[!is.na(distance),]
setorderv(reg.df, c("Index","Diversity")) # rearrange dataframe

# decide which model is best (e.g. linear or polynomial)
z <- reg.df[Index == "Shannon"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
           distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 -0.03645 p > 0.05
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.47 p > 0.05
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.18 p > 0.05
anova(lm0,lm1) # preferred model is lm1, higher R2 and lowered RSS
# polynomials do not add much, and avoid overfitting and choose a parsimonious model
anova(lm1,lm2) # preferred model lm2
anova(lm2,lm3) 
rm(z)


z <- reg.df[Index == "Pielou"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
           distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 0.70
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.85
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.81
anova(lm0,lm1) 
anova(lm1,lm2) # preferred model is lm2, higher R2 and lowered RSS, and significance < 0.01
anova(lm2,lm3) 
rm(z)

lin.ls <- dlply(reg.df, .(Index), function(z){
  setDT(z)
  means <- z[, .(Diversity = mean(Diversity, na.rm = T),
                 distance = mean(distance, na.rm = T)), by = .(sample.type.year)]
  if(unique(z$Index) == "Shannon"){
    lin <- lm(means$Diversity ~ poly(means$distance, 2))
  } else if(unique(z$Index) == "Pielou") {
    lin <-  lm(means$Diversity ~ poly(means$distance, 2))
  }
  
  # check linear assumptions
  #plot(lin) # normality not good
  # large sample sizes, normality does not affect results too much (central limit theorem)
  # homoscedasticity and independence important
  #summary(lin)
  #confint(lin, level = 0.95)
  
  # get data for plotting
  x <- data.frame(x = sort(means$distance))
  pred <- predict(lin, newdata = x, se = T)
  ci <- pred$se.fit[order(means$distance)] * qt(0.95 / 2 + 0.5, pred$df)
  y <- pred$fit[order(means$distance)]
  ymin <- y - ci
  ymax <- y + ci
  
  plot.df <- data.frame(x = sort(means$distance),
                        y = y,
                        ymin = ymin,
                        ymax = ymax,
                        se = pred$se.fit[order(means$distance)])
  
  # extract only colours that are in data frame
  colvec <- colvec[names(colvec) %in% as.character(levels(means$sample.type.year))]
  
  p <- ggplot() +
    theme_cust("pubr") +
    geom_point(data = z, aes(x = distance, y = Diversity), 
               colour = "gray40", alpha = 0.3, size = 2) +
    geom_line(data = plot.df, aes(x = x, y = y), inherit.aes = F, size = 1) +
    geom_line(data = plot.df, aes(x = x, y = ymax), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_line(data = plot.df, aes(x = x, y = ymin), inherit.aes = F, size = 0.7, linetype = "dashed") +
    geom_point(data = means,
               aes(x = distance, y = Diversity, fill = sample.type.year), shape = 21, size = 3) +
    scale_fill_manual(values = colvec, name = "Sample Type") +
    labs(x = "",
         y = paste0(unique(z$Index)))
  
  # get model statistics
  options(scipen = 999) # avoid scientific annotations
  
  fnr <- substitute(italic(R)^2~"="~r2*","~~italic(F)[df]~"="~Fstat,
                    list(r2 = format(summary(lin)$r.squared, digits = 2),
                         Fstat = format(summary(lin)$fstatistic[[1]], digits = 4),
                         df = paste0(format(summary(lin)$fstatistic[[2]], digits = 0),
                                     ",", format(summary(lin)$fstatistic[[3]], digits = 0))))
  pv1 <- summary(lin)$coefficients[2,4]
  pv1 <- if(pv1 < 0.0001){
    "< 0.0001"} else if(pv1 < 0.001){
      "< 0.001"} else if(pv1 < 0.01){
        "< 0.01"} else if(pv1 < 0.05){
          "< 0.05"
        } else {
          paste("=",round(pv1, 2))
        }
  
  if(unique(z$Index) == "Pielou"){
    eq1 <- substitute(italic(y) == a - b %.% italic(x) + b2 %.% italic(x)^2,
                      list(a = format(as.vector(coef(lin)[1]), digits = 2),
                           b = format(as.vector(abs(coef(lin)[2])), digits = 2),
                           b2 = format(as.vector(coef(lin)[3]), digits = 2)))
    
    pv2 <- summary(lin)$coefficients[3,4]
    pv2 <- if(pv2 < 0.0001){
      "< 0.0001"} else if(pv2 < 0.001){
        "< 0.001"} else if(pv2 < 0.01){
          "< 0.01"} else if(pv2 < 0.05){
            "< 0.05"
          } else {
            paste("=",round(pv2, 2))
          }
    ps <- substitute(italic(p)[beta[1]]~pval1*","~italic(p)[beta[2]]~pval2,
                     list(pval1 = pv1,
                          pval2 = pv2))
    (p <- p + 
        annotate("text", x = 0.4, y = 0.9, 
                 label = as.character(as.expression(eq1)), parse = T, size = 2.5)+
        annotate("text", x = 0.4, y = 0.88, 
                 label = as.character(as.expression(fnr)), parse = T, size = 2.5) +
        annotate("text", x = 0.4, y = 0.86, 
                 label = as.character(as.expression(ps)), parse = T, size = 2.5) 
    )
  } else {
    eq1 <- substitute(italic(y) == a - b %.% italic(x),
                      list(a = format(as.vector(coef(lin)[1]), digits = 2),
                           b = format(as.vector(abs(coef(lin)[2])), digits = 2)))
    ps <- substitute(italic(p)~pval1,
                     list(pval1 = pv1))
    
    (p <- p + 
        annotate("text", x = 0.4, y = 6.5, 
                 label = as.character(as.expression(eq1)), parse = T, size = 2.5)+
        annotate("text", x = 0.4, y = 6.3, 
                 label = as.character(as.expression(fnr)), parse = T, size = 2.5) +
        annotate("text", x = 0.4, y = 6.1, 
                 label = as.character(as.expression(ps)), parse = T, size = 2.5) 
    )
  }
  
  #, abs(round(coef(lin)[2], 2)), "*x +",
  #round(coef(lin)[3], 2), "*x"^2*""
  
  list(original = z,
       binned = means,
       lin = lin,
       coef = coef(lin),
       fitted = plot.df,
       plot = p)
  
})

(p <- ggarrange(lin.ls[[1]]$plot, lin.ls[[2]]$plot,
                ncol = 2, common.legend = T, legend = "right"))

(p <- annotate_figure(p, 
                      bottom = text_grob("Distance in ordination space",
                                         just = "centre")))
ggsave("./Figures/Final/Richness_distance_nonlin_reg.png", p,
       width = 22, height = 11, unit = "cm")


####
## RNA
####

# Execute regressions now with RNA diversity
df <- melt(rna.alpha, id.vars = c("DR.names","Data"),
           measure.vars = c("Shannon","Simpson","Pielou","Chao1"),
           variable.name = "Index",
           value.name = "Diversity")

# Only focus on CSS, non rarefied data and two diversity indices
df <- df[Data == "css" & (Index == "Shannon" | Index == "Pielou"),]

# calculate mean of duplicates
df[, .(Diversity = mean(Diversity, na.rm = T)), by = .(DR.names, Index)]

# merge with distance
reg.df <- df[dist.75, c("distance", "sample.type.year") := 
               list(i.dist, i.sample.type.year), on = .(DR.names)]

# remove any NAs, e.g. samples from 2015
reg.df <- reg.df[!is.na(distance),]
setorderv(reg.df, c("Index","Diversity")) # rearrange dataframe

# decide which model is best (e.g. linear or polynomial)
z <- reg.df[Index == "Shannon"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
           distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 0.018
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.33
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.34
anova(lm0,lm1)
anova(lm1,lm2) # preferred model is lm2, higher R2 and lowered RSS
# linear model itself is very bad
anova(lm2,lm3) 
rm(z)


z <- reg.df[Index == "Pielou"]
z <- z[, .(Diversity = mean(Diversity, na.rm = T),
           distance = mean(distance, na.rm = T)), by = .(sample.type.year)]

plot(Diversity ~ distance, data = z) # does not seem linear
lm0 <- lm(z$Diversity ~ 1)  ; summary(lm0)
lm1 <- lm(z$Diversity ~ z$distance); summary(lm1) # adj R2 0.05
lm2 <- lm(z$Diversity ~ poly(z$distance, 2)); summary(lm2) # adj R2 0.53
lm3 <- lm(z$Diversity ~ poly(z$distance, 3)); summary(lm3) # adj R2 0.61
anova(lm0,lm1) 
anova(lm1,lm2)
anova(lm2,lm3) # preferred model is lm2, higher R2 and lowered RSS
# poly 3 has higher R2 but is not significant
rm(z)

lin.ls <- dlply(reg.df, .(Index), function(z){
  setDT(z)
  means <- z[, .(Diversity = mean(Diversity, na.rm = T),
                 distance = mean(distance, na.rm = T)), by = .(sample.type.year)]
  # both are poly 2
  lin <-  lm(means$Diversity ~ means$distance, 2))

# check linear assumptions
#plot(lin) # normality not good
# large sample sizes, normality does not affect results too much (central limit theorem)
# homoscedasticity and independence important
#summary(lin)
#confint(lin, level = 0.95)

# get data for plotting
x <- data.frame(x = sort(means$distance))
pred <- predict(lin, newdata = x, se = T)
ci <- pred$se.fit[order(means$distance)] * qt(0.95 / 2 + 0.5, pred$df)
y <- pred$fit[order(means$distance)]
ymin <- y - ci
ymax <- y + ci

plot.df <- data.frame(x = sort(means$distance),
                      y = y,
                      ymin = ymin,
                      ymax = ymax,
                      se = pred$se.fit[order(means$distance)])

# extract only colours that are in data frame
colvec <- colvec[names(colvec) %in% as.character(levels(means$sample.type.year))]

p <- ggplot() +
  theme_cust("pubr") +
  geom_point(data = z, aes(x = distance, y = Diversity), 
             colour = "gray40", alpha = 0.3, size = 2) +
  geom_line(data = plot.df, aes(x = x, y = y), inherit.aes = F, size = 1) +
  geom_line(data = plot.df, aes(x = x, y = ymax), inherit.aes = F, size = 0.7, linetype = "dashed") +
  geom_line(data = plot.df, aes(x = x, y = ymin), inherit.aes = F, size = 0.7, linetype = "dashed") +
  geom_point(data = means,
             aes(x = distance, y = Diversity, fill = sample.type.year), shape = 21, size = 3) +
  scale_fill_manual(values = colvec, name = "Sample Type") +
  labs(x = "",
       y = paste0(unique(z$Index)))

# get model statistics
options(scipen = 999) # avoid scientific annotations

fnr <- substitute(italic(R)^2~"="~r2*","~~italic(F)[df]~"="~Fstat,
                  list(r2 = format(summary(lin)$r.squared, digits = 2),
                       Fstat = format(summary(lin)$fstatistic[[1]], digits = 4),
                       df = paste0(format(summary(lin)$fstatistic[[2]], digits = 0),
                                   ",", format(summary(lin)$fstatistic[[3]], digits = 0))))
pv1 <- summary(lin)$coefficients[2,4]
pv1 <- if(pv1 < 0.0001){
  "< 0.0001"} else if(pv1 < 0.001){
    "< 0.001"} else if(pv1 < 0.01){
      "< 0.01"} else if(pv1 < 0.05){
        "< 0.05"
      } else {
        paste("=",round(pv1, 2))
      }
if(unique(z$Index) == "Pielou"){
  eq1 <- substitute(italic(y) == a - b %.% italic(x) + b2 %.% italic(x)^2,
                    list(a = format(as.vector(coef(lin)[1]), digits = 2),
                         b = format(as.vector(abs(coef(lin)[2])), digits = 2),
                         b2 = format(as.vector(coef(lin)[3]), digits = 2)))
  
  pv2 <- summary(lin)$coefficients[3,4]
  pv2 <- if(pv2 < 0.0001){
    "< 0.0001"} else if(pv2 < 0.001){
      "< 0.001"} else if(pv2 < 0.01){
        "< 0.01"} else if(pv2 < 0.05){
          "< 0.05"
        } else {
          paste("=",round(pv2, 2))
        }
  ps <- substitute(italic(p)[beta[1]]~pval1*","~italic(p)[beta[2]]~pval2,
                   list(pval1 = pv1,
                        pval2 = pv2))
  (p <- p + 
      annotate("text", x = 0.4, y = 0.9, 
               label = as.character(as.expression(eq1)), parse = T, size = 2.5)+
      annotate("text", x = 0.4, y = 0.88, 
               label = as.character(as.expression(fnr)), parse = T, size = 2.5) +
      annotate("text", x = 0.4, y = 0.86, 
               label = as.character(as.expression(ps)), parse = T, size = 2.5) 
  )
} else {
  eq1 <- substitute(italic(y) == a - b %.% italic(x) + b2 %.% italic(x)^2,
                    list(a = format(as.vector(coef(lin)[1]), digits = 2),
                         b = format(as.vector(abs(coef(lin)[2])), digits = 2),
                         b2 = format(as.vector(coef(lin)[3]), digits = 2)))
  
  pv2 <- summary(lin)$coefficients[3,4]
  pv2 <- if(pv2 < 0.0001){
    "< 0.0001"} else if(pv2 < 0.001){
      "< 0.001"} else if(pv2 < 0.01){
        "< 0.01"} else if(pv2 < 0.05){
          "< 0.05"
        } else {
          paste("=",round(pv2, 2))
        }
  ps <- substitute(italic(p)[beta[1]]~pval1*","~italic(p)[beta[2]]~pval2,
                   list(pval1 = pv1,
                        pval2 = pv2))
  
  (p <- p + 
      annotate("text", x = 0.4, y = 6.5, 
               label = as.character(as.expression(eq1)), parse = T, size = 2.5)+
      annotate("text", x = 0.4, y = 6.3, 
               label = as.character(as.expression(fnr)), parse = T, size = 2.5) +
      annotate("text", x = 0.4, y = 6.1, 
               label = as.character(as.expression(ps)), parse = T, size = 2.5) 
  )
}

#, abs(round(coef(lin)[2], 2)), "*x +",
#round(coef(lin)[3], 2), "*x"^2*""

list(original = z,
     binned = means,
     lin = lin,
     coef = coef(lin),
     fitted = plot.df,
     plot = p)

})

(p <- ggarrange(lin.ls[[1]]$plot, lin.ls[[2]]$plot,
                ncol = 2, common.legend = T, legend = "right"))

(p <- annotate_figure(p, 
                      bottom = text_grob("Distance in ordination space",
                                         just = "centre")))
ggsave

##
set.seed(3)
t <- sample(1:nrow(pb.mat),10)
submat <- pb.mat[t,]

submat <- submat[,colSums(submat) > 0]

meta <- data.frame(Sample = as.character(row.names(sample_df(dna))),
                   sample_df(dna) %>% dplyr::select(sample.type.year, Season, Year, DnaType), 
                   stringsAsFactors = F)
meta <- meta[t,]

# melt to calculate mean variance relationship
melt.mat <- melt.data.table(
  setDT(as.data.frame(pb.mat), keep.rownames = "Sample"),
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

plot.df <- melt.mat[, .(mean = mean(reads, na.rm = T),
                        variance = var(reads, na.rm = T)), by = .(OTU)]

ggplot(plot.df, aes(x = log1p(mean), y = log(variance))) +
  geom_point()

# ASVs with high means also have high variances

#ord.asv <- plot.df$ASV[order(plot.df$mean, decreasing = T)]
#melt.mat$ASV <- factor(melt.mat$ASV, levels = ord.asv)
#melt.mat <- melt.mat[meta, c("sample.type.year", "Season") := list(i.sample.type.year,
#                                                               i.Season), on = .(Sample)]
#ggplot(melt.mat, aes(x = ASV, y = log2(reads + 1), colour = sample.type.year)) +
#  geom_point()


# make mvabund object of community matrix
dna.sp <- mvabund(pb.mat)
dna.sp <- mvabund(submat)

mod <- manyglm(dna.sp ~ meta$sample.type.year * meta$Season, family = "negative.binomial")
# warning but is integer
saveRDS(mod, "./Objects/manyglm.dna.negbinom.log.rds")

# check residuals, it's not optimal, but compared to other families, there is less of a pattern
png(filename="./Figures/General/manyglm_dna_residuals_binom_log.png")
plot(mod)
dev.off()

# test for habitat type and season effect
anova.mod <- anova(mod)
saveRDS(anova.mod, "./Objects/manyglm.dna.negbinom.anova.rds")
print("DONE")
pb.mat <- decostand(pb.mat, "hellinger")
pb.mori <- vegdist(pb.mat, method = "horn")
is.euclid(pb.mori) # FALSE
pb.mori <- sqrt(pb.mori) # make Euclidean
is.euclid(pb.mori) # TRUE
anova <- anova(mod)

## Run bayesian ordination
# test control options, for quick building. Not final
mcmc.control. <- list(n.burnin = 10, 
                      n.iteration = 400, 
                      n.thin = 30, 
                      seed = 3)

fit.lvmbinom <- boral(y = pb.mat, 
                      family = "negative.binomial", 
                      num.lv = 2, 
                      mcmc.control = mcmc.control.,
                      row.eff = "fixed")

####################################################################



dissim.dr <- dissim.dnarna(ter, save.name = "All", output = T)

dist.dr <- dist.dnarna(dna.pcoa[["df"]], save.name = "terr", output = T)
#dist.dr <- dist.dnarna(dnarna.bray[["df"]], save.name = "All_2D", dimensions = 2)

test <- merge(dissim.dr$original.df, dist.dr$indiv.df, by.x = "ID", by.y = "DR.names")
#dist.dr$df[Axis == "Axis.2",]
ggplot(test, aes(x = dist, y = distance.1D, fill = sample.type.year.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_manual(values = colvec) +
  labs(x = "Pair-wise Bray-Curtis dissimilarity", y = "Pair-wise distance in ordination")

p <- ggplot(test, aes(x = dist, y = distance.1D, fill = Season.x)) +
  theme_bw() +
  geom_point(shape = 21) +
  labs(x = "Pair-wise Bray-Curtis dissimilarity", y = "Pair-wise distance in ordination")

ggsave("./Figures/General/terr_distdissimcor.png", p)
cor.test(test$dist, test$distance.1D)