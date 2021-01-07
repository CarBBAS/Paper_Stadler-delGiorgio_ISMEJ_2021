# How do we interpret the mS:mBC ratio? --------------------------------------------------------
# extract sample pairs above 1 and below 1 of the mS:mBC ratio
big <- as.character(diff.df[bcs.ratio > 1,]$DR.names)
small <- as.character(diff.df[bcs.ratio < 1,]$DR.names)

# melt community matrix for tidy data set
commat <- melt.data.table(
  setDT(as.data.frame(otu_mat(pb)), keep.rownames = "Sample"),
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

# add meta variables
commat[setDT(sample_df(pb), keep.rownames = "Sample"), 
       c("DR.names", "DnaType") := list(i.DR.names, i.DnaType), on = .(Sample)]

# keep only the ones that have both DNA and RNA
#commat <- commat[duplicated(DR.names),]
# doesn't make a difference

# cast so that we can calculate difference between DNA and RNA of each sample
temp <- dcast(commat, DR.names + OTU ~ DnaType, value.var = c("reads"))
temp[, diff := DNA - RNA]
temp <- temp[!is.na(diff),]

# calculate summary variables
abs.diff <- temp[, .(mean.diff = mean(abs(diff), na.rm = T),
                     max.diff = max(abs(diff), na.rm = T),
                     var.diff = var(abs(diff), na.rm = T)), by = .(DR.names)]

# calculate abundance difference of shared taxa
shar.diff <- temp[DNA > 0 & RNA > 0,]
shar.diff <- shar.diff[, diff := abs(DNA - RNA)][, .(shar.diff = mean(diff)), by = .(DR.names)]

# calculate how many taxa are shared by sample
n.shared <- temp[DNA > 0 & RNA > 0, ]
n.shared <- n.shared[, .(n = .N), by = .(DR.names)]

# calculate richness difference
dna.rich <- temp[DNA > 0,][, .(n = .N), by = .(DR.names)]
rna.rich <- temp[RNA > 0,][, .(n = .N), by = .(DR.names)]
rich.diff <- dna.rich[rna.rich, n.rna := i.n, on = .(DR.names)][, rich.diff := n - n.rna]

# calculate the abundance difference mean of the OTUs that have been classified as abundant
abun.diff <- temp[DNA >= 47,]
abun.diff <- abun.diff[, .(mean.abun.diff = mean(abs(diff))), by = .(DR.names)]

# extract distances
dists <- dcast(dist.75, DR.names ~ Metric, value.var = "dist")

# calculate replacement
# binary transformation
pa.temp <- setDF(temp); setDT(pa.temp)
pa.temp <- pa.temp[, c("DNA",
                       "RNA") := list(ifelse(DNA > 0, 1, 0),
                                      ifelse(RNA > 0, 1, 0))]
replac <- pa.temp[DNA != RNA,][, .(n = .N), by = .(DR.names)]
ab.replac <- temp[which(pa.temp$DNA != pa.temp$RNA),][, replac.diff := abs(DNA - RNA)][, .(replac.diff = mean(replac.diff)), by = .(DR.names)]

# combine all to one data set
all.met <- abs.diff[n.shared, n.shared := i.n, on = .(DR.names)]
all.met <- all.met[rich.diff, rich.diff := abs(i.rich.diff), on = .(DR.names)]
all.met <- all.met[replac, replac := i.n, on = .(DR.names)]
all.met <- all.met[ab.replac, replac.diff := i.replac.diff, on = .(DR.names)]
all.met <- all.met[shar.diff, on = .(DR.names)]
all.met <- all.met[dists, c("mBC", "mS") := list(i.Bray, i.Sorensen), on = .(DR.names)]
all.met <- all.met[diff.df, ratio := i.bcs.ratio, on = .(DR.names)]

# give >1 and <1 mS:mBC categories
all.met[DR.names %in% big, ratio.cat := ">1"]
all.met[DR.names %in% small, ratio.cat := "<1"]
#all.met <- all.met[!is.na(ratio.cat),]

# add categories
all.met[setDT(sample_df(pb) %>%
                select(DR.names, sample.type.year, Season) %>%
                group_by(DR.names) %>% distinct()), c("sample.type.year","Season") := list(i.sample.type.year,
                                                                                           i.Season),
        on = .(DR.names)]

# Plot regressions
# mBC... with
# Shared n
# Make model
model <- lm(all.met$n.shared ~ all.met$mBC)
model.df <- predict.lm(model, interval = "confidence") %>%
  cbind(all.met %>% select(mBC))

# Shared taxa
(mbc.n <- ggplot(all.met, aes(x = mBC, y = n.shared)) +
    geom_point(aes(fill = sample.type.year, shape = Season),
               size = 3) +
    geom_line(data = model.df, aes(x= mBC, y = fit), colour = "black") +
    geom_line(data = model.df, aes(x= mBC, y = upr), colour = "black", linetype = "dashed") +
    geom_line(data = model.df, aes(x= mBC, y = lwr), colour = "black", linetype = "dashed") +
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
                 label.y = "top", label.x = "right",
                 parse = TRUE) +
    labs(x = expression(paste(italic("m")["BC"])), y = "Number of shared taxa between DNA-RNA") +
    scale_shape_manual(values = c(21,23,25)) +
    scale_fill_manual(values = colvec, name = "Habitat Type") +
    guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
           fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2))))

# mBC... with
# Abundance difference
# Make model
(outliers <- all.met %>%
    rstatix::identify_outliers(mean.diff) %>%
    select(DR.names, is.outlier, is.extreme))

# Exclude extreme outliers
plot.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

model <- lm(plot.df$mean.diff ~ plot.df$mBC)
model.df <- predict.lm(model, interval = "confidence") %>%
  cbind(plot.df %>% select(mBC))

# Shared taxa
mbc.diff <- ggplot(plot.df, aes(x = mBC, y = mean.diff)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_line(data = model.df, aes(x= mBC, y = fit), colour = "black") +
  geom_line(data = model.df, aes(x= mBC, y = upr), colour = "black", linetype = "dashed") +
  geom_line(data = model.df, aes(x= mBC, y = lwr), colour = "black", linetype = "dashed") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["BC"])),
       y = expression(paste("Mean | ", Delta, " OTU CSS reads |"))) +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

# mS... with
# Shared n
# Make model
model <- lm(all.met$n.shared ~ all.met$mS)
model.df <- predict.lm(model, interval = "confidence") %>%
  cbind(all.met %>% select(mS))

# Shared taxa
ms.n <- ggplot(all.met, aes(x = mS, y = n.shared)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_line(data = model.df, aes(x= mS, y = fit), colour = "black") +
  geom_line(data = model.df, aes(x= mS, y = upr), colour = "black", linetype = "dashed") +
  geom_line(data = model.df, aes(x= mS, y = lwr), colour = "black", linetype = "dashed") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["S"])), y = "Number of shared taxa between DNA-RNA") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))


(outliers <- all.met %>%
    rstatix::identify_outliers(mean.diff) %>%
    select(DR.names, is.outlier, is.extreme))

# Exclude extreme outliers
plot.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

ggplot(plot.df, aes(x = n.shared, y = mean.diff)) +
  geom_point(aes(fill = ratio.cat, shape = Season),
             size = 3) +
  #geom_line(data = model.df, aes(x= mS, y = fit), colour = "black") +
  #geom_line(data = model.df, aes(x= mS, y = upr), colour = "black", linetype = "dashed") +
  #geom_line(data = model.df, aes(x= mS, y = lwr), colour = "black", linetype = "dashed") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  #labs(x = expression(paste(italic("m")["S"])), y = "Number of shared taxa between DNA-RNA") +
  scale_shape_manual(values = c(21,23,25)) +
  #scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

ggplot(all.met, aes(x = mS, y = mBC)) +
  geom_point(aes(fill = ratio.cat, shape = Season),
             size = 3) +
  #geom_line(data = model.df, aes(x= mS, y = fit), colour = "black") +
  #geom_line(data = model.df, aes(x= mS, y = upr), colour = "black", linetype = "dashed") +
  #geom_line(data = model.df, aes(x= mS, y = lwr), colour = "black", linetype = "dashed") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  #labs(x = expression(paste(italic("m")["S"])), y = "Number of shared taxa between DNA-RNA") +
  scale_shape_manual(values = c(21,23,25)) +
  #scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

all.met$resid <- resid(lm(all.met$mBC-all.met$mS ~ 0))

model <- lm(all.met$mS ~ all.met$n.shared + all.met$rich.diff + all.met$replac)
summary(model)

model <- lm(all.met$mBC ~ all.met$n.shared + all.met$rich.diff + all.met$shar.diff)
summary(model)

rplc <- ggplot(all.met, aes(x = mS, y = replac)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_smooth(method = "lm", colour = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["S"])), y = "Replacement") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

rich <- ggplot(all.met, aes(x = mS, y = rich.diff)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_smooth(method = "lm", colour = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["S"])), y = "Richness difference") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

shar <- ggplot(all.met, aes(x = mS, y = n.shared)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_smooth(method = "lm", colour = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["S"])), y = "Number of shared OTUs") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

mS.plots <- ggarrange(shar, rich, rplc, ncol = 3, common.legend = T, legend = "right")

(outliers <- all.met %>%
    rstatix::identify_outliers(replac.diff) %>%
    select(DR.names, is.outlier, is.extreme))

# Exclude extreme outliers
plot.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

rplc <- ggplot(plot.df, aes(x = mBC, y = replac.diff)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_smooth(method = "lm", colour = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["BC"])), y = "Abundance difference of unshared OTUs") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

rich <- ggplot(all.met, aes(x = mBC, y = rich.diff)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_smooth(method = "lm", colour = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["BC"])), y = "Richness difference") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

shar <- ggplot(all.met, aes(x = mBC, y = n.shared)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_smooth(method = "lm", colour = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["BC"])), y = "Number of shared OTUs") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

(outliers <- all.met %>%
    rstatix::identify_outliers(shar.diff) %>%
    select(DR.names, is.outlier, is.extreme))

# Exclude extreme outliers
plot.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

shar.d <- ggplot(plot.df, aes(x = mBC, y = shar.diff)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  geom_smooth(method = "lm", colour = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  labs(x = expression(paste(italic("m")["BC"])), y = "Abundance difference of shared OTUs") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

ggarrange(shar, rich, rplc, shar.d, ncol = 4, common.legend = T, legend = "right")



ggplot(plot.df, aes(x = mBC, y = replac.diff)) +
  geom_point(aes(fill = sample.type.year, shape = Season),
             size = 3) +
  #geom_line(data = model.df, aes(x= mS, y = fit), colour = "black") +
  #geom_line(data = model.df, aes(x= mS, y = upr), colour = "black", linetype = "dashed") +
  #geom_line(data = model.df, aes(x= mS, y = lwr), colour = "black", linetype = "dashed") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.. , sep = "~~~")), 
               label.y = "top", label.x = "right",
               parse = TRUE) +
  #labs(x = expression(paste(italic("m")["S"])), y = "Number of shared taxa between DNA-RNA") +
  scale_shape_manual(values = c(21,23,25)) +
  scale_fill_manual(values = colvec, name = "Habitat Type") +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))

ggsave("./Figures/Final/lm_mS_tests.png", mS.plots, height =15, width=35, unit="cm")
ggsave("./Figures/Final/lm_mBC_tests.png", bc.plots, height =15, width=45, unit="cm")


# Kruskal-Wallis -------------------------------------------------------------------------------------
# Statistically test whether there is a difference between these categories
# 1. Mean abundance ----------------------------------------------------------------------------------
# Is there a difference in abundance difference between the categories identified by the mS:mBC ratio?

# Check ANOVA assumptions
# Do we have outliers?
(outliers <- all.met %>%
   group_by(ratio.cat) %>%
   rstatix::identify_outliers(mean.diff) %>%
   select(DR.names, ratio.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]
#aov.df <- all.met
# Check normality
model <- lm(mean.diff ~ ratio.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis instead of one-way ANOVA
(stat.kw <- kruskal.test(mean.diff ~ ratio.cat, data = aov.df)) # significantly different (*** level, p < 0.0001)

# Difference in mean abundance difference among >1 and <1 mS:mBC is statistically significant

(m.ab <- ggplot(aov.df, aes(x = ratio.cat, y = mean.abun.diff)) +
    geom_boxplot(width = 0.3) +
    labs(y = expression(paste("Mean | ", Delta, " OTU CSS reads |")), x = expression(paste(m[S]:m[BC]))) +
    annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(750,730,730), yend = c(750,750,750)) +
    annotate(geom = "text", x = 1.5, y = 800, label = paste("Kruskal-Wallis, p ", abbrev.p(stat.kw$p.value)[1])))

# 2. Mean abundance of dominant taxa ------------------------------------------------------------------
# Is there a difference in abundance of dominant OTUs between the categories identified by the mS:mBC ratio?
# Do we have outliers?
(outliers <- all.met %>%
   group_by(ratio.cat) %>%
   rstatix::identify_outliers(mean.abun.diff) %>%
   select(DR.names, ratio.cat, is.outlier, is.extreme))

# Exclude extreme outliers
aov.df <- all.met[!(DR.names %in% outliers[outliers$is.extreme == T,]$DR.names),]

# Check normality
model <- lm(mean.abun.diff ~ ratio.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(mean.abun.diff ~ ratio.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(dom.diff <- ggplot(aov.df, aes(x = ratio.cat, y = mean.abun.diff)) +
    geom_boxplot(width = 0.3) +
    labs(y = expression(atop(paste("Mean | ", Delta, " OTU CSS reads |"),"of abundant taxa")), x = expression(paste(m[S]:m[BC]))) +
    annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(600,580,580), yend = c(600,600,600)) +
    annotate(geom = "text", x = 1.5, y = 640, label = paste("Kruskal-Wallis, p ", abbrev.p(stat.kw$p.value)[1])))


# 3. Number of shared OTUs ------------------------------------------------------------------
# Is there a difference in number of shared taxa between the categories identified by the mS:mBC ratio?
# Do we have outliers?
(outliers <- all.met %>%
   group_by(ratio.cat) %>%
   rstatix::identify_outliers(n.shared) %>%
   select(DR.names, ratio.cat, is.outlier, is.extreme))

# None are extreme, do not exclude
aov.df <- all.met

# Check normality
model <- lm(n.shared ~ ratio.cat, data = aov.df)
ggqqplot(residuals(model)) # does not look good
rstatix::shapiro_test(residuals(model)) # significant = assumption not fulfilled

# Normality not fulfilled change to Kruskal Wallis
(stat.kw <- kruskal.test(n.shared ~ ratio.cat, data = aov.df)) # significantly different (** level, p < 0.01)

(n.sh <- ggplot(aov.df, aes(x = ratio.cat, y = n.shared)) +
    geom_boxplot(width = 0.3) +
    labs(y = "Mean number of shared OTUs", x = expression(paste(m[S]:m[BC]))) +
    annotate(geom = "segment", x = c(1,1,2), xend = c(2,1,2), y = c(770,750,750), yend = c(770,770,770)) +
    annotate(geom = "text", x = 1.5, y = 800, label = paste("Kruskal-Wallis, p ", abbrev.p(stat.kw$p.value)[1])))

(krusk.plots <- ggarrange(m.ab, dom.diff, n.sh, ncol = 3, align = "hv", labels = "auto"))

ggsave("./Figures/Final/ratio_krusk_tests.png", krusk.plots, width = 22, height = 10, units = "cm")

