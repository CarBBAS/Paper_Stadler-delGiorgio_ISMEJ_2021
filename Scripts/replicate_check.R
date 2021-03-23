# Replicate check ----------------------------------------------------------------------------
# get all samples with a technical (GQ) replicate

dupl.samp <- meta[seq_name %in% row.names(seqtab)[grep("-s2-", row.names(seqtab))],]$dr_match_name

sub.meta <- meta[dr_match_name %in% dupl.samp,]
sub.seq <- seqtab[row.names(seqtab) %in% sub.meta$seq_name,]

sub.meta <- sample_data(sub.meta)

# Assign rownames to be Sample ID's
rownames(sub.meta) <- sub.meta$seq_name
sub.meta <- sub.meta[order(rownames(sub.meta)),]

# Construct phyloseq object
sub.ps <- phyloseq(otu_table(sub.seq, taxa_are_rows = F),
                   sample_data(sub.meta),
                   tax_table(tax))

# Filter only bacteria, omitting chloroplasts and mitochondria
sub.ps <- sub.ps %>%
  subset_taxa(domain == "Bacteria" &
                family  != "Mitochondria" &
                order   != "Chloroplast")


# remove ASVs that do not appear in this dataset
sub.ps <- prune_taxa(taxa_sums(sub.ps) != 0, sub.ps)
sub.ps <- prune_samples(sample_sums(sub.ps) != 0, sub.ps) # remove samples with no reads

# extract ASV matrix
pb.mat <- otu_mat(sub.ps) # convert phyloseq obj to matrix
#pb.mat <- decostand(pb.mat, "hellinger")

pb.bray <- vegan::vegdist(pb.mat, method = "bray")
dist.mat <- melt.dist(pb.bray)
setDT(dist.mat)
dist.mat[meta, "x_match" := i.dr_match_name , on = c("Sample.x" = "seq_name")]
dist.mat[meta, "y_match" := i.dr_match_name , on = c("Sample.y" = "seq_name")]

# keep only those rows where we are comparing replicates
dist.mat[x_match == y_match & dist > 0.8,] # those that are very different are those that compare shallow vs deep
# focus on the actual comparisons
dist.mat[x_match == y_match & dist < 0.8,] # usually below < 0.4, but still dissimilarity is quite high for replicates

# So which of the replicates should we keep?
# Calculate the number of reads
otus <- otu_mat(sub.ps)
shallow <- otus[-grep("DSeq", row.names(otus)),]
deep <- otus[grep("DSeq", row.names(otus)),]
s <- data.frame(apply(shallow, 1, sum)) # s2 has always more reads
d <- data.frame(apply(deep, 1, sum)) # s2 has always more reads

summary(s[seq(1, nrow(s), 2),]) # s2
summary(s[seq(2, nrow(s), 2),]) # no s2
# no s2 more reads

summary(d[seq(1, nrow(d), 2),]) # s2
summary(d[seq(2, nrow(d), 2),]) # no s2
# no s2 more reads

# What about all the other samples?
# phyloseq needs the sample names of the meta data to be the same as the microbial data

ex.meta <- meta[!(dr_match_name %in% dupl.samp),]
ex.seq <- seqtab[row.names(seqtab) %in% ex.meta$seq_name,]
ex.meta <- sample_data(ex.meta)

# Assign rownames to be Sample ID's
rownames(ex.meta) <- ex.meta$seq_name
ex.meta <- ex.meta[order(rownames(ex.meta)),]

# Construct phyloseq object
ps <- phyloseq(otu_table(ex.seq, taxa_are_rows = F),
               sample_data(ex.meta),
               tax_table(tax))

# Filter only bacteria, omitting chloroplasts and mitochondria
pb <- ps %>%
  subset_taxa(domain == "Bacteria" &
                family  != "Mitochondria" &
                order   != "Chloroplast")

# remove ASVs that do not appear in this dataset
pb <- prune_taxa(taxa_sums(pb) != 0, pb)
pb <- prune_samples(sample_sums(pb) != 0, pb) # remove samples with no reads

# extract ASV matrix
otus <- otu_mat(pb)
shallow <- otus[-grep("DSeq", row.names(otus)),]
deep <- otus[grep("DSeq", row.names(otus)),]
as <- data.frame(apply(shallow, 1, sum)) # s2 has always more reads
ad <- data.frame(apply(deep, 1, sum)) # s2 has always more reads

par(mfrow = c(1,3))
hist(as[,1]); hist(s[seq(1, nrow(s), 2),]); hist(s[seq(2, nrow(s), 2),])

hist(ad[,1])
summary(as); summary(d[seq(1, nrow(d), 2),])
# no s2 more reads

summary(ad)
# no s2 more reads


# Continue with analysis ---------------------------------------------------------------------

# phyloseq needs the sample names of the meta data to be the same as the microbial data
meta <- sample_data(meta)

# Assign rownames to be Sample ID's
rownames(meta) <- meta$seq_name
meta <- meta[order(rownames(meta)),]

# Construct phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = F),
               sample_data(meta),
               tax_table(tax))

# How many unclassified OTUs?
t <- ps %>% subset_taxa(taxa_sums(ps) != 0) %>% subset_taxa(is.na(domain))
nrow(tax_mat(t)) # 50602
# of a total
nrow(tax_mat(ps)) # 93710
rm(t)

# Filter only bacteria, omitting chloroplasts and mitochondria
pb <- ps %>%
  subset_taxa(domain == "Bacteria" &
                family  != "Mitochondria" &
                order   != "Chloroplast")






# do a quick PCoA to check difference between duplicates (NMDS does not converge)
pb.mat <- otu_mat(all) # convert phyloseq obj to matrix
pb.mat <- decostand(pb.mat, "hellinger")

# PCoA with Bray-Curtis
pb.bray <- vegdist(pb.mat, method = "bray")
pb.bray <- sqrt(pb.bray) # make Euclidean

ncol(pb.mat) # OTUs
nrow(pb.mat) # Samples

# make PCoA
pb.bray.pcoa <- ape::pcoa(pb.bray) # 372 registers
# plot with custom function (= made to avoid repetitive code)
# custom function is in ./Functions/custom_fun.R
pcoa <- pb.bray.pcoa; physeq <- all; plot.axes <- c(1,2); axes <- plot.axes

# extract scores and variance explained
pdataframe <- data.frame(Sample = as.character(row.names(pcoa$vectors)),
                         pcoa$vectors[,axes[1:length(axes)]],
                         stringsAsFactors = F)  # extract site scores

pb.var <- data.frame(Axes = axes,
                     var = round(100 * pcoa$values$Eigenvalues[axes] / sum(pcoa$values$Eigenvalues), 2),
                     stringsAsFactors = F)

# merge with a selection of meta data
meta <- data.frame(Sample = as.character(row.names(sample_df(physeq))),
                   sample_df(physeq) %>% dplyr::select(sample.type.year, Season, year, 
                                                       dna_type, distance.from.mouth, dr_match_name, replicate), 
                   stringsAsFactors = F)
pdataframe <- merge(pdataframe, meta, by = "Sample")
pdataframe$Sample <- as.character(pdataframe$Sample)


# merge some sample types
pdataframe$sample.type.year <- factor(pdataframe$sample.type.year, levels = c("Soil","Sediment",
                                                                              "Soilwater","Hyporheicwater", 
                                                                              "Wellwater","Stream", "Tributary",
                                                                              "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                              "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                                              "Marine"),
                                      labels = c("Soil","Sediment",
                                                 "Soilwater","Soilwater", 
                                                 "Groundwater","Stream", "Tributary",
                                                 "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
                                                 "Upriver",
                                                 "Reservoirs","Reservoirs", "Reservoirs","Reservoirs", "Downriver",
                                                 "Estuary"))

pdataframe$Season <- factor(pdataframe$Season, levels = c("Spring","Summer","Autumn"),
                            labels = c("Spring","Summer","Autumn"))

pdataframe$DnaType <- factor(pdataframe$dna_type, levels = c("DNA","cDNA"),
                             labels = c("DNA","RNA"))

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
            "#375A8CFF", #Downriver,
            "gray40",
            "red") #Estuary)
# export factors for colouring
sample.factors <- levels(pdataframe$sample.type.year)
names(colvec) <- as.character(sample.factors) # assign names for later easy selection

# main plot
p <- ggplot(pdataframe, aes_string(x = paste0("Axis.", plot.axes[1]), 
                                   y = paste0("Axis.", plot.axes[2]))) +
  theme_cust() +
  geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
  geom_point(aes(fill = sample.type.year, shape = Season, colour = dna_type, size = as.character(replicate))) +
  geom_line(data = pdataframe[pdataframe$dr_match_name %in% pdataframe[pdataframe$replicate == 2,]$dr_match_name &
                                pdataframe$dna_type == "DNA",],
            aes(group = dr_match_name), alpha = 0.5, colour = "red") +
  coord_fixed(1) + # ensure aspect ratio
  scale_shape_manual(values = c(21,23,25)) +
  scale_colour_manual(values = c("black", "white"), name = "Nucleic Acid \nType") +
  scale_size_manual(values = c(2.5, 3.5), name = "Replicate") +
  scale_fill_manual(values = colvec) +
  labs(x = paste0("PC",  plot.axes[1]," [ ", pb.var[pb.var$Axes == plot.axes[1],"var"]," %]"), 
       y = paste0("PC",  plot.axes[2]," [ ", pb.var[pb.var$Axes == plot.axes[2],"var"]," %]")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(shape = guide_legend(order = 2, override.aes=list(size = 2)),
         fill = guide_legend(order = 1, override.aes=list(shape=21, size = 2)))










# Phantom taxa check

# How many taxa only appear in RNA per sample?
temp <- sumdf[OTU == "OTU_1" & dna_type == "cDNA", ]
t<-temp[!duplicated(dr_match_name),]$dr_match_name

temp <- sumdf[dr_match_name %in% t,]
temp <- temp[, .(mean = mean(reads, na.rm = T)), by = .(OTU, dna_type, dr_match_name)]

rna.obs <- dcast(temp, dr_match_name + OTU ~ dna_type, value.var = "mean")
rna.obs[DNA > 0 | cDNA > 0, n.all := .N, .(dr_match_name)]
tt <- rna.obs[cDNA >= 1 & DNA == 0,]

rnaonly <- tt[, .(n = .N,
                  n.all = unique(n.all)), .(dr_match_name)]
rnaonly <- rnaonly[,prop := n * 100 / n.all]
# quite some...

comp.seq <- function (str1, str2) {
  tmp = cbind (strsplit(as.character(str1), ""), strsplit(as.character(str2), ""))
  t(apply (tmp, 1, function(x){x[[1]]==x[[2]]}))
}
# nucleotide difference... different OTUs
otu1 <- "TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTATACAAGACAGGCGTGAAATCCCCGGGCTTAACCTGGGAATGGCGCCTGTGACTGTATAGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGATATGTGGAGGAATACCAATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG"
otu2 <- "TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTATACAAGACAGGCGTGAAATCCCCGGGCTTAACCTGGGAATGGCGCCTGTGACTGTATAGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGATATGTGGAGGAATACCAATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG"
otu3 <- "GACGAACCGTGCAAACGTTATTCGGAATCACTGGGCTTAAAGGGCGCGTAGGCGGGTGATCAAGTCAATGGTGAAATCCTCCAGCTTAACTGGAGAAGTGCCTTTGATACTGATTGTCTAGAGGGAGGTAGGGGCATGTGGAACTTCAGGTGGAGCGGTGAAATGCGTAGATATCTGAAGGAACGCCAGTGGCGAAAGCGATGTGCTGGACCTCTTCTGACGCTGAGGCGCGAAAGCTAGGGGATCAAACGGG"
# just one nucleotide difference
which(comp.seq(otu1, otu2) == F)
which(comp.seq(otu2, otu3) == F)
tax.df <- setDT(data.frame(tax.df), keep.rownames = "OTU")
tax.df[!is.na(genus), n := .N, .(genus)]
tax.df[is.na(genus) & !is.na(family), n := .N, .(family)]
only.in.rna <- unique(tt$OTU)
View(tax.df[OTU %in% only.in.rna,])