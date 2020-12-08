# Hypothetical communitites

library(vegan)

# 100% similar
com1 <- data.frame(OTU1 = c(3,3), OTU2 = c(12, 12), OTU3 = c(0,0), OTU4 = c(30,30), OTU5 = c(3, 3), row.names = c("DNA","RNA"))

bray <- vegdist(com1, method = "bray")
sor <- vegdist(com1, method = "bray", binary = T)

c(bray, sor)

# 100% dissimilar

com2 <- data.frame(OTU1 = c(3,0), OTU2 = c(0, 12), OTU3 = c(0,0), OTU4 = c(30,0), OTU5 = c(3, 0), row.names = c("DNA","RNA"))

bray <- vegdist(com2, method = "bray")
sor <- vegdist(com2, method = "bray", binary = T)

c(bray, sor)

# BC = S
# only happens when the ratio between
# BC = sum of lesser abundance of shared taxa : sum of abundance of all taxa
# Sor = number of shared taxa : sum of the number of observations
# is the same
com3 <- data.frame(OTU1 = c(3,0), OTU2 = c(12, 12), OTU3 = c(0,3), OTU4 = c(30,0), OTU5 = c(0, 0), row.names = c("DNA","RNA"))

bray <- vegdist(com3, method = "bray")
sor <- vegdist(com3, method = "bray", binary = T)

c(bray, sor)

# bray curtis always favours abundant taxa
# BC < S
# a few taxa with big abundance difference in DNA and RNA
# abundances of shared taxa are very similar, the higher the abundance match, the more similar in BC
com4 <- data.frame(OTU1 = c(3,0), OTU2 = c(12, 0), OTU3 = c(0,3), OTU4 = c(120,120), OTU5 = c(0, 0), row.names = c("DNA","RNA"))

bray <- vegdist(com4, method = "bray")
sor <- vegdist(com4, method = "bray", binary = T)

c(bray, sor)

com5 <- data.frame(OTU1 = c(3,0), OTU2 = c(12, 0), OTU3 = c(0,3), OTU4 = c(30,30), OTU5 = c(0, 0), row.names = c("DNA","RNA"))

bray <- vegdist(com5, method = "bray")
sor <- vegdist(com5, method = "bray", binary = T)

c(bray, sor)

# BC > S
com6 <- data.frame(OTU1 = c(3,0), OTU2 = c(12, 12), OTU3 = c(10,3), OTU4 = c(10,0), OTU5 = c(0, 0), row.names = c("DNA","RNA"))

bray <- vegdist(com6, method = "bray")
sor <- vegdist(com6, method = "bray", binary = T)

c(bray, sor)

com7 <- data.frame(OTU1 = c(3,0), OTU2 = c(12, 12), OTU3 = c(0,3), OTU4 = c(60,0), OTU5 = c(0, 0), row.names = c("DNA","RNA"))

bray <- vegdist(com7, method = "bray")
sor <- vegdist(com7, method = "bray", binary = T)

c(bray, sor)



# BC > S
com6 <- data.frame(OTU1 = c(3,0), OTU2 = c(20, 5), OTU3 = c(0,3), OTU4 = c(200,0), OTU5 = c(0, 0), row.names = c("DNA","RNA"))

bray <- vegdist(com6, method = "bray")
sor <- vegdist(com6, method = "bray", binary = T)

c(bray, sor)

com7 <- data.frame(OTU1 = c(4,0), OTU2 = c(12, 2), OTU3 = c(2,0), OTU4 = c(120,60), OTU5 = c(0, 1), row.names = c("DNA","RNA"))

bray <- vegdist(com7, method = "bray")
sor <- vegdist(com7, method = "bray", binary = T)

c(bray, sor)



ex <- data.frame(g = c(6,10), gup = c(7,0), rain = c(4,6))
bray <- vegdist(ex, method = "bray")
sor <- vegdist(ex, method = "bray", binary = T)

c(bray, sor)
