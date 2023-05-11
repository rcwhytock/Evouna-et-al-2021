# Evouna et al 2021
# Created 1 November 2021
# Last modified 21 April 2022

# Load packages
library(vegan)
library(FSA)
library(dendextend)
library(cluster)
library(emmeans)
library(MASS)

#### Load the species x sample matrix for rarefaction ####
ssMatrix <- read.csv("../data/species_x_sample_matrix.csv")
head(ssMatrix)

length(unique(ssMatrix$bar.code.species)) # 105 species

# Load 'abundance' data with number of samples collected on each transect
abundance <- read.csv("../data/abundance_data.csv")
head(abundance)

# Read in site codes to get correct treatment names
siteCodes <- read.csv("../data/siteCodes.csv")
siteCodes

sum(abundance$Samples) # 1336 samples

#### Rarefaction analysis ####

# Transpose for rarefaction
ssMatrix <- t(ssMatrix[, -1])

# Store sample size (min count) for rarefaction
rareSample <- min(rowSums(ssMatrix))

# Calculate rarefied species richness for each site
rareResults <- rarefy(ssMatrix, sample = rareSample)
rareResults <- data.frame(site = names(rareResults),
                          rarefiedRichness = as.numeric(rareResults))

# Plot rarefaction curves # FIGURE 2 #
jpeg(
  filename = "../graphics/Fig2.jpg" ,
  width = 800,
  height = 800,
  quality = 100,
  res = 150
)

rarecurve(
  ssMatrix,
  step = 1,
  label = F,
  col = c(
    "grey50", #other
    "grey50", #other
    "grey50", #other
    "forestgreen", #LF
    "forestgreen", #LF
    "forestgreen", #LF
    "darkorange", #LUS
    "darkorange", #LUS
    "darkorange", #LUS
    "darkorange", #LUS
    "blue", #LBS
    "blue", #LBS
    "blue", #LBS
    "blue", #LBS
    "grey50", #other
    "grey50", #other
    "grey50", #other
    "grey50", #other
    "grey50", #other
    "grey50" #other
  ),
  lty = c(2,
          2,
          2,
          3, 
          3, 
          3, 
          1,#
          1,#
          1,#
          1,#
          4,
          4,
          4,
          4,
          2,
          2,
          2,
          2,
          2,
          2),
  lwd = 2,
  bty = "n",
  xlim = c(0, 200),
  ylim = c(0, 40),
  xlab = "Sample size"
)

legend(
  "bottomright",
  lty = c(2, 3, 1, 4),
  legend = c("Other", "LF", "LUS", "LBS"),
  col = c("grey28", "forestgreen", "darkorange", "blue"),
  bty = "n",
  lwd = 2
)
dev.off()

# Compare rarefied richness for burned savannas only # FIGURE 3b#
rareResults <-
  merge(rareResults, siteCodes) # merge final treatment names

burnedOnly <-
  rareResults[rareResults$TreatmentFinal %in% c("B", "LBS", "MN", "MS"), ]

# Create numeric variable to plot in correct order for paper
burnedOnly$numericTreatment <-
  ifelse(burnedOnly$TreatmentFinal == "LBS", 1, 0)
burnedOnly$numericTreatment <-
  ifelse(burnedOnly$TreatmentFinal == "B",
         2,
         burnedOnly$numericTreatment)
burnedOnly$numericTreatment <-
  ifelse(burnedOnly$TreatmentFinal == "MN",
         3,
         burnedOnly$numericTreatment)
burnedOnly$numericTreatment <-
  ifelse(burnedOnly$TreatmentFinal == "MS",
         4,
         burnedOnly$numericTreatment)

# Test of difference in rarefied species richness in burned savs
mod1 <- glm(rarefiedRichness ~ TreatmentFinal-1, data = burnedOnly) # test difference from 0 for each treatment
hist(resid(mod1))
par(mfrow = c(2,2))
plot(mod1) # ok
summary(mod1)

confints1 <- data.frame(emmeans(mod1, ~TreatmentFinal)) # note that df Inf is ok, see https://cran.r-project.org/web/packages/emmeans/vignettes/FAQs.html#asymp

confints1$TreatmentNumeric <-
  ifelse(confints1$TreatmentFinal == "LBS", 1, 0)
confints1$TreatmentNumeric <-
  ifelse(confints1$TreatmentFinal == "B",
         2,
         confints1$TreatmentNumeric)
confints1$TreatmentNumeric <-
  ifelse(confints1$TreatmentFinal == "MN",
         3,
         confints1$TreatmentNumeric)
confints1$TreatmentNumeric <-
  ifelse(confints1$TreatmentFinal == "MS",
         4,
         confints1$TreatmentNumeric)

# Post hoc pairwase comparison, tukey method for family of 4 estimates
pairsRareBurned <- pairs(emmeans(mod1, ~TreatmentFinal)) 
pairsRareBurned <- data.frame(pairsRareBurned )
write.csv(pairsRareBurned , "../results/TableS3.csv", row.names = F) 

    # Create transparent color for plotting raw data points
mycol <- rgb(0, 0, 0, alpha = 0.35, names = "black")

# Add numeric treatment name for plotting raw data
rareResults$TreatmentNumeric <-
  ifelse(rareResults$TreatmentFinal == "LBS", 1, 0)
rareResults$TreatmentNumeric <-
  ifelse(rareResults$TreatmentFinal == "B",
         2,
         rareResults$TreatmentNumeric)
rareResults$TreatmentNumeric <-
  ifelse(rareResults$TreatmentFinal == "MN",
         3,
         rareResults$TreatmentNumeric)
rareResults$TreatmentNumeric <-
  ifelse(rareResults$TreatmentFinal == "MS",
         4,
         rareResults$TreatmentNumeric)

# Plot
jpeg(
  filename = "../graphics/Fig3b.jpg",
  width = 600,
  height = 550,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

plot(emmean ~ TreatmentNumeric, confints1, axes = F, xlim = c(0.8,4.5), pch = 16,
     ylim = c(0,15), ylab = "Rarefied species richness", xlab = "Sampling area")
points(rarefiedRichness ~ jitter(TreatmentNumeric, factor = 0.2), rareResults, col = mycol, pch = 16)

arrows(x0 = confints1$TreatmentNumeric, 
       x1 = c(2,1,3,4),
       y0 = confints1$lower.CL,
       y1 = confints1$upper.CL,
       code = 3, angle = 90,
       length = 0.05)

axis(1,
     at = c(1:4),
     labels = c("LBS", "B", "MN", "MS"))
axis(2, at = seq(0, 15, 5))
#box()

mtext(
  side = 3,
  text = "(b)",
  adj = 0,
  line = 1
)

# Significant pairwise differences
arrows(
  x0 = 1,
  y0 = 13,
  x1 = 3,
  y1 = 13,
  code = 3,
  angle = 90,
  length = 0.03
)

text(x = 2,
     y = 14,
     substitute(paste(italic('p'), " = 0.039")),
     cex = 0.5)

dev.off()

#### Abundance analysis ####

# Burned savanna abundance first
abundSav <-
  abundance[-which(abundance$Treatment %in% c("LF", "LUS")), ]

sum(abundSav$Samples)
aggregate(abundSav$Samples ~ abundSav$Treatment, FUN = sum)

# Test of difference in rarefied species richness in burned savs
mod2 <- glm.nb(Samples ~ Treatment-1, data = abundSav) # test difference from 0 for each treatment
hist(resid(mod2))
par(mfrow = c(2,2))
plot(mod2) # ok
summary(mod2)

confints2 <- data.frame(emmeans(mod2, ~Treatment)) # note that df Inf is ok, see https://cran.r-project.org/web/packages/emmeans/vignettes/FAQs.html#asymp

confints2$TreatmentNumeric <-
  ifelse(confints2$Treatment == "LBS", 1, 0)
confints2$TreatmentNumeric <-
  ifelse(confints2$Treatment == "B",
         2,
         confints2$TreatmentNumeric)
confints2$TreatmentNumeric <-
  ifelse(confints2$Treatment == "MN",
         3,
         confints2$TreatmentNumeric)
confints2$TreatmentNumeric <-
  ifelse(confints2$Treatment == "MS",
         4,
         confints2$TreatmentNumeric)

# POst hoc pairwase comparison, tukey method for family of 4 estimates
pairsAbundBurned <- pairs(emmeans(mod2, ~Treatment)) 
pairsAbundBurned <- data.frame(pairsAbundBurned)
write.csv(pairsAbundBurned, "../results/TableS2.csv", row.names = F) 


# Plots # FIGURE 3a #
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "LBS", 1, 0)
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "B", 2, abundSav$numericTreatment)
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "MN", 3, abundSav$numericTreatment)
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "MS", 4, abundSav$numericTreatment)

# Plot
jpeg(
  filename = "../graphics/Fig3a.jpg",
  width = 600,
  height = 550,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

plot(exp(emmean) ~ TreatmentNumeric, confints2, axes = F, xlim = c(0.8,4.5), pch = 16,
     ylim = c(0,110), ylab = "Number of samples", xlab = "Sampling area")
points(Samples ~ jitter(numericTreatment, factor = 0.2), abundSav, col = mycol, pch = 16)

arrows(x0 = confints2$TreatmentNumeric, 
       x1 = confints2$TreatmentNumeric,
       y0 = exp(confints2$asymp.LCL),
       y1 = exp(confints2$asymp.UCL),
       code = 3, angle = 90,
       length = 0.05)

axis(1,
     at = c(1:4),
     labels = c("LBS", "B", "MN", "MS"))
axis(2, at = seq(0, 120, 20))
#box()

mtext(
  side = 3,
  text = "(a)",
  adj = 0,
  line = 1
)

# Significant pairwise differences
arrows(
  x0 = 1,
  y0 = 96,
  x1 = 3,
  y1 = 96,
  code = 3,
  angle = 90,
  length = 0.03
)

text(x = 2,
     y = 100,
     substitute(paste(italic('p'), " = 0.001")),
     cex = 0.5)

arrows(
  x0 = 1,
  y0 = 106,
  x1 = 4,
  y1 = 106,
  code = 3,
  angle = 90,
  length = 0.03
)

text(x = 2.5,
     y = 110,
     substitute(paste(italic('p'), " = 0.029")),
     cex = 0.5)

dev.off()


#### Hierarchical clustering ####
row.names(ssMatrix) <-
  siteCodes$Treatment # new treatment names in paper

row.names(ssMatrix)[7] <- "LUS13.1"
row.names(ssMatrix)[8] <- "LUS20.1"
row.names(ssMatrix)[9] <- "LUS13.2"
row.names(ssMatrix)[10] <- "LUS20.2"


dist_mat <- vegdist(ssMatrix, method = 'bray')
hc <- agnes(dist_mat)
dend <- as.dendrogram(hc)

# Figure 4
jpeg(
  filename = "../graphics/Fig4.jpg",
  width = 800,
  height = 800,
  quality = 100,
  res = 150
)


par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
dend <- dend %>%
  color_branches(k = 3) %>%
  set("branches_lwd", c(2, 1, 2)) %>%
  set("branches_lty", c(1, 2, 1))

dend <- color_labels(dend, k = 3)

plot(dend, horiz = T)

dev.off()

#### Lope habitat abundance ####
abundLop <-
  abundance[which(abundance$Treatment %in% c("LF", "LUS", "LBS")), ]

sum(abundLop$Samples)
aggregate(abundLop$Samples ~ abundLop$Treatment, FUN = sum)

# Test of difference in rarefied species richness in burned savs
mod3 <- glm.nb(Samples ~ Treatment-1, data = abundLop) # test difference from 0 for each treatment
hist(resid(mod3))
par(mfrow = c(2,2))
plot(mod3) # ok
summary(mod3)

confints3 <- data.frame(emmeans(mod3, ~Treatment)) # note that df Inf is ok, see https://cran.r-project.org/web/packages/emmeans/vignettes/FAQs.html#asymp

confints3$TreatmentNumeric <-
  ifelse(confints3$Treatment == "LF", 1, 0)
confints3$TreatmentNumeric <-
  ifelse(confints3$Treatment == "LUS",
         2,
         confints3$TreatmentNumeric)
confints3$TreatmentNumeric <-
  ifelse(confints3$Treatment == "LBS",
         3,
         confints3$TreatmentNumeric)

# Post hoc pairwase comparison, tukey method
pairsLopeAbund <- pairs(emmeans(mod3, ~Treatment)) 
pairsLopeAbund <- data.frame(pairsLopeAbund)
write.csv(pairsLopeAbund, "../results/TableS4.csv", row.names = F) 

# Plots
abundLop$numericTreatment <-
  ifelse(abundLop$Treatment == "LF", 1, 0)
abundLop$numericTreatment <-
  ifelse(abundLop$Treatment == "LUS", 2, abundLop$numericTreatment)
abundLop$numericTreatment <-
  ifelse(abundLop$Treatment == "LBS", 3, abundLop$numericTreatment)

# Plot
jpeg(
  filename = "../graphics/Fig5a.jpg",
  width = 600,
  height = 550,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

plot(exp(emmean) ~ TreatmentNumeric, confints3, axes = F, xlim = c(0.8,3.2), pch = 16,
     ylim = c(0,300), ylab = "Number of samples", xlab = "Sampling area")
points(Samples ~ jitter(numericTreatment, factor = 0.2), abundLop, col = mycol, pch = 16)

arrows(x0 = confints3$TreatmentNumeric, 
       x1 = confints3$TreatmentNumeric,
       y0 = exp(confints3$asymp.LCL),
       y1 = exp(confints3$asymp.UCL),
       code = 3, angle = 90,
       length = 0.05)

axis(1,
     at = c(1:3),
     labels = c("LF", "LUS", "LBS"))
axis(2, at = seq(0, 300, 50))
#box()

mtext(
  side = 3,
  text = "(a)",
  adj = 0,
  line = 1
)

# Significant pairwise differences
arrows(
  x0 = 1,
  y0 = 260,
  x1 = 3,
  y1 = 260,
  code = 3,
  angle = 90,
  length = 0.03
)

text(x = 2,
     y = 270,
     substitute(paste(italic('p'), " < 0.001")),
     cex = 0.5)

arrows(
  x0 = 1,
  y0 = 290,
  x1 = 2,
  y1 = 290,
  code = 3,
  angle = 90,
  length = 0.03
)

text(x = 1.5,
     y = 300,
     substitute(paste(italic('p'), " < 0.001")),
     cex = 0.5)

dev.off()

#### Lope rarefied richness ####
lopeOnly <-
  rareResults[which(rareResults$TreatmentFinal %in% c("LF", "LUS", "LBS")),]

lopeOnly <- lopeOnly[,c(1,2,3,4)]

# Create numeric variable to plot in correct order for paper
lopeOnly$numericTreatment <-
  ifelse(lopeOnly$TreatmentFinal == "LF", 1, 0)
lopeOnly$numericTreatment <-
  ifelse(lopeOnly$TreatmentFinal == "LUS",
         2,
         lopeOnly$numericTreatment)
lopeOnly$numericTreatment <-
  ifelse(lopeOnly$TreatmentFinal == "LBS",
         3,
         lopeOnly$numericTreatment)


# Test of difference in rarefied species richness in burned savs
mod4 <- glm(rarefiedRichness ~ TreatmentFinal-1, data = lopeOnly) # test difference from 0 for each treatment
hist(resid(mod4))
par(mfrow = c(2,2))
plot(mod4) # ok
summary(mod4)

confints4 <- data.frame(emmeans(mod4, ~TreatmentFinal)) # note that df Inf is ok, see https://cran.r-project.org/web/packages/emmeans/vignettes/FAQs.html#asymp

confints4$TreatmentNumeric <-
  ifelse(confints4$TreatmentFinal == "LF", 1, 0)
confints4$TreatmentNumeric <-
  ifelse(confints4$TreatmentFinal == "LUS",
         2,
         confints4$TreatmentNumeric)
confints4$TreatmentNumeric <-
  ifelse(confints4$TreatmentFinal == "LBS",
         3,
         confints4$TreatmentNumeric)

# POst hoc pairwase comparison, tukey method for family of 4 estimates
pairsLopeRare <- pairs(emmeans(mod4, ~TreatmentFinal)) 
pairsLopeRare <- data.frame(pairsLopeRare)
write.csv(pairsLopeRare, "../results/TableS5.csv", row.names = F) 

# Plot
jpeg(
  filename = "../graphics/Fig5b.jpg",
  width = 600,
  height = 550,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

plot(emmean ~ TreatmentNumeric, confints4, axes = F, xlim = c(0.8,3.2), pch = 16,
     ylim = c(0,20), ylab = "Rarefied species richness", xlab = "Sampling area")
points(rarefiedRichness ~ jitter(numericTreatment, factor = 0.2), lopeOnly, col = mycol, pch = 16)

arrows(x0 = confints4$TreatmentNumeric, 
       x1 = confints4$TreatmentNumeric,
       y0 = confints4$lower.CL,
       y1 = confints4$upper.CL,
       code = 3, angle = 90,
       length = 0.05)

axis(1,
     at = c(1:3),
     labels = c("LF", "LUS+LIBS", "LBS"))
axis(2, at = seq(0, 20, 5))
#box()

mtext(
  side = 3,
  text = "(b)",
  adj = 0,
  line = 1
)

# Significant pairwise differences
arrows(
  x0 = 1,
  y0 = 13,
  x1 = 3,
  y1 = 13,
  code = 3,
  angle = 90,
  length = 0.03
)

text(x = 2,
     y = 14,
     substitute(paste(italic('p'), " < 0.001")),
     cex = 0.5)

arrows(
  x0 = 2,
  y0 = 15,
  x1 = 3,
  y1 = 15,
  code = 3,
  angle = 90,
  length = 0.03
)

text(x = 2.5,
     y = 16,
     substitute(paste(italic('p'), " < 0.001")),
     cex = 0.5)

dev.off()

#### NMDS ####

# create a similarity matrix
row.names(ssMatrix) <-
  siteCodes$Treatment
lopeMatrix <- ssMatrix[4:14,]

dis <- vegdist(lopeMatrix, method = "bray") # similarity matrix

# calculate the nmds
n <- metaMDS(dis, k = 2) # mds

# Figure 6
jpeg(
  filename = "../graphics/Fig6.jpg",
  width = 800,
  height = 800,
  quality = 100,
  res = 150
)

par(mar = c(4, 4, 1, 1))

plot(n,
     xlim = c(-.6, .45),
     ylim = c(-.51, .4),
     type = "n", bty = "n")

points(
  jitter(n$points, factor = 100),
  pch = c(22,22,22,21,24,21,24,23,23,23,23),
  bg = c(
    "forestgreen",
    "forestgreen",
    "forestgreen",
    "darkorange",
    "darkorange",#
    "darkorange",
    "darkorange",#
    "blue",
    "blue",
    "blue",
    "blue"
  ),
  cex = 1.8
)

# text(n$points, row.names(lopeMatrix),
#      col = c(
#   "chartreuse3",
#   "chartreuse3",
#   "chartreuse3",
#   "darkgoldenrod1",
#   "darkgoldenrod1",
#   "darkgoldenrod1",
#   "darkgoldenrod1",
#   "brown4",
#   "brown4",
#   "brown4",
#   "brown4"
# ),cex = 0.8)

legend(
  "bottomright",
  pch = c(22,24,21,23),
  lty = NULL,
  legend = c("LF", "LUS 20", "LUS 13", "LBS"),
  pt.bg = c("forestgreen", "darkorange", "darkorange", "blue"),
  col = c("darkgreen",  "darkorange", "darkorange", "blue"),
  bty = "n",
  cex = 1
)

dev.off()

