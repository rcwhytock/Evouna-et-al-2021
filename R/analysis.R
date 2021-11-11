# Evouna et al 2021
# Created 1 November 2021
# Last modified 10 Novemeber 2021

# Load packages
library(vegan)
library(FSA)

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
ssMatrix <- t(ssMatrix[,-1])

# Store sample size (min count) for rarefaction
rareSample <- min(rowSums(ssMatrix))

# Calculate rarefied species richness for each site
rareResults <- rarefy(ssMatrix, sample = rareSample)
rareResults <- data.frame(site = names(rareResults),
                          rarefiedRichness = as.numeric(rareResults))

# Plot rarefaction curves # FIGURE 2 #
jpeg(
  filename = "../graphics/rarefiedRichnessCurves.jpg" ,
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
    "cornflowerblue",
    "cornflowerblue",
    "cornflowerblue",
    "chartreuse3",
    "chartreuse3",
    "chartreuse3",
    "darkgoldenrod1",
    "darkgoldenrod1",
    "darkgoldenrod1",
    "darkgoldenrod1",
    "brown4",
    "brown4",
    "brown4",
    "brown4",
    "cornflowerblue",
    "cornflowerblue",
    "cornflowerblue",
    "cornflowerblue",
    "cornflowerblue",
    "cornflowerblue"
  ),
  lty = c(1,
          1,
          1,
          2,
          2,
          2,
          3,
          3,
          3,
          3,
          4,
          4,
          4,
          4,
          1,
          1,
          1,
          1,
          1,
          1),
  lwd = 2,
  bty = "n",
  xlim = c(0, 200),
  ylim = c(0, 40),
  xlab = "Sample size"
)

legend(
  "bottomright",
  lty = c(1, 2, 3, 4),
  legend = c("Other", "LF", "LUS", "LS"),
  col = c("cornflowerblue", "chartreuse3", "darkgoldenrod1", "brown4"),
  bty = "n",
  lwd = 2
)
dev.off()

# Compare rarefied richness for burned savannas only # FIGURE 3b#
rareResults <-
  merge(rareResults, siteCodes) # merge final treatment names

burnedOnly <-
  rareResults[rareResults$TreatmentFinal %in% c("B", "LBS", "MN", "MS"),]

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
kruskal.test(x = burnedOnly$rarefiedRichness,
             g = burnedOnly$TreatmentFina)

dunnTest(rarefiedRichness ~ TreatmentFinal, burnedOnly)

jpeg(
  filename = "../graphics/burnedSavannaRarefied.jpg",
  width = 600,
  height = 600,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

boxplot(
  rarefiedRichness ~ numericTreatment,
  burnedOnly,
  ylab = "Rarefied species richness",
  xlab = "Area",
  axes = F,
  ylim = c(0, 15),
  cex.lab = 0.9
)

axis(1,
     at = c(1:4),
     labels = c("LBS", "B", "MN", "MS"))
axis(2, at = seq(0, 50, 5))
box()

mtext(
  side = 3,
  text = "(b)",
  adj = 0,
  line = 1
)
dev.off()

#### Abundance analysis ####

# Burned savanna abundance first
abundSav <-
  abundance[-which(abundance$Treatment %in% c("LF", "LUS")),]

sum(abundSav$Samples)
aggregate(abundSav$Samples ~ abundSav$Treatment, FUN = sum)

kruskal.test(x = abundSav$Samples,
             g = abundSav$Treatment)

dunnTest(Samples ~ Treatment, abundSav)

# Boxplots # FIGURE 3a #
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "LBS", 1, 0)
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "B", 2, abundSav$numericTreatment)
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "MN", 3, abundSav$numericTreatment)
abundSav$numericTreatment <-
  ifelse(abundSav$Treatment == "MS", 4, abundSav$numericTreatment)


jpeg(
  filename = "../graphics/numberOfEncounters.jpg",
  width = 600,
  height = 600,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

boxplot(
  Samples ~ numericTreatment,
  abundSav,
  ylab = "Number of samples",
  xlab = "Area",
  axes = F,
  ylim = c(0, 80),
  cex.lab = 0.9
)

axis(1,
     at = c(1:4),
     labels = c("LBS", "B", "MN", "MS"))
axis(2, at = seq(0, 80, 10))
box()

# Show significant differences LBS - MN
arrows(
  x0 = 1,
  y0 = 78,
  x1 = 3,
  y1 = 78,
  code = 3,
  angle = 90,
  length = 0.05
)
text(x = 2,
     y = 80,
     substitute(paste(italic('p'), " = 0.026")),
     cex = 0.5)


mtext(
  side = 3,
  text = "(a)",
  adj = 0,
  line = 1
)
dev.off()

#### Hierarchical clustering ####
ssMatrixSc <- scale(ssMatrix)
row.names(ssMatrixSc) <-
  siteCodes$Treatment # new treatment names in paper
dist_mat <- dist(ssMatrixSc, method = 'euclidean')

# Figure 4
hclust_avg <- hclust(dist_mat, method = 'ward.D2')
plot(hclust_avg)

#### Lope unburned savannas abundance ####
abundLop <-
  abundance[which(abundance$Treatment %in% c("LF", "LUS", "LBS")),]

sum(abundLop$Samples)
aggregate(abundLop$Samples ~ abundLop$Treatment, FUN = sum)

kruskal.test(x = abundLop$Samples,
             g = abundLop$Treatment)

dunnTest(Samples ~ Treatment, abundLop)

# Boxplots # FIGURE 5a #
abundLop$numericTreatment <-
  ifelse(abundLop$Treatment == "LF", 1, 0)
abundLop$numericTreatment <-
  ifelse(abundLop$Treatment == "LUS", 2, abundLop$numericTreatment)
abundLop$numericTreatment <-
  ifelse(abundLop$Treatment == "LBS", 3, abundLop$numericTreatment)

jpeg(
  filename = "../graphics/numberOfEncountersLope.jpg",
  width = 600,
  height = 600,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

boxplot(
  Samples ~ numericTreatment,
  abundLop,
  ylab = "Number of samples",
  xlab = "Area",
  axes = F,
  ylim = c(0, 250),
  cex.lab = 0.9
)

axis(1,
     at = c(1:3),
     labels = c("LF", "LUS", "LBS"))
axis(2, at = seq(0, 250, 25))
box()

mtext(
  side = 3,
  text = "(a)",
  adj = 0,
  line = 1
)
dev.off()

#### Lope rarefied richness ####
lopeOnly <-
  rareResults[which(rareResults$TreatmentFinal %in% c("LF", "LUS", "LBS")), ]

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


# Test of difference in rarefied species richness in Lope
kruskal.test(x = lopeOnly$rarefiedRichness,
             g = lopeOnly$TreatmentFina)

dunnTest(rarefiedRichness ~ TreatmentFinal, lopeOnly)

# Figure 5b
jpeg(
  filename = "../graphics/LopeRarefied.jpg",
  width = 600,
  height = 600,
  quality = 100,
  res = 150
)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))

boxplot(
  rarefiedRichness ~ numericTreatment,
  lopeOnly,
  ylab = "Rarefied species richness",
  xlab = "Area",
  axes = F,
  ylim = c(0, 15),
  cex.lab = 0.9
)

axis(1,
     at = c(1:3),
     labels = c("LF", "LBS", "LUS"))
axis(2, at = seq(0, 15, 5))
box()

mtext(
  side = 3,
  text = "(b)",
  adj = 0,
  line = 1
)
dev.off()

#### NMDS ####

# create a similarity matrix
row.names(ssMatrix)
lopeMatrix <- ssMatrix[4:14, ]

dis <- vegdist(lopeMatrix, method = "bray") # similarity matrix

# calculate the mds using the mds() function
n <- metaMDS(dis, k = 2) # mds

# Figure 6
jpeg(
  filename = paste0(fold, "/graphiques/NMDS_Lope.jpg"),
  width = 800,
  height = 800,
  quality = 100,
  res = 150
)

par(mar = c(4, 4, 1, 1))

plot(n,
     xlim = c(-.7, .7),
     ylim = c(-.7, .7),
     col = "white")
text(n$points, row.names(lopeMatrix), cex = 0.6)

dev.off()
