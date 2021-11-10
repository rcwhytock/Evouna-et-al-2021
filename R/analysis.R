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

sum(abundance$Samples) # 1336 samples

#### Rarefaction analysis ####

# Transpose for rarefaction
ssMatrix <- t(ssMatrix[,-1])

# Calculate abundance per study site

# Load new plotting function for rarecurve
source("../R/functions/rarec.R")

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
  lty = c(
    1,
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
    1
  ),
  lwd = 2,
  bty = "n",
  xlim = c(0,200),
  ylim = c(0,40),
  xlab = "Sample size")

legend(
  "bottomright",
  lty = c(1, 2, 3, 4),
  legend = c("Other", "LF", "LUS", "LS"),
  col = c("cornflowerblue", "chartreuse3", "darkgoldenrod1", "brown4"),
  bty = "n",
  lwd = 2
)
dev.off()

# Boxplot rarefied richness for burned savannas only # FIGURE 3b#

# First merge the treatment names from the abundance dataframe
head(abundance)
head(rareResults)

#### Abundance analysis ####



# Unburned savanna abundance first
abundSav <- abundance[-which(abundance$Treatment %in% c("LF", "LUS")),] 

kruskal.test(x = abundSav$Samples, 
             g = abundSav$Treatment)

dunnTest(Samples ~ Treatment, abundSav)

# Boxplots # FIGURE 3a #
abundSav$numericTreatment <- ifelse(abundSav$Treatment == "LBS", 1, 0)
abundSav$numericTreatment <- ifelse(abundSav$Treatment == "B", 2, abundSav$numericTreatment)
abundSav$numericTreatment <- ifelse(abundSav$Treatment == "MN", 3, abundSav$numericTreatment)
abundSav$numericTreatment <- ifelse(abundSav$Treatment == "MS", 4, abundSav$numericTreatment)


jpeg(filename = "../graphics/numberOfEncounters.jpg",width = 600, height = 600,quality = 100, res = 150)

par(mfrow=c(1,1), mar = c(5,5,2,1))

boxplot(
  Samples ~ numericTreatment,
  abundSav,
  ylab = "Number of samples",
  xlab = "Area",
  axes = F,
  ylim = c(0, 80),
  cex.lab = 0.9
)

axis(1, at = c(1:4), labels = c("LBS","B","MN","MS"))
axis(2, at = seq(0,80,10))
box()

# Show significant differences LBS - MN
arrows(x0 = 1,
       y0 = 78,
       x1 = 3,
       y1 = 78,
       code = 3,
       angle = 90,
       length= 0.05)
text(x = 2, y = 80, substitute(paste(italic('p'), " = 0.026")), cex = 0.5)


mtext(side=3,text="(a)",adj=0,line=1)
dev.off()
