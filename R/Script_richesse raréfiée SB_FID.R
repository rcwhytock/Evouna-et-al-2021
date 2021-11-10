#require(xlsx)
require(vegan)
library(pgirmess)
library(FSA)

dat <- read.csv2("../data/rarefied_species_richness.csv")
head(dat) # view data

head(dat)
str(dat) # weird that species richness is imported as a factor...

dat$rarefied.species.richness <- as.numeric(dat$rarefied.species.richness)
dat$habitat <- factor(dat$habitat, levels = c("LS" , "TS", "BS","MS"))
dat$rarefied.species.richness <- as.numeric(dat$rarefied.species.richness)

head(dat)
datComplete <- dat[complete.cases(dat),]

mod1 <- lm(rarefied.species.richness ~ habitat, datComplete)
anov1 <- aov(mod1)
summary(anov1)
TukeyHSD(anov1)

# Non parametric better for these data as 0 inflated and not likely to meet assumptions of normality etc
kruskal.test(x = dat$rarefied.species.richness, 
             g = dat$habitat) # like ANOVA

dunnTest(rarefied.species.richness ~ habitat, dat)

levels(dat$habitat) 

dat$habitatNumeric <- ifelse(dat$habitat == "BS", 2, 0)
dat$habitatNumeric <- ifelse(dat$habitat == "LS", 1, dat$habitatNumeric)
dat$habitatNumeric <- ifelse(dat$habitat == "MS", 3, dat$habitatNumeric)
dat$habitatNumeric <- ifelse(dat$habitat == "TS", 4, dat$habitatNumeric)

jpeg(filename = paste0(fold, "/graphiques/rarefiedRichness.jpg"),width = 500, height = 500, quality = 100, res = 150)
par(mar = c(5,5,2,1))
boxplot(rarefied.species.richness ~ habitatNumeric, xlab = "Area", ylab = "Rarefied species richness", dat, axes = F, cex.lab = 0.9, ylim = c(0,15))
axis(1, at = c(1:4), labels = c("LS","BS","MS","TS"))
axis(2, at = seq(0,16,2), cex.axis = 0.8)
box()
mtext(side=3,text="(b)",adj=0,line=1)

dev.off()

