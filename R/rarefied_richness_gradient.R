#require(xlsx)
require(vegan)
require(labdsv)
require(FSA) # library(pairwiseAdonis)

dat <- read.csv2("../data/rarefied_species_richness.csv")
head(dat) # view data
str(dat)

# Non parametric better for these data as 0 inflated and not likely to meet assumptions of normality, etc
kruskal.test(x = dat$rarefied.species.richness, 
             g = dat$habitat) # like ANOVA

dunnTest(rarefied.species.richness ~ habitat, dat)

dat$habitatNumeric <- ifelse(dat$habitat == "LUS", 2, 0)
dat$habitatNumeric <- ifelse(dat$habitat == "LF", 1, dat$habitatNumeric)
dat$habitatNumeric <- ifelse(dat$habitat == "LS", 3, dat$habitatNumeric)

jpeg(filename = paste0(fold, "/graphiques/rarefiedRichnessLope.jpg"),width = 500, height = 500, quality = 100, res = 150)
par(mar = c(5,5,2,1))

boxplot(rarefied.species.richness ~ habitatNumeric,xlab = "Habitat", ylab = "Rarefied species richness", dat, ylim = c(0,15), axes = F)

axis(1, at = c(1:3), labels = c("LF","LUS","LS"))
axis(2, at = seq(0,16,2), cex.axis = 0.8)
box()

mtext(side=3,text="(b)",adj=0,line=1)
dev.off()




