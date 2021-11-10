# Fidele scripts Dec 2019
#require(xlsx)
require(vegan)
require(labdsv)
library(pairwiseAdonis)
library(pgirmess)
#
fold <-
  "C:/PhD_2019/Chapitre_2_BIOGEOGRPHIE_TERMITES/Données_Chpt2/Script_Final" #fidele's compute
setwd(fold)
getwd()
dir()

dat <-
  read.csv2("données/species x sample matrix for Fidele1_25-09-2021.csv")
head(dat) # view data
sites <-
  read.table("données/sites.txt", header = T, sep = "\t") # sites table

head(sites) # view data

dput(names(dat))

#Garder les noms despeces
spp <- dat$bar.code.species

# Garder les sites de foret
dat <-
  dat[, which(
    names(dat) %in% c(
      "LOP_Ang",
      "LOP_BoM",
      "LOP_Gdeb",
      "LOP_BV1",
      "LOP_Mab",
      "LOP_BV2",
      "LOP_Fou",
      "LOP_Uap",
      "LOP_Cag",
      "LOP_SEG",
      "LOP_StG"
    )
  )]

dput(sites$treatment)
head(dat)

sites <-
  sites[which(sites$treatment %in% c("LF", "LS", "LUS")),]
sites <- droplevels(sites)
head(sites)

for(i in 1:ncol(dat)){
  dat[which(is.na(dat[,i])),i] <- 0
}

dat

names(dat) <- sites$plotname1
head(dat)

dat2 <- data.frame(dat, spp)
dat1z <- aggregate(. ~ spp, dat2, sum)
dat1 <- dat1z[, -1]

rowSums(dat1)

plots <- names(dat)
rich <- colSums((dat1 > 0) * 1)
abund <- colSums(dat1)
df <-
  data.frame(plots, rich, abund) # dataframe containing richness and abundance
head(df)

# Merge the richness and abundance with the sites table
head(sites)
head(df)

# create a similarity matrix
dis <- vegdist(t(dat1), method = "bray") # similarity matrix

# calculate the mds using the mds() function
n <- nmds(dis, k = 2) # mds

jpeg(filename = paste0(fold, "/graphiques/NMDS_Lope.jpg"),width = 800, height = 800, quality = 100, res = 150)
par(mar = c(4,4,1,1))
# plot
plot(n, type = "n", xlim = c(-.7, .7), ylim = c(-.7, .7), cex = 0.8, cex.axis = 0.8) # produce plot
(nn <-
    names(dat1)) # names of the columns of the species-by-samples matrix
nn <- as.character(sites$MDS)

# quelques pointes sont superposer, ajuster avec 'jitter' (randon noise)
text(jitter(n$points, factor = 400), nn, cex = 0.6) # names of sites/plots
dev.off()

