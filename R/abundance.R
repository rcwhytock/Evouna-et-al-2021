# Fidele scripts Dec 2019
#require(xlsx)
require(vegan)
require(labdsv)
library(pairwiseAdonis)
library(pgirmess)
#
fold <- "C:/PhD_2019/Chapitre_2_BIOGEOGRPHIE_TERMITES/Données_Chpt2/Script_Final" #fidele's compute
setwd(fold)
getwd()
dir()

dat <- read.csv2("données/species x sample matrix for Fidele1_25-09-2021.csv")
head(dat) # view data
sites<-read.table("données/sites.txt",header=T,sep="\t") # sites table

head(sites) # view data

# Enlever les sites de foret
dat <- dat[,-which(names(dat) %in% c("LOP_Ang", "LOP_BoM", "LOP_Gdeb"))]
sites <- sites[-which(sites$plotname %in% c("LOP_Ang", "LOP_BoM", "LOP_Gdeb")),]
sites <- droplevels(sites)
head(sites)

# Enlever les sites savanes non brûlées
dat <- dat[,-which(names(dat) %in% c("LOP_BV1", "LOP_Mab", "LOP_BV2","LOP_Fou"))]
sites <- sites[-which(sites$plotname %in% c("LOP_BV1", "LOP_Mab", "LOP_BV2","LOP_Fou")),]
sites <- droplevels(sites)
head(sites)

# replace NA with 0
for (i in 2:ncol(dat)){
  f<-which(is.na(dat[i]))
  dat[f,i]<-0
}
spp<-dat[,1] # species names
length(spp)
length(unique(spp))

dat1x<-dat[,-1] # remove first column (species names)
#nn<-names(dat1)
names(dat1x)<-sites$plotname1
head(dat1x)

dat2<-data.frame(dat1x,spp)
dat1z<-aggregate(.~spp,dat2,sum)
dat1<-dat1z[,-1]

rowSums(dat1)

#row.names(dat4)<-nn5

plots<-names(dat1x)
rich<-colSums((dat1>0)*1)
abund<-colSums(dat1)
df<-data.frame(plots,rich,abund) # dataframe containing richness and abundance
head(df)

# Merge the richness and abundance with the sites table
head(sites)
head(df)

dat3 <- merge(sites,df,by.x="plotname1",by.y="plots")
head(dat3)

dat3$treatmentNumeric <- ifelse(dat3$treatment == "BS", 2, 0)
dat3$treatmentNumeric <- ifelse(dat3$treatment == "LS", 1, dat3$treatmentNumeric)
dat3$treatmentNumeric <- ifelse(dat3$treatment == "MS", 3, dat3$treatmentNumeric)
dat3$treatmentNumeric <- ifelse(dat3$treatment == "TS", 4, dat3$treatmentNumeric)
head(dat3)

# statistics
mod1<-lm(abund~treatment,dat3) # abundance
anova(mod1) # ANOVA table
summary(mod1) # summary
TukeyHSD(aov(abund~treatment,dat3)) # pairwise comparisons (post hoc test)

# Non parametric better for these data as 0 inflated and not likely to meet assumptions of normality etc
kruskal.test(x = dat3$abund, 
             g = dat3$treatment) # like ANOVA

dunnTest(abund ~ treatment, dat3) # like Tukeys


# boxplots
#windows(8,4) # create a window
jpeg(filename = paste0(fold, "/graphiques/numberOfEncounters.jpg"),width = 500, height = 500,quality = 100, res = 150)
par(mfrow=c(1,1), mar = c(5,5,2,1))
# split the window into 2 panels
# abundance (number of encounters or hits)
levels(dat3$treatment)
#dat3$treatment <- factor(dat3$treatment, levels = c("LS","TS","BS","MS"))
boxplot(abund~treatmentNumeric,dat3,ylab="Number of samples",xlab="Area", axes = F, ylim = c(0,80), cex.lab = 0.9)
axis(1, at = c(1:4), labels = c("LS","BS","MS","TS"))
axis(2, at = seq(0,80,10))
box()

# Show significant differences LS - TS
arrows(x0 = 1,
       y0 = 70,
       x1 = 4,
       y1 = 70,
       code = 3,
       angle = 90,
       length= 0.05)
text(x = 2.5, y = 73, substitute(paste(italic('p'), " = 0.02")), cex = 0.5)

# difference LS - MS
arrows(x0 = 1,
       y0 = 77,
       x1 = 3,
       y1 = 77,
       code = 3,
       angle = 90,
       length= 0.05)
text(x = 2, y = 80, substitute(paste(italic('p'), " = 0.003")), cex = 0.5)

mtext(side=3,text="(a)",adj=0,line=1)
dev.off()

