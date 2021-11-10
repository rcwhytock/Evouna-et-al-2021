# Fidele scripts Dec 2019
#require(xlsx)
require(vegan)
require(labdsv)
library(pairwiseAdonis)

#
fold <- "C:/PhD_2019/Chapitre_2_BIOGEOGRPHIE_TERMITES/Données_Chpt2/Script_Final" #fidele's compute
setwd(fold)
getwd()
dir()

dat <- read.csv2("données/species x sample matrix for Fidele1_25-09-2021.csv")
head(dat) # view data
sites<-read.table("données/sites.txt",header=T,sep="\t") # sites table


# Enlever les sites savanes Mouila
dat <- dat[,-which(names(dat) %in% c("Mouila1", "Mouila2", "Mouila3"))]
sites <- sites[-which(sites$plotname %in% c("Mouila1", "Mouila2", "Mouila3")),]
sites <- droplevels(sites)

# Enlever les sites savanes Tchibanga
dat <- dat[,-which(names(dat) %in% c("Tchib1", "Tchib2", "Tchib3"))]
sites <- sites[-which(sites$plotname %in% c("Tchib1", "Tchib2", "Tchib3")),]
sites <- droplevels(sites)

# Enlever les sites savanes Batéké
dat <- dat[,-which(names(dat) %in% c("BAT_Bon", "BAT_Lec", "BAT_Mve"))]
sites <- sites[-which(sites$plotname %in% c("BAT_Bon", "BAT_Lec", "BAT_Mve")),]
sites <- droplevels(sites)

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
dat3<-merge(sites,df,by.x="plotname1",by.y="plots")

# boxplots
#windows(8,4) # create a window
par(mfrow=c(1,1)) # split the window into 2 panels
# abundance (number of encounters or hits)
levels(dat3$treatment)
dat3$treatment <- factor(dat3$treatment, levels = c("LF", "LUS", "LS"))

jpeg(filename = paste0(fold, "/graphiques/AbondacegradientLope.jpg"),width = 500, height = 500, quality = 100, res = 150)
par(mar = c(5,5,2,1))

boxplot(abund~treatment,dat3,ylab="Number of samples",xlab="Habitat", ylim = c(0,190), axes = F)


axis(1, at = c(1:3), labels = c("LF","LUS","LS"))
axis(2, at = seq(0,190,20), cex.axis = 0.8)
box()
mtext(side=3,text="(a)",adj=0,line=1)

arrows(x0 = 1,
       y0 = 180,
       x1 = 3,
       y1 = 180,
       code = 3,
       angle = 90,
       length= 0.05)
text(x = 2, y = 188, labels = substitute(paste(italic('p'), " = 0.01")), cex = 0.5)

dev.off()

jpeg(filename = paste0(fold, "/graphiques/richessegradientLope.jpg"),width = 500, height = 500, quality = 100, res = 150)
par(mar = c(5,5,2,1))
# number of species
boxplot(rich~treatment,dat3,ylab="Species richness",xlab="Habitat")
mtext(side=3,text="a)",adj=0,line=1)
dev.off()

# statistics
mod1<-lm(abund~treatment,dat3) # abundance
anova(mod1) # ANOVA table
summary(mod1) # summary
TukeyHSD(aov(abund~treatment,dat3)) # pairwise comparisons (post hoc test)

# Non parametric better for these data as 0 inflated and not likely to meet assumptions of normality etc
kruskal.test(x = dat3$abund, 
             g = dat3$treatment) # like ANOVA

dunnTest(abund ~ treatment, dat3) # like Tukeys



mod2<-lm(rich~treatment,dat3) # no of species
anova(mod2) # ANOVA table
summary(mod2) # summary
TukeyHSD(aov(rich~treatment,dat3)) # pairwise comparisons



