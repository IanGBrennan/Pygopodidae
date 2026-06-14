library(googlesheets4)
library(dplyr)
library(ggplot2)

# one major problem with the SquamBase data:
# max TotalLength - max SVL is not accurate estimate for tails (see notes in SquamBase)
# need to exclude things that don't seem right!

# need to clarify the 'reduced' vs. 'limbless' classifications
# some groups are in both (Delma, Chamaesaura)

################################################################################

# read in the raw data

gs4_deauth()
mdata <- read_sheet("https://docs.google.com/spreadsheets/d/1spQTYO9bKFsNZ2rfoNRCdyKFWCamxMxOaAHOCxR1G_I/edit?usp=sharing")
mdata$genus_species <- paste(mdata$genus,"_",mdata$species,sep="")

################################################################################

# we need to clean up and remove bad records of snakes from the SquamBase dataset

sqb <- dplyr::filter(mdata, Source == "Meiri 2024") # separate the SquamBase records
allelse <- dplyr::filter(mdata, !Source == "Meiri 2024") # separate all other records
# remove the squambase species that are already in the database
sqb <- sqb[-which(sqb$genus_species %in% allelse$genus_species),]
# recombine the datasets
mdata <- rbind(allelse,sqb)

################################################################################

# select the data that interests us
mdata <- dplyr::select(mdata, group, family, genus, species, genus_species, HL, SVL, TL, Source, limbs)

# remove nonconforming species and make two additional "traits"
zdata <- mdata[complete.cases(mdata[,c('SVL','TL')]),]
zdata$Total <- zdata$SVL+zdata$TL
zdata$TL.SVL <- zdata$TL/zdata$SVL

################################################################################

# there are going to be a bunch of species with questionable data and 
# we need to filter those out

# species we probably need to have a second look at:
#View(dplyr::filter(zdata, group=="Snakes" & TL.SVL > 0.75))

families <- unique(zdata$family)
check.into <- NULL
for (k in 1:length(unique(families))){
  curr.fam <- dplyr::filter(zdata, family==families[[k]])
  cutoff.upper <- median(curr.fam$TL.SVL) + (2*sd(curr.fam$TL.SVL))
  cutoff.lower <- median(curr.fam$TL.SVL) - (2*sd(curr.fam$TL.SVL))
  check.hi <- dplyr::filter(curr.fam, TL.SVL > cutoff.upper); 
  check.hi$check <- "high"
  check.hi$diff <- check.hi$TL.SVL - cutoff.upper
  check.lo <- dplyr::filter(curr.fam, TL.SVL < cutoff.lower); 
  check.lo$check <- "low"
  check.lo$diff <- check.lo$TL.SVL - cutoff.lower
  check.into <- rbind(check.into, check.hi, check.lo)
}
# anguids are cool because they have both short tailed and long tailed things!
View(check.into)

# flag the worst records from the SquamBase
bad.sqb <- dplyr::filter(check.into, Source=="Meiri 2024" & check=="high")
# and remove them from our dataset
zdata <- dplyr::anti_join(zdata, bad.sqb)

################################################################################

# VISUALIZE the sample data

# plot the tail/body ratio as a function of body length
ggplot() +
  geom_point(data = transform(zdata, group = NULL), color = "grey85", aes(x=log(SVL), y=TL/SVL)) +
  geom_point(data = zdata, aes(x=log(SVL), y=TL/SVL, color=group)) +
  facet_wrap(~group) + theme_bw() + theme(legend.position = "none")
zdata[which(zdata$TL.SVL == max(zdata$TL.SVL)),]

ggplot() +
  geom_point(data = transform(zdata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = zdata, aes(x=log(Total), y=TL/SVL, color=group)) +
  facet_wrap(~group) + theme_bw() + theme(legend.position = "none")

# plot the tail/body ratio as a linear function sorted by family/group
ggplot() +
  geom_point( data=zdata, aes(x=log(Total), y=TL/SVL, color=family), alpha=0.5, size=0.5) +
  geom_smooth(data=zdata, aes(x=log(Total), y=TL/SVL, color=family), method=lm, se=F, lineend="round") +
  facet_wrap(~group) + theme_bw() + theme(legend.position="none")


# plot the tail/body ratio as a function of total length
ggplot() +
  geom_point(data = transform(zdata, limbs = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = zdata, aes(x=log(Total), y=TL/SVL, fill=limbs), shape=21, size=2) +
  facet_wrap(~limbs) + 
  scale_fill_brewer(palette="Paired") +
  theme_bw()

################################################################################

plot.pygo <- ggplot() +
  geom_point(data = transform(zdata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL), size=2) +
  geom_point(data = dplyr::filter(zdata,family=="Pygopodidae"), aes(x=log(Total), y=TL/SVL,fill=genus),shape=21,size=4) +
  #scale_fill_brewer(palette="Paired") +
  theme_bw() + theme(legend.position = "none")

plot.snake <- ggplot() +
  geom_point(data = transform(zdata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = dplyr::filter(zdata,group=="Snakes"), aes(x=log(Total), y=TL/SVL,fill="red"),shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

plot.lr <- ggplot() +
  geom_point(data = transform(zdata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = dplyr::filter(zdata,!group=="Snakes" & limbs=="reduced" & !family=="Pygopodidae"), 
             aes(x=log(Total), y=TL/SVL),fill="orange",shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

plot.ll <- ggplot() +
  geom_point(data = transform(zdata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = dplyr::filter(zdata,!group=="Snakes" & limbs=="limbless" & !family=="Pygopodidae"), 
             aes(x=log(Total), y=TL/SVL),fill="lightBlue",shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

# establish a plot layout
layout <- "
AAAABB
AAAACC
AAAADD
"
plot.pygo + plot.lr + plot.ll + plot.snake + plot_layout(design=layout)

################################################################################

# MAKE the species data

# Make species averages
ydata <- zdata %>%
  group_by(genus_species, group, family, limbs) %>%
  summarise_at(c('SVL','TL','Total','TL.SVL'), mean) %>%
  mutate(genus = strex::str_before_first(genus_species,"_"))

################################################################################

# VISUALIZE the species data

ggplot() +
  geom_point(data = transform(ydata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = ydata, aes(x=log(Total), y=TL/SVL, color=group)) +
  facet_wrap(~group) + theme_bw() + theme(legend.position = "none")

# plot the tail/body ratio as a linear function sorted by family/group
ggplot() +
#  geom_point(data = transform(ydata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL), size=0.5) +
  geom_point( data=ydata, aes(x=log(Total), y=TL/SVL, color=family), alpha=0.5, size=0.5) +
  geom_smooth(data=ydata, aes(x=log(Total), y=TL/SVL, color=family), method=lm, se=F, lineend="round") +
  facet_wrap(~group) + theme_bw() + theme(legend.position="none")

plot.pygo <- ggplot() +
  geom_point(data = transform(ydata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL), size=2) +
  geom_point(data = dplyr::filter(ydata,family=="Pygopodidae"), aes(x=log(Total), y=TL/SVL,fill=genus),shape=21,size=4) +
  #scale_fill_brewer(palette="Paired") +
  theme_bw() + theme(legend.position = "none")

plot.snake <- ggplot() +
  geom_point(data = transform(ydata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = dplyr::filter(ydata,group=="Snakes"), aes(x=log(Total), y=TL/SVL,fill="red"),shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

plot.limbs <- ggplot() +
  geom_point(data = transform(ydata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = dplyr::filter(ydata, !limbs %in% c("limbless","reduced")), 
             aes(x=log(Total), y=TL/SVL),fill="orange",shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

plot.ll <- ggplot() +
  geom_point(data = transform(ydata, group = NULL), color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data = dplyr::filter(ydata,!group=="Snakes" & limbs%in%c("reduced","limbless") & !family=="Pygopodidae"), 
             aes(x=log(Total), y=TL/SVL),fill="lightBlue",shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")
  
# establish a plot layout
layout <- "
AAAABB
AAAACC
AAAADD
"
plot.pygo + plot.limbs + plot.ll + plot.snake + plot_layout(design=layout)

dplyr::filter(ydata, !limbs %in% c("limbless","reduced") & TL/SVL<0.5 & log(Total) > 5.5)



ggplot(ydata) +
  geom_point(aes(x=log(Total), y=TL/SVL, fill=group), shape=21, size=3) +
  scale_fill_brewer(palette="Paired") +
  theme_bw()

ggplot(ydata) +
  geom_point(aes(x=log(SVL), y=log(TL), color=group, size=TL/SVL)) +
  scale_color_brewer(palette="Paired") +
  theme_bw()

# we can idenfity the center of the TL/SVL and Total length
ydata$TLSVL.c <- (ydata$TL/ydata$SVL) - mean(ydata$TL/ydata$SVL)
ydata$Total.c <- log(ydata$Total) - mean(log(ydata$Total))
ggplot(ydata) +
  geom_point(aes(x=Total.c, y=TLSVL.c, color=group)) +
  scale_color_brewer(palette="Paired") +
  theme_bw()


ydata$TL.SVL <- ydata$TL/ydata$SVL
ydata[which(ydata$TL.SVL == max(ydata$TL.SVL)),]
head(ydata[order(ydata$TL.SVL, decreasing=T),])

################################################################################

# Plotting relationship between tail/body length and total length

# I WANT TO SEPARATE THE SNAKES OUT SO THAT THEY ARE IN A BLOCK
# AND LIZARDS ARE IN A SEPARATE BLOCK

# we can plot the relationship by family and estimate the linear fit
ggplot() +
  geom_point(data = transform(ydata, family = NULL), 
             color = "grey85", aes(x=log(Total), y=TL/SVL)) +
  geom_point(data=ydata, aes(x=log(Total), y=TL/SVL, 
                             fill=family), shape=21, alpha=1) +
  stat_smooth(data=ydata, geom="line", aes(x=log(Total), y=TL/SVL), 
              method=lm, se=F, color="white", alpha=0.75) +
#  scale_color_manual(values=brewer.pal(length(),"Spectral")) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~family)

# loop through to capture the slope and intercept of each lm
families <- unique(ydata$family)
lm.stat <- NULL
for(t in 1:length(families)){
  curr.fam <- dplyr::filter(ydata, family==families[[t]])
  if(nrow(curr.fam)<2){next}
  lm.res <- lm(data=curr.fam, (TL/SVL)~log(Total))
  lm.res <- data.frame(family=families[[t]], slope=lm.res$coefficients[[2]], intercept=lm.res$coefficients[[1]], n=nrow(curr.fam))
  lm.stat <- rbind(lm.stat, lm.res)
}
lm.stat %>%
  filter(n>4) %>%
  arrange(slope)
# Lacertidae, Cordylidae, and Pygopodidae are extreme

###########################################################

###########################################################

###########################################################

###########################################################

###########################################################

# Messing around with the 'Elongation Index' of Title et al.

# OLD: SVL / (2 * sqrt(mass/(SVL*pi)))
# OLD: SVL / 2 * r
# NEW: 
# NEW: 

105/2*sqrt(20.45/(105*3.14))

105/(2*sqrt(20.45/(105*3.14)))

sbase <- read.csv("/Users/ianbrennan/Downloads/geb13812-sup-0003-tables1.csv", check.names=F)
sbase <- dplyr::select(sbase, "Species name (Binomial)", "sub-order", "infraorder", "Family", 
                       "maximum total length (mm)", "maximum SVL (mm)")
colnames(sbase) <- c("Genus_species", "SubOrder", "Infraorder", "Family", "maxTL", "maxSVL")
sbase$Genus_species <- gsub(" ","_",sbase$Genus_species)

sbase2 <- dplyr::filter(sbase, !is.na(maxTL))
sbase3 <- dplyr::filter(sbase2, !is.na(maxSVL))

overlaps <- intersect(sbase3$Genus_species, ydata$genus_species)
noy <- setdiff(sbase3$Genus_species, overlaps)
sbase4 <- dplyr::filter(sbase3, Genus_species %in% noy)
#sbase4$TailLength <- sbase4$maxTL - sbase4$maxSVL

write.csv(sbase4, file="~/Desktop/Squambase_NoHBT.csv")



###########################################################

sbase <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/SquamBase.csv", check.names=F)
sbase <- dplyr::select(sbase, "Species name (Binomial)", "sub-order", "infraorder", "Family", 
                       "maximum total length (mm)", "maximum SVL (mm)", )

all.data <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/Singularity.csv")
all.data2 <- dplyr::filter(all.data, !is.na(completeSVL))
all.data3 <- dplyr::filter(all.data2, !is.na(mass))
title.ei <- dplyr::select(all.data3, treename, family, genus, mass, completeSVL, elongationIndex)
title.ei$OEI <- title.ei$completeSVL/(2*sqrt(title.ei$mass/(title.ei$completeSVL*3.14159)))
title.ei$BW <- 2*sqrt(title.ei$mass/(title.ei$completeSVL*3.14159))


setdiff(mdata$genus_species, all.data$treename)
setdiff(ydata$genus_species, all.data3$treename)

hbt.over <- ydata[which(ydata$genus_species %in% all.data3$treename),]
hbt.over$mass <- sapply(hbt.over$genus_species, function(x) all.data3[which(all.data3$treename==x),"mass"])
hbt.over$old_EIndex <- sapply(hbt.over$genus_species, function(x) all.data3[which(all.data3$treename==x),"elongationIndex"])
hbt.over$EIndex2 <- hbt.over$Total/(2*sqrt(hbt.over$mass/(hbt.over$Total*3.14159)))
hbt.over$widthSVL <- 2*sqrt(hbt.over$mass/(hbt.over$SVL*3.14159))
hbt.over$widthTotal <- 2*sqrt(hbt.over$mass/(hbt.over$Total*3.14159))
hbt.over$radius <- sqrt(hbt.over$mass/(3.14159*hbt.over$SVL))
hbt.over$OEI <- hbt.over$SVL/(2*sqrt(hbt.over$mass/(hbt.over$SVL*3.14159)))

hbt.over <- dplyr::filter(hbt.over, EIndex < 5000)
hbt.over$genus <- sapply(hbt.over$genus_species, function(x) strsplit(x,"_")[[1]][1])

ggplot() +
  geom_point(data=hbt.over, aes(x=EIndex2, y=old_EIndex, color=family))

ggplot() +
  geom_point(data=hbt.over, aes(x=EIndex, y=EIndex2, color=family))


head(hbt.over[order(hbt.over$EIndex,decreasing=T),])
head(hbt.over[order(hbt.over$old_EIndex,decreasing=F),])
head(hbt.over[order(hbt.over$EIndex2,decreasing=T),])


# We want to compare the Title et al. Elongation Index (OEI)
# against a more legitimate estimate of the same index (NEI)
# Title et al. estimated body width by extrapolating diameter of 
# an idealized cylinder assuming volume was mass
# Instead, we actually measured body width directly
oei <- dplyr::select(title.ei, family, genus, treename, mass, completeSVL, OEI)
colnames(oei) <- c("family","genus","genus_species","mass","SVL","OEI")

ndat <- read_sheet("https://docs.google.com/spreadsheets/d/1spQTYO9bKFsNZ2rfoNRCdyKFWCamxMxOaAHOCxR1G_I/edit?usp=sharing")
ndat <- dplyr::select(ndat, family, genus, species, HL, SVL, BW, TL)

ndat$genus_species <- paste(ndat$genus,"_",ndat$species,sep="")

ndat2 <- ndat[complete.cases(ndat[,c('SVL','BW','TL')]),]
# Make species averages
ndat3 <- ndat2 %>%
  group_by(genus_species, family) %>%
  summarise_at(c('SVL','TL','BW'), mean) %>%
  mutate(NEI = SVL/BW, VNEI = (SVL+TL)/BW) 

ggplot(ndat3) +
  geom_point(aes(x=NEI,y=TL/SVL,color=family)) +
  ggtitle("Elongation Index (SVL/BW) vs. TailLength/SVL")
ggplot(ndat3) +
  geom_point(aes(x=VNEI,y=TL/SVL,color=family)) +
  ggtitle("New Elongation Index (TotalL/BW) vs. TailLength/SVL")

testo <- full_join(oei, ndat3, by="genus_species")
testa <- dplyr::select(testo, family.x, genus_species, OEI, NEI)
testa <- testa[complete.cases(testa),]

ggplot(testa) +
  geom_point(aes(x=OEI,y=NEI,color=family.x)) +
  ggtitle("New Elongation Index vs. Old Elongation Index") 
  

