# pygopods are an enigmatic radiation that fills Australia's invertebrate-eating snake problem

setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet")

# check to make sure Varanus mabitang is included

################################################################################

# DIETARY NICHE BREADTH

################################################################################

# start by identifying what sampling we've got

DATAFILE = "Squamata_Diets_Families.csv"
DATA.OLD = "/Users/ianbrennan/Documents/Science/Title_supplementary_Snake_Singularity/data/3.trait_preparation_diversification/dietModeling/diet.csv"
TREEFILE = "best_ultrametric_fulltree_ddBD_revision.tre"
NAMEMAP = "dietname-to-treename.csv"

names_map = read.csv(NAMEMAP, header=TRUE)
diet = read.csv(DATAFILE, header=TRUE)
diet.old = read.csv(DATA.OLD, header=T)
tree = ape::read.tree(TREEFILE)
alldat = read.csv('alldat.csv')


# keep track of the species that we have added to the diet dataset
new.sp <- setdiff(diet$predator, diet.old$predator)

# check if the new.sp are already in the NAMEMAP
missing.sp <- new.sp[which(!new.sp %in% names_map$dietname)]
# if "character(0) then all the new species are accounted for
# this does not guarantee they are in the tree though!

# check which of the missing.sp can be added easily to the NAMEMAP and do that
add.easy <- missing.sp[which(missing.sp %in% tree$tip.label)]
names_map <- rbind(names_map, data.frame(dietname=add.easy,treename=add.easy))
write.csv(names_map,file="dietname-to-treename.csv", quote=F, row.names=F)
# and which missing.sp need to be hand added because they aren't in the tree
add.hard <- missing.sp[which(!missing.sp %in% tree$tip.label)]
# we'll have to find suitable stand-ins from the tree for these

# check which new.sp are not in the tree you need to find suitable stand-ins from the tree for these
add.harder <- new.sp[which(!new.sp %in% tree$tip.label)]
diet_not_tree = setdiff(unique(diet$predator), tree$tip.label)

# get a report of what group the missing species belong to
snake <- NULL; family <- NULL
for(k in 1:length(diet_not_tree)){
  curr.taxon <- diet_not_tree[[k]]
  if(curr.taxon %in% alldat$treename){
    snake <- append(snake, ifelse(alldat[which(alldat$treename==curr.taxon),"clade_Serpentes"]==1,"snake","lizard"))
    family <- append(family, alldat[which(alldat$treename==curr.taxon),"family"])
  }
  if(!curr.taxon %in% alldat$treename){
    snake <- append(snake, "unknown")
    family <- append(family, "unknown")
  }
}
diet_not_tree <- data.frame(taxon=diet_not_tree, group=snake, family=family)

# add the treename to the diet dataframe
diet$treename = names_map[match(diet$predator, names_map[,1]), 2]

# check if we still have any NAs
subset(diet, is.na(diet$treename))

################################################################################




################################################################################




# kinda following protocol in Figure3.R

diet.data <- read.csv("/Users/ianbrennan/Documents/Science/Title_supplementary_Snake_Singularity/data/3.trait_preparation_diversification/dietModeling/diet-matrix.csv", row.names = 1)

total.items <- rowSums(diet.data)
diet.percent <- diet.data/total.items
testo <- rowSums(diet.percent)

dietProportionsStan <- sweep(log(diet.percent), 1, rowMeans(log(diet.percent))) # clr

diet.pca <- prcomp(diet.percent)

diet.pca$x
plot(diet.pca$x[,2] ~ diet.pca$x[,1])

diet.res <- data.frame(diet.pca$x)
diet.res["Acanthophis_antarcticus",1:3]
diet.res["Lialis_burtonis",1:3]
diet.res["Paradelma_orientalis",1:3]
diet.res["Nucras_intertexta",1:3]


plot.rotations(diet.pca, plot=T, which.pca="prcomp")


#############################################################
source("innovate_elaborate.R")
alldat <- read.csv('alldat.csv')

# establish the groups of interest
snakes <- alldat[which(alldat[, 'Lizard_Snake_Pygopod'] == 1), 'treename']
snakes <- snakes[snakes %in% tree$tip.label]; snakes <- snakes[1:length(snakes)]

lizards <- alldat[which(alldat[, 'Lizard_Snake_Pygopod'] == 0), 'treename']
lizards <- lizards[lizards %in% tree$tip.label]

pygopod <- alldat[which(alldat[, 'Lizard_Snake_Pygopod'] == 2), 'treename']
pygopod <- pygopod[pygopod %in% tree$tip.label]

geckos <- alldat[which(alldat[, 'clade_Gekkota'] == 1), 'treename']
geckos <- geckos[geckos %in% tree$tip.label]

#nolimb <- alldat[which(alldat[, 'numberLimbs'] < 4), 'treename']
nolimb <- alldat[which(alldat[, 'numberDigits'] < 8), 'treename']
nolimb <- nolimb[nolimb %in% tree$tip.label]

# read in the tree
tree <- ape::read.tree("best_ultrametric_fulltree_ddBD_revision.tre")
dpn <- read.csv("diet-proportions-phylo.csv", row.names=1)
dietProportions <- data.matrix(dpn)
# dietProportions <- log(sweep(dietProportions[,-31], 1, dietProportions[,31], "/")) # alr
dietProportionsStan <- data.frame(sweep(log(dietProportions), 1, rowMeans(log(dietProportions)))) # clr
dietProportionsStan$group <- "lizard"
dietProportionsStan[rownames(dietProportionsStan)%in%snakes,"group"] <- "snake"

# get dietary breadth as the inverse Simpson index
breadths <- data.frame(breadth = apply(dpn, 1, vegan::diversity, index="invsimpson"))
plot(density(breadths$breadth))
breadths$family <- sapply(rownames(breadths), function(x) alldat[which(alldat$treename==x),"family"])
breadths$lsp <- sapply(rownames(breadths), function(x) alldat[which(alldat$treename==x),"Lizard_Snake_Pygopod"])

breadth.cat <- rep("limbed lizard", nrow(breadths))
breadth.cat[rownames(breadths) %in% nolimb] <- "limbless lizard"
breadth.cat[rownames(breadths) %in% snakes] <- "snake"
breadth.cat[rownames(breadths) %in% pygopod] <- "pygopod"
breadths$group <- breadth.cat

# devtools::install_github('erocoar/gghalves')
library(gghalves)
install.packages("see")
library(see)

ggplot(data=breadths, aes(x=group, y=breadth, fill=group)) +
  geom_violinhalf(position = position_nudge(x = .2, y = 0), width=2, trim=T) +
  geom_jitter(alpha = 0.1, width = 0.15) +
  stat_summary(fun = "median", geom = "point", color = "blue", size = 3) +
  ylim(c(0,15)) +
  theme_classic() + theme(legend.position = "none")


# Summarize the breadth scores by family
family.breadth <- breadths %>%
  group_by(family) %>%
  summarise(across(c(breadth), list(mean=mean,variance=var), .names = "{fn}"))
# Plot breadth as a function of family
ggplot() +
#  geom_violin(data=breadths, aes(x=family,y=breadth)) +
  geom_point(data=breadths, aes(x=family,y=breadth)) +
  geom_boxplot(data=breadths, aes(x=family,y=breadth)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

####
# what is the relationship between average breadth and variance per family
# is there a way to point out groups that have species with narrow and 
# non-overlapping diets?
# I don't think variance is really telling us what I'm looking for.
####

ggplot(data=family.breadth[complete.cases(family.breadth),],
       aes(x=mean,y=variance,color=family,label=family)) +
  geom_point() +
  geom_text(hjust=0, vjust=0) +
  theme_classic() + theme(legend.position="none")



####################################################

# what if we wanted to determine the average pairwise distance for each *species*

all.dists <- dispRity::pairwise.dist(dpn, diag=T, upper=T)
all.dists <- vegan::vegdist(dpn, method="euclidean", diag=T, upper=T)

m <- matrix(0, nrow = 1442, ncol = 1442)
m[upper.tri(m, diag = F)] <- all.dists
m[lower.tri(m, diag = F)] <- all.dists

mean.dists <- apply(m, 1, mean)

breadths$mean.dist <- mean.dists

ggplot(data=breadths,
       aes(x=breadth,y=mean.dist,color=family)) +
  geom_point() +
  theme_classic() + theme(legend.position="bottom")

####################################################

# If we want to randomly downsample the snake space to match the lizard sampling
# s.df <- dplyr::filter(dietProportionsStan, group=="snake")
# s.df <- dplyr::sample_n(s.df, 309)
# l.df <- dplyr::filter(dietProportionsStan, group=="lizard")
# dietProportionsStan <- rbind(s.df, l.df)

dpn.pca = prcomp(dietProportionsStan[,1:32])
#plot.rotations(dpn.pca, plot=T, which.pca="prcomp")
#autoplot(dpn.pca, loadings=T, loadings.label=T)
dpn.x <- data.frame(dpn.pca$x)
dpn.x$group <- "lizard"; dpn.x[rownames(dpn.pca$x) %in% snakes,"group"] <- "snake"
ggplot(dpn.x) + geom_point(aes(x=PC1, y=PC2, fill=group), shape=21, color="black", size=3) + theme_classic()
summary(dpn.pca)
biplot(dpn.pca)


ptCols <- rep('slateblue1', nrow(dpn.pca$x))
ptCols[rownames(dpn.pca$x) %in% snakes] <- 'coral2'
ptCols[rownames(dpn.pca$x) %in% geckos] <- 'goldenrod1'
ptCols[rownames(dpn.pca$x) %in% pygopod] <- '#61cf61'
dpn.x$colors <- ptCols

table(ptCols)

library(ggplot2)

# plot the spread of diets following the color schemes above
ggplot(dpn.x) +
  #geom_point(aes(x=PC1, y=PC2), fill=dpn.x$colors, shape=21, color="black", size=3) +
  geom_point(aes(x=PC1, y=PC2), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=PC1, y=PC7), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=PC6, y=PC18), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=-7.54,y=0.98),fill="white", shape=21, size=4) +
  #geom_point(aes(x=-0.81,y=-0.09),fill="pink", shape=21, size=4) +
  theme_classic()


fishes <- dietProportionsStan[which(dietProportionsStan$fishes==max(dietProportionsStan$fishes)),]
ggplot() +
  geom_point(data=dpn.x, aes(x=PC1,y=PC2), fill="lightgrey", color="darkgrey", shape=21, size=3) +
  geom_point(data=dplyr::filter(dpn.x, rownames(dpn.x)%in%rownames(fishes)), 
             aes(x=PC1,y=PC2), fill="green", shape=21, size=3) + 
  theme_classic()


library(patchwork)

limbs <- rep('limbed', nrow(dpn.x))
limbs[rownames(dpn.x) %in% nolimb] <- 'reduced'
limbs[rownames(dpn.x) %in% snakes] <- 'snake'
limbs[rownames(dpn.x) %in% pygopod] <- 'pygopod'
limbs[grep("Pygopus",rownames(dpn.x))] <- 'Pygopus'
limbs[grep("Aprasia",rownames(dpn.x))] <- 'Aprasia'
limbs[grep("Delma",rownames(dpn.x))] <- 'Delma'
limbs[grep("Pletholax",rownames(dpn.x))] <- 'Pletholax'
limbs[grep("Ophidiocephalus",rownames(dpn.x))] <- 'Ophidiocephalus'
limbs[grep("Lialis",rownames(dpn.x))] <- 'Lialis'
limbs[grep("Paradelma",rownames(dpn.x))] <- 'Paradelma'

dpn.x$limbs <- limbs

ggplot() +
  geom_point(data = transform(dpn.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2)) +
  geom_point(data = dpn.x, aes(x=PC1, y=PC2, fill=limbs),size=3,shape=21) +
  facet_wrap(~limbs) + theme_bw() + theme(legend.position = "none")

plot.pygo <- ggplot() +
  geom_point(data = transform(dpn.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(dpn.x, limbs %in% c("Aprasia","Delma","Paradelma","Pletholax",
                                                      "Pygopus","Ophidicephalus","Lialis")), aes(x=PC1, y=PC2, fill=limbs),shape=21,size=4) +
  theme_bw() + theme(legend.position = "none")

plot.snake <- ggplot() +
  geom_point(data = transform(dpn.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(dpn.x, limbs=="snake"), aes(x=PC1, y=PC2),shape=21,size=3,fill="darkBlue") +
  theme_bw() + theme(legend.position = "none")

plot.limbs <- ggplot() +
  geom_point(data = transform(dpn.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(dpn.x, limbs=="limbed"), aes(x=PC1, y=PC2, fill=Family),shape=21,size=3,fill="orange") +
  theme_bw() + theme(legend.position = "none")

plot.ll <- ggplot() +
  geom_point(data = transform(dpn.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(dpn.x, limbs=="reduced"), aes(x=PC1, y=PC2, fill=Family),shape=21,size=3,fill="lightBlue") +
  theme_bw() + theme(legend.position = "none")

# establish a plot layout
layout <- "
AAAABB
AAAACC
AAAADD
"
plot.pygo + plot.limbs + plot.ll + plot.snake + plot_layout(design=layout)

#############################################################

ptCols <- rep('goldenrod1', nrow(dpn.pca$x))
ptCols[rownames(dpn.pca$x) %in% nolimb] <- 'lightBlue'
ptCols[rownames(dpn.pca$x) %in% snakes] <- 'coral2'
ptCols[rownames(dpn.pca$x) %in% pygopod] <- '#61cf61'
dpn.x$colors <- ptCols

table(ptCols)

library(ggplot2)
ggplot(dpn.x) +
  #geom_point(aes(x=PC1, y=PC2), fill=dpn.x$colors, shape=21, color="black", size=3) +
  geom_point(aes(x=PC3, y=PC4), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=PC1, y=PC7), fill=dpn.x$colors, shape=21, color="black", size=3) +
  theme_classic()

#dplyr::filter(dpn.x, PC2 < -3 && colors == "")

# inspect the limbless dataset
all.snakes <- snakes[1:2094] # there's something weird with 'snakes' and it won't work
dpn %>%
  dplyr::mutate(taxon = rownames(dpn)) %>%
  dplyr::filter(! taxon %in% all.snakes &
                  taxon %in% nolimb &
                ! taxon %in% pygopod)



#############################################################

diet.tree <- drop.tip(tree, setdiff(tree$tip.label, row.names(dpn.x)))
mvMORPH::mvBM(tree=diet.tree, data=dpn.x[,1:6], model="BM1")
bm.pc1 <- mvMORPH::mvBM(tree=diet.tree, data=dplyr::select(dpn.x,PC1), model="BM1")
squ.thetas <- sapply(1:10, function(x) mvMORPH::mvBM(tree=diet.tree, data=dplyr::select(dpn.x,x), model="BM1")$theta[[1]])

#############################################################

# function to get the euclidean distance between two points
euclidean <- function(a, b) sqrt(sum((a - b)^2))

families <- unique(alldat$family)
dpn.x$Family <- sapply(rownames(dpn.x), function(x) alldat[which(alldat$treename==x),"family"])

fd.fam <- NULL
for(k in 1:length(families)){
  curr.fam <- dplyr::filter(dpn.x, Family == families[[k]])
  nsamps <- nrow(curr.fam)
  if(nsamps<3){next}
  curr.uni <- dplyr::distinct(curr.fam)
  curr.uni <- curr.uni[,1:31]
  if(nsamps<3){next}
  #euclidean(curr.uni[1,], curr.uni[2,])
  fd <- dispRity::func.div(curr.uni)
  fe <- dispRity::func.eve(curr.uni)
  pd <- dispRity::pairwise.dist(curr.uni); pd <- mean(pd)
  fd.fam <- rbind(fd.fam, data.frame(family=families[[k]],fd=fd,fe=fe,pd=pd,no.species=nsamps))
}
plot(fd.fam$fd ~ fd.fam$pd)
plot(fd.fam$fe ~ fd.fam$pd)

# Combine the dietary breadth and pairwise distance data
pair.breadth <- dplyr::full_join(family.breadth, fd.fam)

library(ggrepel)
ggplot(pair.breadth, aes(x=mean, y=pd, color=family, label=family)) +
#  geom_text_repel() +
  geom_point() + xlab("mean breadth") + ylab("mean pairwise distance") +
  geom_text(hjust=0, vjust=0) +
  theme_classic() + theme(legend.position="none")

ggplot(pair.breadth, aes(x=mean, y=fd, color=family, label=family)) +
  geom_point() + xlab("mean breadth") + ylab("functional diversity") +
  geom_text(hjust=0, vjust=0) +
  theme_classic() + theme(legend.position="none")



ggplot(dpn.x) +
  geom_point(aes(x=PC1, y=PC2, fill=Family), shape=21, color="black", size=3) +
  facet_wrap(~Family) +
  theme_classic() + theme(legend.position = "none")

ggplot(dpn.x) +
  geom_point(aes(x=PC3, y=PC4, fill=Family), shape=21, color="black", size=3) +
  facet_wrap(~Family) +
  theme_classic() + theme(legend.position = "none")

ggplot(dpn.x) +
  geom_point(aes(x=PC5, y=PC6, fill=Family), shape=21, color="black", size=3) +
  facet_wrap(~Family) +
  theme_classic() + theme(legend.position = "none")

# Pygopod are just behind Varanids and Colubrids in showing the largest average pairwise distance in diets among squamate families
fd.fam[order(fd.fam$pd,decreasing=T),]
# but are not remarkable in functional diversity or functional evenness (as expected)
fd.fam[order(fd.fam$fd,decreasing=T),]
fd.fam[order(fd.fam$fe,decreasing=T),]

####
match.tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dpn.x)))

#####################################################################################
#####################################################################################
#####################################################################################

# I THINK THIS IS THE WAY TO GO, GETTING THE SIGMA/THETA VALUES FROM A NONSNAKE TREE 
nosnake.tree <- drop.clade(match.tree, tip=c("Anilios_nigrescens","Morelia_viridis"))
nosnake.tree <- drop.tip(nosnake.tree, "NA")
plot(nosnake.tree)
nsim <- 100 # set the number of simulations we want to do

nosnake.dPS <- subset(dietProportionsStan, rownames(dietProportionsStan)%in%nosnake.tree$tip.label)
nosnake.dPS <- subset(dpn, rownames(dpn)%in%nosnake.tree$tip.label)

sigmas <- NULL; thetas <- NULL
for(k in 1:32){
  #  curr.fit <- mvMORPH::mvBM(tree=match.tree, data=dplyr::select(dpn.x, colnames(dpn.x)[[k]]), model="BM1", echo=F, diagnostic=F)
  # curr.fit <- mvMORPH::mvBM(tree=nosnake.tree, data=dplyr::select(nosnake.dPS, colnames(nosnake.dPS)[[k]]), model="BM1", echo=F, diagnostic=F)
  curr.fit <- bounded_bm(tree=nosnake.tree, x=setNames(nosnake.dPS[,k],rownames(nosnake.dPS)),lims=c(0,1),
             levs=200, parallel=T)
  #sigmas <- append(sigmas, curr.fit$sigma)
  sigmas <- append(sigmas, curr.fit$sigsq)
  #thetas <- append(thetas, curr.fit$theta)
  thetas <- append(thetas, curr.fit$x0)
  print(k)
}
all.sim.pca <- NULL
all.sim.tr <- NULL
for(y in 1:nsim){
  sims <- NULL
  for(k in 1:32){
    curr.sim <- phytools::fastBM(match.tree, a=thetas[[k]], sig2=sigmas[[k]], bounds=c(0,1))
    # need to standardize them??
    sims <- cbind(sims, curr.sim)
  }
  # translate estimates into proportions per row
  sims <- prop.table(sims, margin=1)
  sims <- data.frame(sweep(log(sims), 1, rowMeans(log(sims)))) # clr
  colnames(sims) <- colnames(nosnake.dPS)[1:32]
  all.sim.tr[[y]] <- sims
  sim.pca <- data.frame(prcomp(sims)$x)
  all.sim.pca[[y]] <- sim.pca
  print(y)
}

pd.df <- data.frame(family=families)
for(y in 1:nsim){
  pds <- NULL
  for(k in 1:length(families)){
    fam.names <- rownames(dplyr::filter(dpn.x, Family == families[[k]]))
    #  fam.names <- rownames(dplyr::filter(dpn.x, Family == "Pygopodidae"))
    curr.fam <- all.sim.pca[[y]][rownames(all.sim.pca[[y]])%in%fam.names,]
    nsamps <- nrow(curr.fam)
    pd <- dispRity::pairwise.dist(curr.fam); pd <- mean(pd)
    pds <- append(pds, pd)
  }
  pd.df <- cbind(pd.df, pds)
}
#pd.df <- pd.df[complete.cases(pd.df),]
pd.qt <- data.frame(t(apply(pd.df[,2:(nsim+1)], 1, quantile, c(0.05,0.95),na.rm=T)))
colnames(pd.qt) <- c("lower","upper")
pd.qt <- cbind(family=families, pd.qt)

all.res <- dplyr::full_join(pd.qt,fd.fam); View(all.res)
all.res <- all.res[complete.cases(all.res),]

pd.res <- NULL
for(j in 1:nrow(all.res)){
  curr.res <- all.res[j,]
  if(curr.res$pd < curr.res$lower){pd.res<-append(pd.res,"lower");next}
  if(curr.res$pd > curr.res$upper){pd.res<-append(pd.res,"greater");next}
  if(curr.res$pd > curr.res$lower && curr.res$pd < curr.res$upper){pd.res<-append(pd.res,"expected")}
}
all.res$result <- pd.res
all.res$diff.upper <- all.res$pd - all.res$upper
View(all.res)

df$Species <- factor(df$Species, levels = rev(levels(df$Species)))
sort(factor(all.res$family, levels=rev(factor(all.res$family))))

ggplot(all.res) +
  geom_segment( aes(x=family, xend=family, y=lower, yend=upper), color="grey", lineend="round", lwd=2, alpha=0.75) +
  #geom_point( aes(x=family, y=lower), color="grey", size=2 ) +
  #geom_point( aes(x=family, y=upper), color="grey", size=2 ) +
  geom_point( aes(x=family, y=pd, color=result), size=4, alpha=0.75 ) +
  coord_flip() + theme_bw() + theme(legend.position="bottom")


#####################################################################################
#####################################################################################
#####################################################################################

# fit a Brownian Motion model to each individual trait and extract the theta and sigma estimates
sigmas <- NULL; thetas <- NULL
for(k in 1:32){
#  curr.fit <- mvMORPH::mvBM(tree=match.tree, data=dplyr::select(dpn.x, colnames(dpn.x)[[k]]), model="BM1", echo=F, diagnostic=F)
  curr.fit <- mvMORPH::mvBM(tree=match.tree, data=dplyr::select(dietProportionsStan, colnames(dietProportionsStan)[[k]]), model="BM1", echo=F, diagnostic=F)
  sigmas <- append(sigmas, curr.fit$sigma)
  thetas <- append(thetas, curr.fit$theta)
  print(k)
}
# sigmas.all <- sigmas; thetas.all <- thetas

sims <- NULL
for(k in 1:32){
  curr.sim <- phytools::fastBM(match.tree, a=thetas[[k]], sig2=sigmas[[k]])
  sims <- cbind(sims, curr.sim)
}
colnames(sims) <- colnames(dietProportionsStan)[1:32]

sim.pca <- data.frame(prcomp(sims)$x)

pds <- NULL; 
for(k in 1:length(families)){
  fam.names <- rownames(dplyr::filter(dpn.x, Family == families[[k]]))
  #  fam.names <- rownames(dplyr::filter(dpn.x, Family == "Pygopodidae"))
  curr.fam <- sims[rownames(sim.pca)%in%fam.names,]
  nsamps <- nrow(curr.fam)
  pd <- dispRity::pairwise.dist(curr.fam); pd <- mean(pd)
  pds <- append(pds, pd)
}
pd.fam <- cbind(pd.fam, pds)
print(j)

pd.fam <- NULL
for (j in 1:500){
  
  sims <- NULL
  for(k in 1:31){
    curr.sim <- phytools::fastBM(match.tree, a=thetas[[k]], sig2=sigmas[[k]])
    sims <- cbind(sims, curr.sim)
  }
  colnames(sims) <- colnames(dpn.x)[1:31]
  
  pds <- NULL; 
  for(k in 1:length(families)){
    fam.names <- rownames(dplyr::filter(dpn.x, Family == families[[k]]))
    #  fam.names <- rownames(dplyr::filter(dpn.x, Family == "Pygopodidae"))
    curr.fam <- sims[rownames(sims)%in%fam.names,]
    nsamps <- nrow(curr.fam)
    pd <- dispRity::pairwise.dist(curr.fam); pd <- mean(pd)
    pds <- append(pds, pd)
  }
  pd.fam <- cbind(pd.fam, pds)
  print(j)
}

rownames(pd.fam) <- families
pd.fam <- pd.fam[complete.cases(pd.fam),]
pd.quants <- data.frame(t(apply(pd.fam, 1, quantile, c(0.1,0.90)))); 
colnames(pd.quants) <- c("lower","upper")
pd.quants <- pd.quants[sort(rownames(pd.quants)),]
pd.quants$family <- rownames(pd.quants)


testo <- dplyr::full_join(fd.fam, pd.quants)
testo <- testo[complete.cases(testo),]; View(testo)

# get 95% quantiles from the 100 simulations

pygo.names <- rownames(dplyr::filter(dpn.x, Family == "Pygopodidae"))

#############################################################

lizard.dietProportionsStan <- dietProportionsStan[which(rownames(dietProportionsStan) %in% c(lizards,pygopod)),]
liz.pca = prcomp(lizard.dietProportionsStan[,1:34])
plot.rotations(liz.pca, plot=T, which.pca="prcomp")
autoplot(liz.pca, loadings=T, loadings.label=T)
liz.x <- data.frame(liz.pca$x)
plot(liz.x$PC2~liz.x$PC1)
summary(liz.pca)
biplot(liz.pca)


ptCols <- rep('slateblue1', nrow(liz.pca$x))
ptCols[rownames(liz.pca$x) %in% pygopod] <- '#61cf61'
liz.x$colors <- ptCols

table(ptCols)

library(ggplot2)
ggplot(liz.x) +
  geom_point(aes(x=PC1, y=PC2), fill=liz.x$colors, shape=21, color="black", size=3) +
  theme_classic()

#############################################################

diet.liz <- drop.tip(tree, setdiff(tree$tip.label, row.names(liz.x)))
liz.thetas <- mvMORPH::mvBM(tree=diet.liz, data=liz.x[,1:6], model="BM1")

bm.pc1.liz <- mvMORPH::mvBM(tree=diet.liz, data=dplyr::select(liz.x,PC1), model="BM1")
mvMORPH::mvBM(tree=diet.liz, data=dplyr::select(liz.x,PC1), model="BM1", param=list(sigma=0.59))

liz.thetas <- sapply(1:10, function(x) mvMORPH::mvBM(tree=diet.liz, data=dplyr::select(liz.x,x), model="BM1")$theta[[1]])


#############################################################

# function to get the euclidean distance between two points
euclidean <- function(a, b) sqrt(sum((a - b)^2))

families <- unique(alldat$family)
liz.x$Family <- sapply(rownames(liz.x), function(x) alldat[which(alldat$treename==x),"family"])

liz.fam <- NULL
for(k in 1:length(families)){
  curr.fam <- dplyr::filter(liz.x, Family == families[[k]])
  if(nrow(curr.fam)<3){next}
  curr.uni <- dplyr::distinct(curr.fam)
  curr.uni <- curr.uni[,1:31]
  if(nrow(curr.uni)<3){next}
  #euclidean(curr.uni[1,], curr.uni[2,])
  fd <- dispRity::func.div(curr.uni)
  fe <- dispRity::func.eve(curr.uni)
  pd <- dispRity::pairwise.dist(curr.uni); pd <- mean(pd)
  liz.fam <- rbind(liz.fam, data.frame(family=families[[k]],fd=fd,fe=fe,pd=pd))
}
plot(liz.fam$fd ~ liz.fam$pd)
plot(liz.fam$fe ~ liz.fam$pd)

ggplot(liz.x) +
  geom_point(aes(x=PC1, y=PC2, fill=Family), shape=21, color="black", size=3) +
  facet_wrap(~Family) +
  theme_classic() + theme(legend.position = "none")

ggplot(liz.x) +
  geom_point(aes(x=PC3, y=PC4, fill=Family), shape=21, color="black", size=3) +
  facet_wrap(~Family) +
  theme_classic() + theme(legend.position = "none")

ggplot(liz.x) +
  geom_point(aes(x=PC5, y=PC6, fill=Family), shape=21, color="black", size=3) +
  facet_wrap(~Family) +
  theme_classic() + theme(legend.position = "none")

# Pygopods show the largest average pairwise distance in diets among lizard families
liz.fam[order(liz.fam$pd, decreasing=T),]
# but are not exceptional for evenness or functional diversity
liz.fam[order(liz.fam$fd, decreasing=T),]
liz.fam[order(liz.fam$fe, decreasing=T),]


#############################################################

limbs <- rep('limbed', nrow(liz.x))
limbs[rownames(liz.x) %in% nolimb] <- 'reduced'
limbs[rownames(liz.x) %in% pygopod] <- 'pygopod'
liz.x$limbs <- limbs

ggplot() +
  geom_point(data = transform(liz.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2)) +
  geom_point(data = liz.x, aes(x=PC1, y=PC2, color=limbs),size=3) +
  facet_wrap(~limbs) + theme_bw() + theme(legend.position = "none")

plot.pygo <- ggplot() +
  geom_point(data = transform(liz.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(liz.x, limbs=="pygopod"), aes(x=PC1, y=PC2, fill=limbs),shape=21,size=4) +
  theme_bw() + theme(legend.position = "none")

plot.snake <- ggplot() +
  geom_point(data = transform(liz.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(liz.x, limbs=="snake"), aes(x=PC1, y=PC2, fill=limbs),shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

plot.limbs <- ggplot() +
  geom_point(data = transform(liz.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(liz.x, limbs=="limbed"), aes(x=PC1, y=PC2, fill=limbs),shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

plot.ll <- ggplot() +
  geom_point(data = transform(liz.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = dplyr::filter(liz.x, limbs=="reduced"), aes(x=PC1, y=PC2, fill=limbs),shape=21,size=3) +
  theme_bw() + theme(legend.position = "none")

# establish a plot layout
layout <- "
AAAABB
AAAACC
AAAADD
"
plot.pygo + plot.limbs + plot.ll + plot.snake + plot_layout(design=layout)

ggplot() +
  geom_point(data = transform(liz.x, limbs = NULL), color = "grey85", aes(x=PC1, y=PC2),size=2) +
  geom_point(data = liz.x, aes(x=PC1, y=PC2, fill=limbs),shape=21,size=3) +
  facet_wrap(~limbs) + theme_bw() + theme(legend.position = "none")

#############################################################

gecko.dietProportionsStan <- dietProportionsStan[which(rownames(dietProportionsStan) %in% c(geckos)),]
gek.pca = prcomp(gecko.dietProportionsStan)
plot.rotations(gek.pca, plot=T, which.pca="prcomp")
autoplot(gek.pca, loadings=T, loadings.label=T)
gek.x <- data.frame(gek.pca$x)
plot(gek.x$PC2~gek.x$PC1)
summary(gek.pca)
biplot(gek.pca)


ptCols <- rep('orangered', nrow(gek.pca$x))
ptCols[rownames(gek.pca$x) %in% pygopod] <- '#61cf61'
gek.x$colors <- ptCols

table(ptCols)

library(ggplot2)
ggplot(gek.x) +
  geom_point(aes(x=PC1, y=PC2), fill=gek.x$colors, shape=21, color="black", size=3) +
  theme_classic()















#############################################################



library(ggfortify)
df <- iris[1:4]
pca_res <- prcomp(df, scale. = TRUE)

autoplot(pca_res)






dpt <- read.csv("/Users/ianbrennan/Documents/Science/Title_supplementary_Snake_Singularity/data/3.trait_preparation_diversification/dietModeling/diet-proportions-nonphylo.csv", row.names=1)
dpt.m <- data.matrix(dpt)
dpt.s <- sweep(log(dpt.m), 1, rowMeans(log(dpt.m)))
dpt.pca = prcomp(dpt.s)
plot.rotations(dpt.pca, plot=T, which.pca="prcomp")
dpt.x <- data.frame(dpt.pca$x)
plot(dpt.x$PC2~dpt.x$PC1)
summary(dpt.pca)



diet.tb <- NULL
for(j in 1:nrow(diet.data)){
  curr.tax <- unlist(diet.data[j,])
  curr.diet <- curr.tax[which(curr.tax>0)]
  for(y in 1:length(curr.diet)){
    diet.tb <- rbind(diet.tb, data.frame(predator=rownames(diet.data)[[j]],
                                         prey=names(curr.diet)[[y]],
                                         count=curr.diet[[y]]))
  }
}
write.csv(diet.tb, file="~/Desktop/Squamata_Diets.csv", quote=F)

pygo.names <- c("Delma_impar",
                "Delma_inornata",
                "Delma_fraseri",
                "Delma_nasuta",
                "Lialis_burtonis",
                "Lialis_jicari",
                "Pygopus_lepidopodus",
                "Pygopus_nigriceps",
                "Aprasia_haroldi",
                "Aprasia_inaurita",
                "Aprasia_parapulchella",
                "Aprasia_pulchella",
                "Aprasia_repens",
                "Aprasia_smithi",
                "Aprasia_striolata",
                "Paradelma_orientalis",
                "Ophidocephalus_taeniatus")

tt <- read.tree("/Users/ianbrennan/Desktop/Pygopodidae_Diet/best_ultrametric_fulltree_ddBD_revision.tre")

pygo.names%in%tt$tip.label

# make a function to visualize the contribution of traits to the PCs
# this function assumes the pca.obj is the output of 
plot.rotations <- function(pca.obj, plot=T, which.pca=c("prcomp","princomp","lda")){
  require(reshape2)
  require(ggplot2)
  
  ## PRCOMP
  if(which.pca=="prcomp"){
    # make the pca rotation matrix a dataframe
    pca.rot <- data.frame(pca.obj$rotation)
    # change the column names to reflect the relative contributions to the variance  
    colnames(pca.rot) <- sapply(1:ncol(pca.rot), function(x) paste0(colnames(pca.rot)[[x]],"_",round(summary(pca.obj)$importance[2,x]*100,2)))
    # create a variable with the rownames
    pca.rot$trait <- rownames(pca.rot)
    # reshape the data from wide to long format
    pca.rot <- reshape2::melt(pca.rot, id.vars="trait", variable.name="PC")
    # plot the resulting data
    rot.plot <- ggplot(pca.rot,aes(x=trait,y=value)) +
      geom_col() + theme_bw() +
      facet_wrap(~PC) +
      coord_flip() + scale_x_discrete(limits=rev)
  }
  
  ## PRINCOMP
  if(which.pca=="princomp"){
    # make the pca rotation matrix a dataframe
    pca.rot <- data.frame(loadings(pca.obj)[])
    # change the column names to reflect the relative contributions to the variance 
    cum.var <- pca.obj$sdev^2/sum(pca.obj$sdev^2)
    colnames(pca.rot) <- sapply(1:ncol(pca.rot), function(x) paste0(colnames(pca.rot)[[x]],"_",round(cum.var[[x]]*100,2)))
    # create a variable with the rownames
    pca.rot$trait <- rownames(pca.rot)
    # reshape the data from wide to long format
    pca.rot <- reshape2::melt(pca.rot, id.vars="trait", variable.name="Comp.")
    # plot the resulting data
    rot.plot <- ggplot(pca.rot,aes(x=trait,y=value)) +
      geom_col() + theme_bw() +
      facet_wrap(~Comp.) +
      coord_flip() + scale_x_discrete(limits=rev)
  }
  
  ## MASS::LDA
  if(which.pca=="lda"){
    # make the pca rotation matrix a dataframe
    lda.rot <- data.frame(pca.obj$scaling)
    prop.trace <- (pca.obj$svd^2)/sum(pca.obj$svd^2)
    # change the column names to reflect the relative contributions to the variance  
    colnames(lda.rot) <- sapply(1:ncol(lda.rot), function(x) paste0(colnames(lda.rot)[[x]],"_",round(prop.trace[[x]]*100,2)))
    # create a variable with the rownames
    lda.rot$trait <- rownames(lda.rot)
    # reshape the data from wide to long format
    lda.rot <- reshape2::melt(lda.rot, id.vars="trait", variable.name="LD")
    # plot the resulting data
    rot.plot <- ggplot(lda.rot,aes(x=trait,y=value)) +
      geom_col() + theme_bw() +
      facet_wrap(~LD) +
      coord_flip() + scale_x_discrete(limits=rev)
  }
  
  if(plot==T){print(rot.plot)}
  return(rot.plot)
}








sq <- read.csv("/Users/ianbrennan/Desktop/Pygopodidae_Diet/Squamata_Diets_Families.csv")
dttr <- read.csv("/Users/ianbrennan/Desktop/Pygopodidae_Diet/dietname-to-treename.csv")
fazza <- NULL
for(j in 1:nrow(sq)){
  if(sq$predator[[j]] %in% dttr$dietname){
    curr.name <- dttr[which(dttr$dietname==sq$predator[[j]]),"treename"]
    if(curr.name %in% alldat$treename){fazza[[j]] <- alldat[which(alldat$treename==curr.name),"family"]}
    if(!curr.name %in% alldat$treename){fazza[[j]] <- "NOT_IN_ALLDAT"}
    }
  #if(sq$predator[[j]] %in% alldat$treename){fazza[[j]]<-alldat[which(alldat$treename==sq$predator[[j]]),"family"]}
  if(!sq$predator[[j]] %in% dttr$dietname){fazza[[j]]<-"UNKNOWN"}
}
sq$family <- unlist(fazza)
#sq$family <- sapply(sq$predator, function(x) alldat[which(alldat$treename==x),"family"])
write.csv(sq, file="/Users/ianbrennan/Desktop/Pygopodidae_Diet/Squamata_Diets_Families.csv", row.names=F, quote=F)

which(class(sq$family)=="list")

sqf <- read.csv("/Users/ianbrennan/Desktop/Pygopodidae_Diet/Squamata_Diets_Families.csv")
overlaps <- intersect(sqf$predator, tree$tip.label)
getdrop <- data.frame(lost=setdiff(sqf$predator, overlaps))
getdrop$genus <- sapply(getdrop$lost, function(x) strsplit(x,"_")[[1]][1])


dttr <- read.csv("/Users/ianbrennan/Desktop/Pygopodidae_Diet/dietname-to-treename.csv")
setdiff(sq$predator, dttr$dietname)

setdiff(sq$predator, rownames(dpn))

