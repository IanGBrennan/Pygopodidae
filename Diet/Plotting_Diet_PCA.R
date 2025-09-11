# pygopods are an enigmatic radiation that fills Australia's invertebrate-eating snake problem

setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet")

################################################################################

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

# keep track of the species that we have added to the diet dataset
new.sp <- setdiff(diet$predator, diet.old$predator)

# check if the new.sp are already in the NAMEMAP
missing.sp <- new.sp[which(!new.sp %in% names_map$dietname)]

# check which of the missing.sp can be added easily to the NAMEMAP and do that
add.easy <- missing.sp[which(missing.sp %in% tree$tip.label)]
names_map <- rbind(names_map, data.frame(dietname=add.easy,treename=add.easy))
write.csv(names_map,file="dietname-to-treename.csv", quote=F, row.names=F)
# and which missing.sp need to be hand added because they aren't in the tree
add.hard <- missing.sp[which(!missing.sp %in% tree$tip.label)]
# we'll have to find suitable stand-ins from the tree for these

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
dpn <- read.csv("diet-proportions-phylo.csv", row.names=1)
dietProportions <- data.matrix(dpn)
# dietProportions <- log(sweep(dietProportions[,-31], 1, dietProportions[,31], "/")) # alr
dietProportionsStan <- data.frame(sweep(log(dietProportions), 1, rowMeans(log(dietProportions)))) # clr
dietProportionsStan$group <- "lizard"
dietProportionsStan[rownames(dietProportionsStan)%in%snakes,"group"] <- "snake"

# If we want to randomly downsample the snake space to match the lizard sampling
# s.df <- dplyr::filter(dietProportionsStan, group=="snake")
# s.df <- dplyr::sample_n(s.df, 309)
# l.df <- dplyr::filter(dietProportionsStan, group=="lizard")
# dietProportionsStan <- rbind(s.df, l.df)

dpn.pca = prcomp(dietProportionsStan[,1:31])
#plot.rotations(dpn.pca, plot=T, which.pca="prcomp")
#autoplot(dpn.pca, loadings=T, loadings.label=T)
dpn.x <- data.frame(dpn.pca$x)
dpn.x$group <- "lizard"; dpn.x[rownames(dpn.pca$x) %in% snakes,"group"] <- "snake"
ggplot(dpn.x) + geom_point(aes(x=PC1, y=PC2, fill=group), shape=21, color="black", size=3) + theme_classic()
summary(dpn.pca)
biplot(dpn.pca)





alldat <- read.csv('alldat.csv')
tree <- ape::read.tree("best_ultrametric_fulltree_ddBD_revision.tre")

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


ptCols <- rep('slateblue1', nrow(dpn.pca$x))
ptCols[rownames(dpn.pca$x) %in% snakes] <- 'coral2'
ptCols[rownames(dpn.pca$x) %in% geckos] <- 'goldenrod1'
ptCols[rownames(dpn.pca$x) %in% pygopod] <- '#61cf61'
dpn.x$colors <- ptCols

table(ptCols)

library(ggplot2)
ggplot(dpn.x) +
  #geom_point(aes(x=PC1, y=PC2), fill=dpn.x$colors, shape=21, color="black", size=3) +
  geom_point(aes(x=PC1, y=PC2), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=PC1, y=PC7), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=PC6, y=PC18), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=-7.54,y=0.98),fill="white", shape=21, size=4) +
  #geom_point(aes(x=-0.81,y=-0.09),fill="pink", shape=21, size=4) +
  theme_classic()

#############################################################

ptCols <- rep('pink', nrow(dpn.pca$x))
ptCols[rownames(dpn.pca$x) %in% nolimb] <- 'lightBlue'
ptCols[rownames(dpn.pca$x) %in% snakes] <- 'coral2'
ptCols[rownames(dpn.pca$x) %in% pygopod] <- '#61cf61'
dpn.x$colors <- ptCols

table(ptCols)

library(ggplot2)
ggplot(dpn.x) +
  geom_point(aes(x=PC1, y=PC2), fill=dpn.x$colors, shape=21, color="black", size=3) +
  #geom_point(aes(x=PC1, y=PC3), fill=dpn.x$colors, shape=21, color="black", size=3) +
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
  if(nrow(curr.fam)<3){next}
  curr.uni <- dplyr::distinct(curr.fam)
  curr.uni <- curr.uni[,1:31]
  if(nrow(curr.uni)<3){next}
  #euclidean(curr.uni[1,], curr.uni[2,])
  fd <- dispRity::func.div(curr.uni)
  fe <- dispRity::func.eve(curr.uni)
  pd <- dispRity::pairwise.dist(curr.uni); pd <- mean(pd)
  fd.fam <- rbind(fd.fam, data.frame(family=families[[k]],fd=fd,fe=fe,pd=pd))
}
plot(fd.fam$fd ~ fd.fam$pd)
plot(fd.fam$fe ~ fd.fam$pd)

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

# Pygopods show the largest average pairwise distance in diets among squamate families
fd.fam[order(fd.fam$pd,decreasing=T),]
# but are not remarkable in functional diversity or functional evenness (as expected)
fd.fam[order(fd.fam$fd,decreasing=T),]
fd.fam[order(fd.fam$fe,decreasing=T),]

#############################################################

lizard.dietProportionsStan <- dietProportionsStan[which(rownames(dietProportionsStan) %in% c(lizards,pygopod)),]
liz.pca = prcomp(lizard.dietProportionsStan)
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

