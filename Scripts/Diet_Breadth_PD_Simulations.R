library(phytools)


setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet")

################################################################################

# read in the tree
tree <- ape::read.tree("best_ultrametric_fulltree_ddBD_revision.tre")
dpn <- read.csv("diet-proportions-phylo.csv", row.names=1)
dietProportions <- data.matrix(dpn)
dietProportionsStan <- data.frame(sweep(log(dietProportions), 1, rowMeans(log(dietProportions)))) # clr
dpn.pca = prcomp(dietProportionsStan[,1:32])
dpn.x <- data.frame(dpn.pca$x)

# get dietary breadth as the inverse Simpson index
breadths <- data.frame(breadth = apply(dpn, 1, vegan::diversity, index="invsimpson"))

match.tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dpn.x)))


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
all.sim.raw <- NULL
for(y in 1:nsim){
  sims <- NULL
  for(k in 1:32){
    curr.sim <- phytools::fastBM(match.tree, a=thetas[[k]], sig2=sigmas[[k]], bounds=c(0,1))
    # need to standardize them??
    sims <- cbind(sims, curr.sim)
  }
  # translate estimates into proportions per row
  sims <- prop.table(sims, margin=1)
  all.sim.raw[[y]] <- sims
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



# get the breadth for every individual
sim.breadths <- lapply(all.sim.raw, function(x) data.frame(breadth = apply(x, 1, vegan::diversity, index="invsimpson")))
# flatten the simulated breadths into a single dataframe
sim.breadths <- do.call(cbind, sim.breadths)
# get confidence intervals for each species
sim.bqt <- data.frame(t(apply(sim.breadths, 1, quantile, c(0.05,0.95))))
colnames(sim.bqt) <- c("lower","upper")
breadths <- cbind(breadths, sim.bqt)
breadths[which(breadths$breadth<breadths$lower),]
breadths[which(breadths$breadth>breadths$lower),]
breadths[which(breadths$breadth>breadths$upper),]

