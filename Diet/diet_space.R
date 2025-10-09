# let's see if we can identify other taxa that are close to our focal
# species in multidimensional diet space

diet.dist <- data.frame(as.matrix(dist(dpn.x[,1:34], method="euclidean", diag=T)))


# make a function to identify nearby points
closePoints <- function(dist.df, sample, hits){
  all.hits <- dist.df[,sample]
  names(all.hits) <- rownames(dist.df)
  top.hits <- head(sort(all.hits),hits)
  return(top.hits)
}

closePoints(diet.dist, "Paradelma_orientalis", 25)
closePoints(diet.dist, "Lialis_burtonis", 25)
closePoints(diet.dist, "Lialis_jicari", 25)
closePoints(diet.dist, "Aprasia_inaurita", 25)
closePoints(diet.dist, "Pletholax_gracilis", 25)
closePoints(diet.dist, "Pygopus_lepidopodus", 25)
closePoints(diet.dist, "Pygopus_nigriceps", 25)

# dietary specialization of Pygopus is likely overstated (quite similar to Delma and even Ophidiocephalus)
# diet of Aprasia similar to blindsnakes
# diet of Lialis similar to colubrids