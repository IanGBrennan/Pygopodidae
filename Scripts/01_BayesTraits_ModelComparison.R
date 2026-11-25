
setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/BayesTraits/Results")
############################################################################

# get all the directories
dir.names <- list.dirs(full.names=F)

all.marg <- NULL
for(j in 2:length(dir.names)){
  # set the current directory
  curr.dir <- file.path(getwd(),dir.names[j])
  
  # grab all the stones files
  stfiles <- dir(curr.dir, "Stones.txt", full.names = T)
  
  # extract the names of the traits
  trait.names <- sapply(sapply(stfiles, function(x) strex::str_after_last(x,"/")), function(y) strex::str_before_last(y,"_TRAIT"))
  
  # make a dataframe to hold the marginal likelihoods of the runs
  margLik.df <- NULL
  for(k in 1:length(stfiles)){
    stone <- read.delim(stfiles[[k]]) 
    margLik <- as.numeric(stone[nrow(stone),])
    margLik.df <- rbind(margLik.df, data.frame(trait = trait.names[[k]],
                                               margLik = margLik))
  }
  
  if(j==1){all.marg <- margLik.df}
  if(j>1){all.marg <- cbind(all.marg, margLik.df[,2])}
}
# rename the results
colnames(all.marg) <- dir.names[2:length(dir.names)]

# write the results to a file
write.csv(all.marg, file=file.path(getwd(),"BayesTraits_MarginalLikelihoods_RESULTS.csv"), row.names=F, quote=F)

############################################################################

