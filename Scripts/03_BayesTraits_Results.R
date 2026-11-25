setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae")

############################################################################

source("Scripts/fabric_Functions.R")

############################################################################

load("Data/BayesTraits/Pygopodidae_BT_Input.RData")
load("Data/Ancestral_Trait_Estimates.RData")
# if you're here just for plotting, load here:
# load("Data/Amphibolurinae_fabric_AllChainsSummary.RData")

emp <- data.frame(TL.SVL=bt.data$TL.SVL); rownames(emp) <- bt.data$name_in_tree
emp.anc <- rbind(emp, anc.df)

############################################################################

# generate an md5sum key to match against nodes
# any single POST-PROCESSED (FPP) fabric file (using the same tree) will 
# be appropriate and will provide corresponding md5 and taxa lists
nodekey <- fabricMD5("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/BayesTraits/Results/fabric_Run01/TL.SVL_TRAIT.txt.VarRates.txt.csv", bt.tree)

############################################################################

# process a single trait and plot the fabric results
tail.res <- fabricProcess(fpp.path = "/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/BayesTraits/Results/fabric_Run01/TL.SVL_TRAIT.txt.VarRates.txt.csv",
                        phy = bt.tree, trait.name = "TL.SVL")
fabricPlot(process.obj = tail.res, phy = bt.tree, trait.name = "TL/SVL", plot.type="phylogram")
# circles denote shifts in the evolvability (v), red increases, blue decreases
# colored branches denote directional trends (b) in traits, red increases, blue decreases

# or plot it differently:
tail.fab <- fabricPlotVB(process.obj = tail.res, phy = bt.tree, plot.type = "fan", 
                       line.width = 1, font.size = 0.1, trait = "TL/SVL")
# circles denote directional trends (b) in traits, red bigger, blue smaller
# colored branches denote shifts int eh evolvability (v), red increases, blue decreases

# or plot the path from the MRCA to a clade of interest
# here colored edges denote rate changes (background is grey), indicated by a white circle
# and arrows show directional changes
fabricTraitgram(process.obj = chain.sum, phy = bt.tree, trait.name = "TL.SVL",
                trait = emp.anc, tip.spread = c("Takydromus_formosanus","Takydromus_sexlineatus"),
                focus = "clade")
# compared to Takydromus which show a trend towards just longer tails
# pygopods show a variety of trends
fabricTraitgram(process.obj = chain.sum, phy = bt.tree, trait.name = "TL.SVL",
                trait = emp.anc, tip.spread = c("Delma_fraseri","Aprasia_smithi"),
                focus = "clade")



############################################################################


# identify all the fabric directories
fab.dir <- dir("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/BayesTraits/Results/fabric", 
               pattern="Run", include.dirs=T, full.names=T)

# loop across each directory and summarize the run information
# NOTE: you need to have used the fabric post processor on these runs already!
all.chains <- NULL
for (p in 1:length(fab.dir)){
  
  # identify all the fabric files
  fab.files <- dir(fab.dir[[p]], ".csv", full.names = T) 
  
  # extract the trait names
  trait.names <- sapply(sapply(fab.files, function(x) strex::str_after_last(x, "/")), function(y) strex::str_before_first(y, "_"))
  
  # loop through each but keep a joint output of all the b & v shifts
  full.shifts <- NULL; ancestors <- NULL
  for (j in 1:length(fab.files)){
    cur.res <- fabricProcess(fpp.path = fab.files[[j]],
                             phy = bt.tree, trait.name = trait.names[[j]])
    full.shifts <- rbind(full.shifts, cur.res)
  }
  # specify which chain the result comes from
  full.shifts$chain <- p
  # combine all shifts results with other chains
  all.chains <- rbind(all.chains, full.shifts)
}

# summarize the shifts into data about individual branches
chains.summary <- NULL
unique.combo <- expand(all.chains, nesting(node, parameter, trait))
for (h in 1:nrow(unique.combo)){
  curr.combo <- dplyr::filter(all.chains, 
                              node==unique.combo$node[[h]],
                              parameter==unique.combo$parameter[[h]],
                              trait==unique.combo$trait[[h]])
  curr.res <- data.frame(node = curr.combo$node[[1]],
                         scale.mean = mean(curr.combo$scale),
                         scale.sd   = sd(curr.combo$scale),
                         parameter = curr.combo$parameter[[1]],
                         trait = curr.combo$trait[[1]],
                         #chains = paste0(min(curr.combo$chain),":",max(curr.combo$chain))
                         node.parent = curr.combo$node.parent[[1]],
                         chains = nrow(curr.combo))
  chains.summary <- rbind(chains.summary, curr.res)
}
# each row of 'chains.summary' will correspond to one branch for one trait

# exclude shifts from the posterior that weren't found in all chains
chain.sum <- dplyr::filter(chains.summary, chains == 2)

# save these data to file
save(all.chains, chains.summary, chain.sum,
     file = "Data/BayesTraits/Pygopodidae_fabric_AllChainsSummary.RData")
