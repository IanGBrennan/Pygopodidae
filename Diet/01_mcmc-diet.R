library(macroevolution) # devtools::install_github("blueraleigh/macroevolution")
setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet")


DATAFILE = "Squamata_Diets_Families.csv"
TREEFILE = "best_ultrametric_fulltree_ddBD_revision.tre"
NAMEMAP = "dietname-to-treename.csv"


names_map = read.csv(NAMEMAP, header=TRUE)
diet = read.csv(DATAFILE, header=TRUE)

diet$predator = names_map[match(diet$predator, names_map[,1]), 2]

phy = read.newick(TREEFILE)
#phy.tre = ape::read.tree(TREEFILE)

common_spp = intersect(tiplabels(phy), unique(diet[,1]))
#diet_not_tree = setdiff(unique(diet$predator), phy.tre$tip.label)

diet = diet[diet$predator %in% common_spp, ]
phy = keep.tip(phy, common_spp)

for (NSTATE in c(20L, 50L, 100L, 1000L))
{

    mcmc = make.rcm.dmm(phy, diet, NSTATE)
    mcmc(niter=2^15, thin=2^4, output.file=sprintf("mcmc-%s.out", NSTATE))

}
