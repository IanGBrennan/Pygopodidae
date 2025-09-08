# Double checking that all species with diet data are represented in the tree

nfile <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Squamata_Diets_Families.csv")
ofile <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/dietname-to-treename.csv")
ntree <- read.tree("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/best_ultrametric_fulltree_ddBD_revision.tre")

newies <- setdiff(nfile$predator, ofile$dietname)
data.frame(newies)

in.tree <- data.frame(sapply(newies, function(x) x%in%ntree$tip.label))
colnames(in.tree) <- "in_tree"
in.tree

ntree$tip.label[grep("Dalophia",ntree$tip.label)]; ofile$dietname[grep("Dalophia",ofile$dietname)]
ntree$tip.label[grep("Monopeltis",ntree$tip.label)]; ofile$dietname[grep("Monopeltis",ofile$dietname)]
ntree$tip.label[grep("Dibamus",ntree$tip.label)]; ofile$dietname[grep("Dibamus",ofile$dietname)]
sort(ntree$tip.label[grep("Varanus",ntree$tip.label)])
sort(ntree$tip.label[grep("Liolaemus",ntree$tip.label)]); ofile$dietname[grep("Liolaemus",ofile$dietname)]
sort(ntree$tip.label[grep("Eumeces",ntree$tip.label)]); ofile$dietname[grep("Eumeces",ofile$dietname)]
sort(ntree$tip.label[grep("Ctenophorus",ntree$tip.label)]); ofile$dietname[grep("Ctenophorus",ofile$dietname)]
