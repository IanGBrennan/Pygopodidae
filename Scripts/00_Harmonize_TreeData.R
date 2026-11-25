library(googlesheets4)

setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae")

gs4_deauth()
#mdata <- read_sheet("https://docs.google.com/spreadsheets/d/1spQTYO9bKFsNZ2rfoNRCdyKFWCamxMxOaAHOCxR1G_I/edit?usp=sharing")
#mdata$genus_species <- paste(mdata$genus,"_",mdata$species,sep="")

mdata <- read.csv("Data/Head_Body_Tail_Nov2026_TreeNames.csv")

################################################################################

# let's match up the HBT data with the squamate phylogeny

#stree <- read.tree("Diet/best_ultrametric_fulltree_ddBD_revision.tre")
stree <- read.tree("Trees/Title_Tree.tre")

overlaps <- intersect(mdata$name_in_tree, stree$tip.label)
missing.zdt <- setdiff(zdata$genus_species, overlaps)
missing.tre <- setdiff(stree$tip.label, overlaps)

#name_in_tree <- NULL
#for(k in 1:nrow(mdata)){
#  if(mdata$genus_species[[k]] %in% overlaps){name_in_tree <- append(name_in_tree, mdata$genus_species[[k]])}
#  if(mdata$gen)
#}

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
mdata <- dplyr::select(mdata, group, family, genus, species, genus_species, name_in_tree, HL, SVL, TL, Source, limbs)

# remove nonconforming species and make two additional "traits"
zdata <- mdata[complete.cases(mdata[,c('SVL','TL')]),]
zdata$SVL <- as.numeric(zdata$SVL); zdata$TL <- as.numeric(zdata$TL)
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
  print(families[[k]])
  curr.fam <- dplyr::filter(zdata, family==families[[k]])
  cutoff.upper <- median(curr.fam$TL.SVL) + (2*sd(curr.fam$TL.SVL))
  cutoff.lower <- median(curr.fam$TL.SVL) - (2*sd(curr.fam$TL.SVL))
  check.hi <- dplyr::filter(curr.fam, TL.SVL > cutoff.upper);
  if(nrow(check.hi)>0){
    check.hi$check <- "high"
    check.hi$diff <- check.hi$TL.SVL - cutoff.upper
  }
  check.lo <- dplyr::filter(curr.fam, TL.SVL < cutoff.lower);
  if(nrow(check.lo)>0){
    check.lo$check <- "low"
    check.lo$diff <- check.lo$TL.SVL - cutoff.lower
  }
  check.into <- rbind(check.into, check.hi, check.lo)
}
# anguids are cool because they have both short tailed and long tailed things!
View(check.into)

# flag the worst records from the SquamBase
bad.sqb <- dplyr::filter(check.into, Source=="Meiri 2024" & check=="high" & group=="Snakes")
# and remove them from our dataset
zdata <- dplyr::anti_join(zdata, bad.sqb)

################################################################################

# MAKE the species data

# Make species averages
ydata <- zdata %>%
  group_by(name_in_tree, group, limbs) %>%
  summarise_at(c('SVL','TL','Total','TL.SVL'), mean) %>%
  mutate(genus = strex::str_before_first(name_in_tree,"_"))

################################################################################

# make matching tree and data

final.overlaps <- intersect(ydata$name_in_tree, stree$tip.label); length(final.overlaps)
final.tree <- drop.tip(stree, setdiff(stree$tip.label,final.overlaps)); Ntip(final.tree)
final.data <- dplyr::filter(ydata, name_in_tree %in% final.overlaps); nrow(final.data)

bt.data <- final.data
bt.tree <- final.tree
# Name the nodes
bt.tree$node.label <- paste("n", (Ntip(bt.tree)+1):(Ntip(bt.tree) + Nnode(bt.tree)), sep="")

# save these as a file we can later load
save(bt.data, bt.tree, file = "Data/BayesTraits/Pygopodidae_BT_Input_Nov2026.RData")

################################################################################

# GENERATE BAYESTRAITS INPUTS

############################################################################

# Generate control files for BayesTraits
source("Scripts/VarRates_Functions.R")

# make a directory to store the trait data
system(paste("mkdir", file.path(getwd(),"Data","BayesTraits")))
system(paste("mkdir", file.path(getwd(),"Data","BayesTraits","Control_Files")))

# let's make control files with trait specific sigma/alpha priors
#trait.names <- names(final.data)
trait.names <- "TL.SVL"

# generate control files for Brownian Motion
for(j in 1:length(trait.names)){
  BTControl(phy = final.tree,
            burnin = "1000000",
            iterations = "10000000",
            sigma = "estimate",
            model = "BM",
            stones = c("100","100000"),
            file.name = paste0("Data/BayesTraits/Control_Files/BM_Control_",trait.names[[j]],".txt"),
            traits = setNames(final.data[[trait.names[[j]]]],final.data$name_in_tree))
}

# generate control files for VarRates
for(j in 1:length(trait.names)){
  BTControl(phy = final.tree,
            burnin = "1000000",
            iterations = "10000000",
            sigma = "estimate",
            model = "VarRates",
            stones = c("100","100000"),
            file.name = paste0("Data/BayesTraits/Control_Files/VarRates_Control_",trait.names[[j]],".txt"),
            traits = setNames(final.data[[trait.names[[j]]]],final.data$name_in_tree))
}

# generate control files for fabric
for(j in 1:length(trait.names)){
  BTControl(phy = final.tree,
            burnin = "1000000",
            iterations = "10000000",
            sigma = "estimate",
            model = "Fabric",
            fabric.limits = "none",
            beta.prior = "VR_LS_BL weibull 1.1 1.5",
            stones = c("100","100000"),
            file.name = paste0("Data/BayesTraits/Control_Files/fabric_Control_",trait.names[[j]],".txt"),
            traits = setNames(final.data[[trait.names[[j]]]],final.data$name_in_tree))
}

# generate control files for Early Burst (Delta)
for(j in 1:length(trait.names)){
  BTControl(phy = final.tree,
            burnin = "1000000",
            iterations = "10000000",
            sigma = "estimate",
            model = "Delta",
            stones = c("100","100000"),
            file.name = paste0("Data/BayesTraits/Control_Files/EB_Control_",trait.names[[j]],".txt"),
            traits = setNames(final.data[[trait.names[[j]]]],final.data$name_in_tree))
}

# generate control files for OU
for(j in 1:length(trait.names)){
  BTControl(phy = final.tree,
            burnin = "1000000",
            iterations = "10000000",
            sigma = "estimate",
            model = "OU",
            stones = c("100","100000"),
            file.name = paste0("Data/BayesTraits/Control_Files/OU_Control_",trait.names[[j]],".txt"),
            traits = setNames(final.data[[trait.names[[j]]]],final.data$name_in_tree))
}

############################################################################

# make a directory to store the trait data
system(paste("mkdir", file.path(getwd(),"Data","BayesTraits","Traits")))

# write each trait to a separate file
#for (k in 1:length(final.data)){
for (k in trait.names){
  trait <- final.data %>% ungroup %>% select(name_in_tree, k)
  tname <- k
  filename <- paste0(file.path(getwd(),"Data","BayesTraits","Traits", tname), "_TRAIT.txt")
  write.table(trait, file=filename, sep=" ", quote=F, col.names=F)
}

############################################################################

writeNexus(final.tree, file="Trees/Title_BT.tre")

# Generate control files for BayesTraits
source("Scripts/VarRates_Functions.R")

# To make life easier, we can create a single shell script to
# run BayesTraits across a series of files in parallel
bt.dir <- "/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/BayesTraits"
trait.files <- dir(file.path(bt.dir,"Traits"), "_TRAIT.txt", full.names = T)
control.files <- dir(file.path(bt.dir,"Control_Files"), full.names = T)

BT.shell.numbers(in.files = trait.files, 
                 ctl.files =  control.files,
                 tree.name = "Trees/Title_BT.tre")
# make sure to adjust the path to the control file to match each trait
# the shell script can be run in parallel if you've got GNU parallel
# $ parallel -j 8 --bar :::: <path_to>/BTShell.sh

############################################################################

source("Scripts/Analysis_Helper_Scripts/fabric_Functions.R")

############################################################################

load("Data/Amphibolurinae_Data.RData")
# if you're here just for plotting, load here:
# load("Data/Amphibolurinae_fabric_AllChainsSummary.RData")

############################################################################

# generate an md5sum key to match against nodes
# any single POST-PROCESSED (FPP) fabric file (using the same tree) will 
# be appropriate and will provide corresponding md5 and taxa lists
nodekey <- fabricMD5("/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/fabric_Runs_NoPriors/Run01/UpperLeg_Run01_TRAIT.txt.VarRates.txt.csv", agam.tree)

############################################################################

# process a single trait and plot the fabric results
ul.res <- fabricProcess(fpp.path = "/Applications/BayesTraitsV4.1.2/Amphibolurinae_Fabric_Traits/Run01/UpperLeg_Run01_TRAIT.txt.VarRates.txt.csv",
                                                phy = agam.tree, trait.name = "Upper Leg")
fabricPlot(process.obj = ul.res, phy = agam.tree, trait.name = "Upper Leg", plot.type="phylogram")
# circles denote shifts in the evolvability (v), red increases, blue decreases
# colored branches denote directional trends (b) in traits, red increases, blue decreases

# or plot it differently:
ul.fab <- fabricPlotVB(process.obj = ul.res, phy = agam.tree, plot.type = "phylogram", 
             line.width = 1, font.size = 0.3, trait = "Upper Leg")

############################################################################

################################################################################

# THIS BELOW SECTION IS MAKING SURE WE ARE 
# (1) REPLACING THE PYGO PORTION OF THE TREE WITH OURS
# (2) NOT THROWING AWAY LOTS OF SPECIES ACCIDENTALLY

hbt.names <- read.csv("Data/Head_Body_Tail_June2026_with_tree_names.csv")
hbt.names.u <- unique(hbt.names$genus_species)

ftree <- read.tree("Trees/Title_Tree.tre")

stree.small <- drop.tip(ftree, tip=unique(hbt.names$name_in_tree))

curr.tree <- stree.small
for (j in 1:length(hbt.names.u)){
  curr.name <- hbt.names.u[[j]]
  which.rows <- which(hbt.names$genus_species==curr.name)
  curr.row <- dplyr::filter(hbt.names, hbt.names$genus_species==curr.name)
  if(!is.na(curr.row$name_in_tree[[1]])){next}
  if(is.na(curr.row$name_in_tree[[1]])){
    options <- curr.tree$tip.label[grep(curr.row$genus[[1]],curr.tree$tip.label)]
    if(length(options)==0){next}
    if(length(options)>0){
      hbt.names[which.rows,"name_in_tree"] <- options[[1]]
      curr.tree <- drop.tip(curr.tree, options[[1]])
    }
  }
}

write.csv(hbt.names, file="Data/Head_Body_Tail_June2026_with_tree_names_2.csv",
          quote=F, row.names=F)



aim.tree <- read.tree("Diet/best_ultrametric_fulltree_ddBD_revision.tre")
pygo <- extract.clade(aim.tree, getMRCA(aim.tree, c("Delma_inornata","Aprasia_smithi")))

aim.tree2 <- drop.tip(aim.tree, pygo$tip.label[1:38])
pygo.mc <- read.tree("~/Desktop/Pygo.tre")
pygo.mc <- drop.tip(pygo.mc, c("Carphodactylidae_Saltuarius_cornutus_A002601",
                               "Pygopodidae_Lialis_burtonis_ABTC70465",
                               "Pygopodidae_Delma_tincta_QMJ80136",
                               "Pygopodidae_Delma_tincta_SAMAR40238",
                               "Pygopodidae_Delma_tincta_WAMR157124",             
                               "Pygopodidae_Delma_tincta_WAMR114391",
                               "Pygopodidae_Aprasia_repens_WAMR127404"))

pygo.mc$edge.length <- pygo.mc$edge.length*100
pygo.mc$root.edge <- NULL
pygo.mc$tip.label <- strex::str_before_last(pygo.mc$tip.label,"_")
pygo.mc$tip.label <- strex::str_after_first(pygo.mc$tip.label,"_")

grep("Ophidiocephalus",aim.tree2$tip.label)
aim.tree2$tip.label[[6785]] <- "Ophidiocephalus_taeniatus_2"

aim.tree3 <- bind.tree(aim.tree2, pygo.mc, where=6785, 
                       position=max(nodeHeights(pygo.mc)))
aim.tree3 <- phytools::force.ultrametric(aim.tree3)
write.tree(aim.tree3, file="~/Desktop/TEST.tre")




# 2. Identify elements that turn into NA when forced to be numeric
non_numeric_matrix <- apply(zdata[,c("SVL","TL")], 2, function(x) is.na(as.numeric(x)))

# 3. Get the exact row and column coordinates
locations <- which(non_numeric_matrix, arr.ind = TRUE)
print(locations)

for(k in 1:nrow(zdata)){
  print(k)
  zdata[k,"SVL"] + zdata[k,"TL"]
}

