
# SOME RECORDS ARE DISAPPEARING
# LIKE VARANUS GILLENI, LOOK INTO THIS


# read in the full Cavalcanti dataset
all.diet <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Cavalcanti_Table_S2_imputed.csv")

# select the columns we're interested in
diet <- all.diet[,c(2, #"Lizard_species",
                    3, #"Family",
                    15, #"Reference_Diet",
                    16, #"Journal_Diet",
                    17, #"N_diet
                    grep("VOL",colnames(all.diet)))]

# make a new column to keep the reference information
diet$Reference <- paste(diet$Reference_Diet,diet$Journal_Diet,sep="_")

# read in the key that translates Cavalcanti categories to Title et al.
recat <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Cavalcanti_ReCategorization.csv")
tcats <- unique(recat$Title_cat)
newcats <- data.frame(Genus_species = diet$Lizard_species, Family = diet$Family,
                      n = diet$N_diet, Reference = diet$Reference)
# loop through and summarize the data per row
# translating the volumes to full stomach contents sensu Title et al.
for(j in 1:length(tcats)){
  curr.cat <- dplyr::select(diet, recat[which(recat$Title_cat==tcats[[j]]),"Caval_cat"])
  newcats[,tcats[[j]]] <- data.frame(apply(curr.cat, 1, sum))
}

# round the fields to the closest integer
newcats[,5:29] <- round((newcats[,5:29]/100)*newcats$n,0) 

# pivot the data into longer format to match the Squamata_Diets_Families.csv
nd2 <- newcats %>%
  pivot_longer(
    cols = spiders:misc_vertebrates,
#    names_to = "week",
#    names_prefix = "wk",
    values_to = "count",
    values_drop_na = TRUE
  ) %>%
  filter(., !count==0)

# Combine duplicate rows
nd2.sum <- nd2 %>%
  group_by(Genus_species, Family, name, Reference) %>%
  summarise(total_value = sum(count))

colnames(nd2.sum) <- c("predator","family","prey","source","count")
nd2.sum <- select(nd2.sum, predator, prey, count, family, source)

# save the output
write.csv(nd2.sum, file="/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Cavalcanti_ReCategorization_Longer.csv",
          row.names = F, quote=F)

library(googlesheets4)
gs4_deauth()
full.data <- read_sheet("https://docs.google.com/spreadsheets/d/16pFZgxQ5OdrH1rkYTB5w_TEitdsYNFBIOjNXt-DwDQQ/edit?gid=0#gid=0")

# all the species which aren't in the Title et al database
#cav.new <- dplyr::filter(nd2.sum, predator %in% unique(full.data$predator))
# we need to do this by identifying species that aren't sampled at all
# and for species already sampled, we want to add records of new diet data
# e.g. for Varanus gilleni, previously we had 'lizards' as the only prey item
# with new information we can also add 'orthopterans'
cav.add <- NULL
for(j in 1:nrow(nd2.sum)){
  if(!nd2.sum[j,"predator"] %in% unique(full.data$predator)){cav.add <- rbind(cav.add, nd2.sum[j,])}
  if(nd2.sum[j,"predator"] %in% unique(full.data$predator)){
    curr.full <- dplyr::filter(full.data, predator == nd2.sum[j,"predator"][[1]])
    if(!nd2.sum[j,"prey"] %in% curr.full$prey){cav.add <- rbind(cav.add, nd2.sum[j,])}
    else{next}
  }
}

# save the output
write.csv(cav.add, file="/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Cavalcanti_NewData_ToAdd.csv",
          row.names = F, quote=F)

# compare the old diet data file to the new one to get a list of the species we added
data.final <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Squamata_Diets_Families.csv")
data.prev  <- read.csv("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Squamata_Diets_Families_OLD.csv")

# here's all the newest species
new.names <- data.frame(dietname = setdiff(unique(data.final$predator), unique(data.prev$predator)))
besttree <- ape::read.tree("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/best_ultrametric_fulltree_ddBD_revision.tre")
new.names$treename <- besttree$tip.label[match(new.names$dietname,besttree$tip.label)]
write.csv(new.names, file="/Users/ianbrennan/Documents/GitHub/Pygopodidae/Diet/Cavalcanti_spp_Added.csv")


