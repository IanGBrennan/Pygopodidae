
setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae")

library(ggplot2)
library(dplyr)

################################################################################

# DATA PREPARATION

################################################################################

mdata <- read.csv("Data/Head_Body_Tail_TreeNames.csv")

# select the data that interests us
mdata <- dplyr::select(mdata, group, family, genus, species, genus_species, SVL, TL, Source, limbs)
mdata$SVL <- as.numeric(mdata$SVL)

# remove nonconforming species and make two additional "traits"
zdata <- mdata[complete.cases(mdata[,c('SVL','TL')]),]
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
  curr.fam <- dplyr::filter(zdata, family==families[[k]])
  cutoff.upper <- median(curr.fam$TL.SVL) + (2*sd(curr.fam$TL.SVL))
  cutoff.lower <- median(curr.fam$TL.SVL) - (2*sd(curr.fam$TL.SVL))
  check.hi <- dplyr::filter(curr.fam, TL.SVL > cutoff.upper); 
  check.hi$check <- "high"
  check.hi$diff <- check.hi$TL.SVL - cutoff.upper
  check.lo <- dplyr::filter(curr.fam, TL.SVL < cutoff.lower); 
  check.lo$check <- "low"
  check.lo$diff <- check.lo$TL.SVL - cutoff.lower
  check.into <- rbind(check.into, check.hi, check.lo)
}
# anguids are cool because they have both short tailed and long tailed things!
View(check.into)

# flag the worst records from the SquamBase
bad.sqb <- dplyr::filter(check.into, Source=="Meiri 2024" & check=="high")
# and remove them from our dataset
zdata <- dplyr::anti_join(zdata, bad.sqb)

# Make species averages
ydata <- zdata %>%
  group_by(genus_species, group, family, limbs) %>%
  summarise_at(c('SVL','TL','Total','TL.SVL'), mean) %>%
  mutate(genus = strex::str_before_first(genus_species,"_"))

################################################################################

# Fit linear models to each family and extract the coefficients

families <- unique(ydata$family)
lm.family <- NULL
for(j in families){
  curr.fam <- dplyr::filter(ydata, family == j)
  if(nrow(curr.fam)<3){next}
  if(nrow(curr.fam)>=3){
    lm.res <- lm(curr.fam$TL.SVL ~ curr.fam$Total)
    r2 <- summary(lm.res)$adj.r.squared
    pval <- summary(lm.res)$coef[2,4]
    std.err <- summary(lm.res)$coef[2,2]
    lm.family <- rbind(lm.family, data.frame(
                            group = curr.fam$group[[1]],
                            family = j,
                            r2 = r2,
                            intercept = lm.res$coefficients[[1]],
                            slope = lm.res$coefficients[[2]],
                            std_err = std.err,
                            p_value = pval,
                            n_spp = nrow(curr.fam)))
    
  }
}
lm.family <- lm.family[order(lm.family$slope),]

# Fit linear models to each group and extract the coefficients

groups <- unique(ydata$group)
lm.group <- NULL
for(k in groups){
  curr.fam <- dplyr::filter(ydata, group == k)
  if(nrow(curr.fam)<3){next}
  if(nrow(curr.fam)>=3){
    lm.res <- lm(curr.fam$TL.SVL ~ curr.fam$Total)
    r2 <- summary(lm.res)$adj.r.squared
    pval <- summary(lm.res)$coef[2,4]
    std.err <- summary(lm.res)$coef[2,2]
    lm.group <- rbind(lm.group, data.frame(
      group = k,
      #family = k,
      r2 = r2,
      intercept = lm.res$coefficients[[1]],
      slope = lm.res$coefficients[[2]],
      std_err = std.err,
      p_value = pval,
      n_spp = nrow(curr.fam)))
  }
}
lm.group <- lm.group[order(lm.group$slope),]

# Save both results to a data file
save(lm.family, lm.group, file="Data/LM_TLSVL_byGroup.RData")
