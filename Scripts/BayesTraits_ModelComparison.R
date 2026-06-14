setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae/Data/BayesTraits/Results/OU/OU_Run01")
stfiles <- dir(getwd(), "Stones.txt", full.names = T)

trait.names <- sapply(sapply(stfiles, function(x) strex::str_after_last(x,"/")), function(y) strex::str_before_last(y,"_TRAIT"))

margLik.df <- NULL
for(k in 1:length(stfiles)){
 stone <- read.delim(stfiles[[k]]) 
 margLik <- as.numeric(stone[nrow(stone),])
 margLik.df <- rbind(margLik.df, data.frame(trait = trait.names[[k]],
                                margLik = margLik))
}
write.csv(margLik.df, file="/Users/ianbrennan/Documents/Science/Seals/Data/BayesTraits/Results/BayesTraits_MarginalLikelihoods_OU_Run03.csv", row.names=F, quote=F)


# 1. Define the palette function
ramp_blu <- colorRampPalette(c("grey", "#BEE2CF", "#8AC1C4", "#559FB9", "#397BA8", "#1B5897", "#213970"))
ramp_red <- colorRampPalette(c("grey", "#FAD4AC", "#F0A881", "#E57A61", "#C75255", "#A72547", "#751133"))


# 2. Generate a specific number of colors
colors_blu <- ramp_blu(100)
colors_red <- ramp_red(100)

# 3. Plot the colors to view them
image(1:200, 1, as.matrix(1:200), col = c(rev(colors_blu),colors_red), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
