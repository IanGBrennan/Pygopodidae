library(vegan)

setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae")

###############################################################################

# read in the all-data frame from Title et al.
alldat <- read.csv('Diet/alldat.csv')

# read in our phylogenetically smoothed diet proportions
dpn <- read.csv("Diet/diet-proportions-phylo.csv", row.names=1)

# add a family column to the diet data
dpn$family <- sapply(rownames(dpn), function(x) alldat[which(alldat$treename==x),"family"])

# Dietary data
diet <- dpn[,1:32]

# Family
family <- dpn$family

###############################################################################

# Use the Bray-Curtis method for creating a distance matrix from proportional data
diet_dist <- vegdist(diet, method="bray")
# turn the distances into a matrix
diet_mat <- as.matrix(diet_dist)

# Establish what we think is a sensible minimum sample size
min.samp <- 5 # excludes things with less than 5 species

###############################################################################

# Create a function to calculate mean pairwise distances
mean_pairwise <- function(indices, dist_mat) {
  # Fewer than min.samp species = cannot calculate pairwise distance
  if (length(indices) < min.samp) {return(NA_real_)}
  d <- dist_mat[indices, indices, drop = FALSE]
  # Upper triangle, excluding diagonal
  d <- d[upper.tri(d)]
  mean(d, na.rm = TRUE)
}

###############################################################################

# Estimate the mean pairwise distance for each family
family.df <- NULL
for(k in families){
  inds <- which(as.character(dpn$family)==k)
  pair.dist <- mean_pairwise(inds, diet_mat)
  family.df <- rbind(family.df, data.frame(family=k, pair.dist=pair.dist))
}
family.df[order(family.df$pair.dist),]


###############################################################################

# Create a function to generate the permutation test
pairwise_breadth_test <- function(fam1, fam2,
                                  dat, dist_mat,
                                  nperm = 999,
                                  min.spp = 5) {
  
  # Indices of species in each family
  ind1 <- which(as.character(dat$family) == fam1)
  ind2 <- which(as.character(dat$family) == fam2)
  
  # Need at least min.spp species per family
  if (length(ind1) < min.spp || length(ind2) < min.spp) {
    return(data.frame(
      Family_1 = fam1,
      Family_2 = fam2,
      N_1 = length(ind1),
      N_2 = length(ind2),
      Breadth_1 = NA_real_,
      Breadth_2 = NA_real_,
      Difference = NA_real_,
      P = NA_real_
    ))
  }
  
  # Observed breadth
  breadth1 <- mean_pairwise(ind1, dist_mat)
  breadth2 <- mean_pairwise(ind2, dist_mat)
  
  observed_difference <- breadth1 - breadth2
  
  # Combine the species from the two families
  combined <- c(ind1, ind2)
  
  # Original group labels
  labels <- c(
    rep(fam1, length(ind1)),
    rep(fam2, length(ind2))
  )
  
  # Permutation distribution
  perm_difference <- numeric(nperm)
  for (i in seq_len(nperm)) {
    # Shuffle family labels while retaining group sizes
    shuffled_labels <- sample(labels)
    perm_ind1 <- combined[shuffled_labels == fam1]
    perm_ind2 <- combined[shuffled_labels == fam2]
    perm_breadth1 <- mean_pairwise(perm_ind1, dist_mat)
    perm_breadth2 <- mean_pairwise(perm_ind2, dist_mat)
    perm_difference[i] <- (perm_breadth1 - perm_breadth2)
  }
  
  # Two-sided permutation P-value
  p <- (sum(abs(perm_difference) >= abs(observed_difference)) + 1) / (nperm + 1)
  
  data.frame(
    Family_1 = fam1,
    Family_2 = fam2,
    N_1 = length(ind1),
    N_2 = length(ind2),
    Breadth_1 = breadth1,
    Breadth_2 = breadth2,
    Difference = observed_difference,
    P = p
  )
}

###############################################################################

# Make an index of all family pair combinations
families <- unique(as.character(dpn$family))
pairs <- combn(families, 2, simplify = FALSE)

###############################################################################

# Run all the pairwise family comparisons
results <- bind_rows(
  lapply(pairs, function(x) {
    pairwise_breadth_test(
      fam1 = x[1],
      fam2 = x[2],
      dat = dpn,
      dist_mat = diet_mat,
      nperm = 999,
      min.spp = 5
    )
  }
  )
)

###############################################################################

# Correct the p-values for multiple comparisons
results$P_adjusted <- p.adjust(results$P, method = "BH")

# Difference:
# positive means Family 1 has greater dietary breadth
# negative means Family 2 has greater dietary breadth


###############################################################################

# Identify the significant results
significant <- results %>%
  filter(P_adjusted < 0.05)

###############################################################################

# Filter down to just comparisons with the Pygopodidae
results.pygo <- dplyr::filter(results, Family_1=="Pygopodidae" | Family_2=="Pygopodidae")
results.pygo <- results.pygo[complete.cases(results.pygo),]
results.pygo[order(results.pygo$P_adjusted),]

###############################################################################

library(tidyr)
library(ggplot2)

# Reduce the results to exclude NA comparisons
results.clean <- results[complete.cases(results),]

# Order the families phylogenetically for plotting (obnoxious, yes)
families.clean <- c("Pygopodidae",
                    "Diplodactylidae",
                    "Gekkonidae",
                    "Phyllodactylidae",
                    "Sphaerodactylidae",
                    "Scincidae",
                    "Cordylidae",
                    "Amphisbaenidae",
                    "Lacertidae",
                    "Gymnophthalmidae",
                    "Teiidae",
                    "Chamaeleonidae",
                    "Agamidae",
                    "Phrynosomatidae",
                    "Iguanidae",
                    "Dactyloidae",
                    "Tropiduridae",
                    "Liolaemidae",
                    "Leiosauridae",
                    "Varanidae",
                    "Leptotyphlopidae",
                    "Typhlopidae",
                    "Pythonidae",
                    "Boidae",
                    "Uropeltidae",
                    "Pareidae",
                    "Viperidae",
                    "Homalopsidae",
                    "Elapidae",
                    "Lamprophiidae",
                    "Colubridae")


# Create an empty matrix for the p-values
p_mat <- matrix(
  NA_real_,
  nrow = length(families.clean),
  ncol = length(families.clean),
  dimnames = list(families.clean, families.clean)
)

# Create an empty matrix for the difference values
diff_mat <- matrix(
  NA_real_,
  nrow = length(families.clean),
  ncol = length(families.clean),
  dimnames = list(families.clean, families.clean)
)

# Fill both matrices
for (i in seq_len(nrow(results.clean))) {
  
  # choose the families for comparison
  f1 <- results.clean$Family_1[i]
  f2 <- results.clean$Family_2[i]
  
  # extract the parameter estimates
  p <- results.clean$P_adjusted[i]
  d <- results.clean$Difference[i]
  
  # make the P-value is symmetrical
  p_mat[f1, f2] <- p
  p_mat[f2, f1] <- p
  
  # Difference is directional, so reverse the sign
  # when switching Family 1 and Family 2
  diff_mat[f1, f2] <- d
  diff_mat[f2, f1] <- -d
}

# Make the diagonals NA
diag(p_mat) <- NA; 
diag(diff_mat) <- NA; 

# Upper triangle to NA (you might want this for plotting)
# p_mat[upper.tri(p_mat)] <- NA
# diff_mat[upper.tri(diff_mat)] <- NA

# Diagonal = 1 (a family compared with itself) and difference = 0
diag(p_mat) <- 1
diag(diff_mat) <- 0

# make the upper/lower triangles inverse
diff_mat[f2, f1] <- -d

# Put all the data into a long table format for ggplot
heatmap_df <- expand.grid(Family_1 = families.clean, Family_2 = families.clean) %>%
  mutate(
    P_adjusted = mapply(
      function(x, y) p_mat[x, y],
      Family_1,
      Family_2
    ),
    Difference = mapply(
      function(x, y) diff_mat[x, y],
      Family_1,
      Family_2
    )
  ) %>%
  mutate(
    Significant = P_adjusted < 0.05,
    Difference_plot = ifelse(
      Significant,
      Difference,
      NA_real_
    )
  )

###############################################################################

# Plot the results as a heatmap
ggplot(heatmap_df, aes(x = Family_1, y = Family_2, fill = Difference_plot)) +
  geom_tile(colour = "grey90", linewidth = 0.2) +
  coord_fixed() +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    na.value = "white",
    name = "Breadth difference\n(Family 1 − Family 2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
    title = "Pairwise differences in dietary breadth among families",
    subtitle = "Only significant comparisons (BH-adjusted P < 0.05) are colored")
# non-significant pairs are white
# significant comparisons are colored by difference
  # positive (red) means Family 1 has greater dietary breadth
  # negative (blue) means Family 2 has greater dietary breadth

#############