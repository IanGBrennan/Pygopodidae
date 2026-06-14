#!/usr/bin/env Rscript
# Match genus_species values in the HBT data to tip labels in the tree.
# If the exact species name exists in the tree, copy it to name_in_tree.
# If not, attempt a synonym lookup and then a fuzzy spelling search.

library(ape)

csv_file <- "Data/Head_Body_Tail_June2026.csv"
tree_file <- "Diet/best_ultrametric_fulltree_ddBD_revision.tre"
syn_file <- "Diet/dietname-to-treename.csv"
out_file <- "Data/Head_Body_Tail_June2026_with_tree_names.csv"

message("Reading CSV: ", csv_file)
df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)

message("Reading tree: ", tree_file)
tree <- read.tree(tree_file)
tips <- tree$tip.label

synonyms <- NULL
if (file.exists(syn_file)) {
  message("Reading synonym map: ", syn_file)
  synonyms <- read.csv(syn_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("dietname", "treename") %in% names(synonyms))) {
    stop("Expected synonym file to contain columns 'dietname' and 'treename'.")
  }
}

if (!"genus_species" %in% names(df)) {
  stop("CSV does not contain a genus_species column.")
}

match_type <- rep(NA_character_, nrow(df))
name_in_tree <- rep(NA_character_, nrow(df))
suggested_tree_names <- rep(NA_character_, nrow(df))

# Exact matches
exact_idx <- which(df$genus_species %in% tips)
name_in_tree[exact_idx] <- df$genus_species[exact_idx]
match_type[exact_idx] <- "exact"

# synonym map matches, if available
if (!is.null(synonyms)) {
  syn_map <- setNames(synonyms$treename, synonyms$dietname)
  syn_idx <- which(is.na(name_in_tree) & df$genus_species %in% names(syn_map))
  if (length(syn_idx) > 0) {
    mapped <- syn_map[df$genus_species[syn_idx]]
    valid <- !is.na(mapped) & mapped %in% tips
    if (any(valid)) {
      fix_idx <- syn_idx[valid]
      name_in_tree[fix_idx] <- mapped[valid]
      match_type[fix_idx] <- "synonym"
    }
  }
}

# fuzzy matching for remaining rows
remaining <- which(is.na(name_in_tree))
for (i in remaining) {
  query <- df$genus_species[i]
  candidates <- tips[agrep(query, tips, max.distance = 0.15, value = FALSE)]
  if (length(candidates) == 0) {
    distances <- adist(query, tips)
    best <- which.min(distances)
    if (length(best) == 1 && distances[best] <= 2) {
      candidates <- tips[best]
    }
  }
  if (length(candidates) == 1) {
    name_in_tree[i] <- candidates
    match_type[i] <- "fuzzy"
  } else if (length(candidates) > 1) {
    suggested_tree_names[i] <- paste(candidates, collapse = ";")
    match_type[i] <- "ambiguous_fuzzy"
  } else {
    match_type[i] <- "missing"
  }
}

# Add columns
out_df <- df
out_df$name_in_tree <- name_in_tree
out_df$match_type <- match_type
out_df$suggested_tree_names <- suggested_tree_names

write.csv(out_df, out_file, row.names = FALSE, quote = TRUE)

message("Wrote output: ", out_file)
message("Summary:")
print(table(match_type, useNA = "ifany"))
