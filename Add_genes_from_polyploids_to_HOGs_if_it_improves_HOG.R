setwd("/my/path/to/Analysis/Analysis_with_corrected_expanded_HOGs/")
load("/my/path/to/A2TEA_finished.RData")

polyploids <- c("Species_A", "Species_B")
diploids <- c("Species_C", "Species_D", "Species_E", "Species_F")

library(dplyr)

# it outputs a data frame per gene with the HOG names and the corrected_HOG names

# If a gene from diploid is in a setdiff, stop checking the subsequent setdiffs and keep the previous set. Do this by not comparing the setdiffs current to previous,
# but always current to set_1 and if the negative condition is met, get the previous comparison (set_prev to set_1) as the best_set 

# ! important: depending on which expansions you are interested in, change the hypothesis you are subsetting accordingly ! #

hypothesis = "hypothesis_3"

get_extended_HOG_info <- function(index, name) {
  expanded_OG <- HYPOTHESES.a2tea[[hypothesis]]@expanded_OGs[[index]]
  setdiff_list <- list()
  
  # Get original set_1 genes
  set_1 <- tryCatch(expanded_OG@add_OG_analysis$set_1@genes$value, error = function(e) NULL)
  if (is.null(set_1)) return(NULL)
  
  # Tag original genes
  original_df <- data.frame(
    gene = set_1,
    gene_type = "original",
    stringsAsFactors = FALSE
  )
  
  # Compute valid setdiffs
  for (j in 2:6) {
    set_curr <- tryCatch(expanded_OG@add_OG_analysis[[paste0("set_", j)]], error = function(e) NULL)
    set_prev <- tryCatch(expanded_OG@add_OG_analysis[[paste0("set_", j - 1)]], error = function(e) NULL)
#    set_1 <- tryCatch(expanded_OG@add_OG_analysis[["set_1"]], error = function(e) NULL)
    # calculate set_diff always to set_1 
    #if (!is.null(set_curr) && !is.null(set_prev)) {
    if (!is.null(set_curr)) {
      genes_diff <- setdiff(set_curr@genes$value, set_1)
      setdiff_list[[paste0("setdiff_", j, "_", "1")]] <- genes_diff
    }
  }
  
  chosen_set <- "set_1"
  added_genes <- character()
  
  # loop over all setdiffs
  for (k in 6:2) {
    diff_name <- paste0("setdiff_", k, "_", "1")
    genes <- setdiff_list[[diff_name]]
    
    # test if only genes from a polyploid species are added
    if (!is.null(genes) && length(genes) > 0) {
      species_info <- HOG_DE.a2tea[HOG_DE.a2tea$gene %in% genes, "species"]
      is_poly <- species_info %in% polyploids
      is_diplo <- species_info %in% diploids
      valid <- all(is_poly & !is_diplo, na.rm = TRUE)
      
      # give the name of the current set
      if (valid) {
        chosen_set <- paste0("set_", k)
        added_genes <- genes
        break
      }
    }
  }
  
  # create a df for all HOGs where there are genes added
  # Only create the added_df if there are added genes
  if (length(added_genes) > 0) {
    added_df <- data.frame(
      gene = added_genes,
      gene_type = "added",
      stringsAsFactors = FALSE
    )
  } else {
    added_df <- data.frame()  # Empty data.frame if no added genes
  }
  
  # Combine original and added genes (if any)
  combined_df <- bind_rows(original_df, added_df)
  
  # Add metadata columns
  result <- combined_df %>%
    left_join(HOG_DE.a2tea %>% select(gene, species) %>% distinct(), by = "gene") %>%
    mutate(
      HOG_name = name,
      best_set = chosen_set,
      new_HOG = paste0("corrected_", name)
    ) %>%
    select(HOG_name, best_set, gene, species, gene_type, new_HOG)
  
  return(result)
}

# Prepare vectors
indices <- seq_along(HYPOTHESES.a2tea[[hypothesis]]@expanded_OGs)
hog_names <- names(HYPOTHESES.a2tea[[hypothesis]]@expanded_OGs)

# Apply function across all HOGs
extended_HOGs_list <- mapply(get_extended_HOG_info, indices, hog_names, SIMPLIFY = FALSE)
extended_HOGs_df <- bind_rows(extended_HOGs_list)

# Export to CSV
write.csv(extended_HOGs_df, file = paste0(hypothesis, "_extended_HOGs_with_original_and_added_genes.csv"), row.names = FALSE)

# this takes several minutes to run, but it works. At least the 2 tested examples are properly corrected.