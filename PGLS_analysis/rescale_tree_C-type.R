################################################################################
############################# Load required packages ###########################
################################################################################

# Install if needed:
# install.packages(c("ape","phytools","OUwie","geiger","nlme"))

library(ape)        # phylogenetic tree manipulation
library(phytools)   # stochastic character mapping (simmap)
library(OUwie)      # evolutionary models with multiple regimes
library(geiger)     # phylogenetic utilities
library(nlme)       # will be used later for PGLS

################################################################################
############################ Read tree and dataset #############################
################################################################################

tree <- read.tree(here("data", "ultrametric_tree.nwk"))
dat  <- read.csv(here("data","c_type_input_data.csv"), stringsAsFactors = FALSE)

################################################################################
############################# Check species matching ###########################
################################################################################

# Identify species present in the dataset but not in the tree
setdiff(dat$species, tree$tip.label)

# If the result is not character(0), species names must be corrected
# (e.g., underscores, capitalization, spelling differences).


################################################################################
############################# Define evolutionary regimes ######################
################################################################################

# Create a named vector of regimes for stochastic mapping
# Names = species, values = regime
states <- setNames(dat$regime, dat$species)

# Optional: collapse rare regimes (fewer than 5 species) into "Other"
tab <- table(states)
rare <- names(tab)[tab < 5]
states[states %in% rare] <- "Other"

# Update the regime column in the dataset (after collapsing)
dat$regime <- states[match(dat$species, names(states))]


################################################################################
########################### Stochastic mapping (simmap) ########################
################################################################################

# Number of stochastic maps to account for uncertainty in regime transitions
nsim <- 10

# Generate stochastic maps using an equal-rates model (ER)
simmap_trees <- make.simmap(tree, states, model="ER", nsim=nsim)


################################################################################
########################### Prepare data for OUwie #############################
################################################################################

ouwie_dat <- data.frame(
  species = dat$species,
  regime  = dat$regime,
  trait   = dat$gene_count,
  stringsAsFactors = FALSE
)


################################################################################
########################### Estimate regime-specific σ² ########################
################################################################################

# Store model results
results <- vector("list", nsim)

for(i in seq_len(nsim)){

  message("Running OUwie on simmap ", i, "/", nsim)

  # Fit a Brownian Motion model with multiple evolutionary rates (BMS)
  results[[i]] <- tryCatch(
    OUwie(simmap_trees[[i]], ouwie_dat, model="BMS",
          simmap.tree=TRUE, diagn=FALSE),
    error = function(e) e
  )
}


################################################################################
########################### Extract σ² estimates ################################
################################################################################

# Function to extract sigma² from each model
extract_sigma2 <- function(res){

  if(inherits(res, "error")) return(NULL)
  if(is.null(res$solution)) return(NULL)

  return(res$solution["sigma.sq", ])
}

# Apply extraction to all stochastic maps
sigma_list <- lapply(results, extract_sigma2)

# Combine into a matrix (rows = simmaps, columns = regimes)
sigma_mat <- do.call(rbind, sigma_list)

# Display matrix of evolutionary rates
print(sigma_mat)

# Compute the median σ² for each regime
median_sigma <- apply(sigma_mat, 2, median, na.rm=TRUE)

print(median_sigma)


################################################################################
########################### Rescale the phylogenetic tree ######################
################################################################################

# Use one stochastic map as reference for branch segmentation
simt <- simmap_trees[[1]]

# Function to rescale branch lengths using regime-specific σ²
rescale_tree <- function(simmap_tree, sigma_by_regime){

  newtree <- simmap_tree

  # Iterate over all branches
  for(i in seq_len(nrow(simmap_tree$edge))){

    # simmap stores branch segments in $maps
    maps <- simmap_tree$maps[[i]]
    # Example: Diptera=0.12, Other=0.04

    # Multiply each branch segment by its corresponding σ²
    scaled_segments <- mapply(
      function(length, regime) length * sigma_by_regime[regime],
      maps,
      names(maps)
    )

    # New branch length = sum of scaled segments
    newtree$edge.length[i] <- sum(scaled_segments)
  }

  return(newtree)
}

# Apply rescaling
tree_rescaled <- rescale_tree(simt, median_sigma)


################################################################################
########################### Save the rescaled tree #############################
################################################################################

write.tree(tree_rescaled, file="tree_rescaled_c-type_BMS.nwk")