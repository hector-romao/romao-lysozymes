################################################################################
############################ Load required packages ############################
################################################################################

library(ape)
library(phytools)
library(OUwie)

################################################################################
############################ Read tree and dataset #############################
################################################################################

tree <- read.tree(here("data", "ultrametric_tree.nwk"))
dat  <- read.csv(here("data","i_type_input_data.csv"), stringsAsFactors = FALSE)

################################################################################
############################ Check species matching ############################
################################################################################

setdiff(dat$species, tree$tip.label)   # should return character(0)

################################################################################
############################ Define evolutionary regimes #######################
################################################################################

states <- setNames(dat$regime, dat$species)

# Collapse rare regimes (<5 species)
tab <- table(states)
states[states %in% names(tab)[tab < 5]] <- "Other"

dat$regime <- states[match(dat$species, names(states))]

################################################################################
############################ Stochastic mapping ################################
################################################################################

nsim <- 10
simmap_trees <- make.simmap(tree, states, model="ER", nsim=nsim)

################################################################################
############################ Estimate σ² using OUwie ###########################
################################################################################

ouwie_dat <- data.frame(
  species = dat$species,
  regime  = dat$regime,
  trait   = dat$gene_count,
  stringsAsFactors = FALSE
)

results <- lapply(simmap_trees, function(tr){
  tryCatch(
    OUwie(tr, ouwie_dat, model="OU1", simmap.tree=TRUE, diagn=FALSE),
    error=function(e) NULL
  )
})

# Extract sigma² from the last valid model
valid_models <- results[!sapply(results, is.null)]
sigma_value <- valid_models[[1]]$solution["sigma.sq", ]

print(sigma_value)

################################################################################
############################ Rescale the tree ##################################
################################################################################

simt <- simmap_trees[[1]]

tree_rescaled <- simt 
tree_rescaled$edge.length <- simt$edge.length * sigma_value

################################################################################
############################ Save the rescaled tree ############################
################################################################################

write.tree(tree_rescaled, file="tree_rescaled_i-type_OU.nwk")