################################################################################
############################ Load required packages ############################
################################################################################

library(ape)
library(nlme)
library(phytools)
library(here)

################################################################################
############################ Read tree and datasets ############################
################################################################################

pgls_tree_c_type <- read.tree(here("PGLS_analysis","data","tree_rescaled_c-type_BMS.nwk")) 
pgls_tree_i_type <- read.tree(here("PGLS_analysis","data","tree_rescaled_i-type_OU.nwk"))

c_data_larvae <- read.csv(here("PGLS_analysis","data","pgls_c_type_input_data_larvae.csv"), row.names = 1)
i_data_larvae <- read.csv(here("PGLS_analysis","data","pgls_i_type_input_data_larvae.csv"), row.names = 1)
c_data_adult <- read.csv(here("PGLS_analysis","data","pgls_c_type_input_data_adults.csv"), row.names = 1)
i_data_adult <- read.csv(here("PGLS_analysis","data","pgls_i_type_input_data_adults.csv"), row.names = 1)

################################################################################
############################ Prepare variables #################################
################################################################################

prepare_data <- function(df){

  df$food_source <- relevel(as.factor(df$food_source), ref="Plants")

  # log-transform to avoid log(0)
  df$log_gene <- log(df$gene_count + 1)

  # species column required by corBrownian
  df$species <- rownames(df)

  return(df)
}

c_data_larvae <- prepare_data(c_data_larvae)
i_data_larvae <- prepare_data(i_data_larvae)
c_data_adult <- prepare_data(c_data_adult)
i_data_adult <- prepare_data(i_data_adult)


################################################################################
############################ Brownian correlation ##############################
################################################################################


corBM_rescaled_c_type <- corBrownian(1, phy=pgls_tree_c_type, form = ~ species)
corBM_rescaled_i_type <- corBrownian(1, phy=pgls_tree_i_type, form = ~ species)

################################################################################
############################ Fit PGLS models ###################################
################################################################################

c_fit_larvae <- gls(log_gene ~ food_source, data=c_data_larvae,
             correlation = corBM_rescaled_c_type, method="REML")

i_fit_larvae <- gls(log_gene ~ food_source, data=i_data_larvae,
             correlation = corBM_rescaled_i_type, method="REML")
c_fit_adult <- gls(log_gene ~ food_source, data=c_data_adult,
             correlation = corBM_rescaled_c_type, method="REML")

i_fit_adult <- gls(log_gene ~ food_source, data=i_data_adult,
             correlation = corBM_rescaled_i_type, method="REML")

summary(c_fit_larvae)
summary(i_fit_larvae)
summary(c_fit_adult)
summary(i_fit_adult)
