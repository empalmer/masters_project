
## Read in filtered data 
# This data is: 
# Samples filtered to have sampling depth > 100
# OTUs sampled to only be kingdom Bacteria 
# Samples filtered to only be fecal samples 
# OTUs filtered to be present in > 10% of samples 
# Summed to the genus level 

filtered_data <- read_rds(here::here("Data","American_Gut","genus_threshold_10.rds"))


## get the taxa list for calculating the correlation matrix. 
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))
# get the dimensions 
taxa_count_list <- tax2cor(filtered_data, taxa_levels = c("k","p","c","o","f","g"))
taxa_count_list


#extract just the counts - the dimension list to give to the cor2zcor fxn. 
dim <- map(taxa_count_list,2)
dim


# Calculate the R matrix. 
# Currently, see 03_manual_R_calculation