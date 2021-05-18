
# Load data and packages ---- 
library(tidyverse)
source(here::here("R","american_gut_analysis","taxonomy_to_correlation_list.R"))
source(here::here("GEE_Chen2020","cor2zcor.r"))



# Try out with full data ----- 
#gut_full_sum_species_level <- read_rds(here::here("Data","American_Gut","sum_species_level.rds"))
# get the dimensions 
#taxa_count_list <- tax2cor(gut_full_sum_species_level, taxa_levels = c("k","p","c","o","f","g","s"))
#taxa_count_list


#extract just the counts - the dimension list to give to the cor2zcor fxn. 
#dim <- map(taxa_count_list,2)
#dim


# run the R matrix 
#R <- cor2zcor(n = 1, dim = dim, par = 1)[[1]]
#R

#write_rds(R, here::here("R","ame"))

# Fails after 2 hours with 
# Error: vector memory exhausted (limit reached?)
# Execution halted


# filtered data ---- 

## load in ----- 
filtered_data <- read_rds(here::here("Data","American_Gut","sum_genus_filtered.rds"))

## get taxa list ---- 
taxa_count_list <- tax2cor(filtered_data, taxa_levels = c("k","p","c","o","f","g"))
taxa_count_list


#extract just the counts - the dimension list to give to the cor2zcor fxn. 
dim <- map(taxa_count_list, 2)
dim


# run the R matrix ---- 
R <- cor2zcor(n = 1, dim = dim, par = 1)[[1]]
R

write_rds(R, here::here("R","american_gut_analysis","R_mat_filter.rds"))
