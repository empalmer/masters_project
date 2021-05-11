# Load data and packages ------ 
library(tidyverse)
source(here::here("R","american_gut_analysis","taxonomy_to_correlation_list.R"))
source(here::here("GEE_Chen2020","cor2zcor.r"))

# Load data ---- 
filtered_data <- read_rds(here::here("Data","American_Gut","gut_bacteria_site.rds"))

# Fill in blanks ---- 
gut_full_blank <- blank_taxa_fill(filtered_data)

# sum at genus level  ---- 
gut_full_sum_genus_level <- sum_at_taxa_level(gut_full_blank,c("k","p","c","o","f","g"))


# write rds ---- 
write_rds(gut_full_sum_genus_level , here::here("Data","American_Gut","sum_genus_filtered.rds"))