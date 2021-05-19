# load in functions 
library(tidyverse)
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))


# Load in data 
# -log10 transformed relative abundance
log_ra <- read_rds(here::here("Data","temp_data","log_ra_100.rds"))
# zcor for 1 observation
zcor1 <-  read_rds(here::here("Data","temp_data","zcor_100_1.rds"))
# correlation matrix 
r_mat <- read_rds(here::here("Data","temp_data","R_100.rds"))

#calculate the zcor_reduced 
zcor_reduced <- adjust_zcor(log_ra, 162, 94, r_mat, zcor1, "ra")

# convert sample id to factor for geepack
log_ra <- log_ra %>% 
    mutate(sample_id2 = factor(sample_id))

# Model fit using geepack 
library(geepack)
geepack_fit_1 <- geeglm(log_ra ~ age_present, family = gaussian, data = log_ra, id = sample_id2, 
                        corstr = "independence")


# If we have a 'userdefined' correlation structure from zcor
# At least, it will not run for me
geepack_fit_zcor <- geeglm(log_ra ~ age_present, family = gaussian, data = log_ra, id = sample_id2,
                           corstr = "userdefined", zcor = zcor_reduced)



# Check correct dimension of zcor
clusz <- log_ra %>% 
    group_by(sample_id) %>%
    summarize(n = sum(!is.na(ra)))  %>% pull(n)
sum(clusz * (clusz - 1) / 2)
nrow(zcor_reduced)


