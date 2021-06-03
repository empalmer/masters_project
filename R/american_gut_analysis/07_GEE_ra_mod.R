# load in functions 
library(tidyverse)
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))


# Load in data 
# -log10 transformed relative abundance
log_ra <- read_rds(here::here("Data","temp_data","log_ra_antibiotic_100.rds"))

# R and zcor for 1 sample:
r_zcor <- read_rds(here::here("Data","temp_data","r_zcor_antibiotic_100.rds"))
R1 <- r_zcor[[1]]
zcor1 <- r_zcor[[2]]


# calculate the zcor_reduced 
# 164 OTUs, #100 samples 
# This takes ~ 1 min to run. 
zcor_reduced <- adjust_zcor(log_ra, 164, 100, R1, zcor1, "ra")


# Model fit using geepack 
library(geepack)

# If we have a 'userdefined' correlation structure from zcor
# takes MANY HOURS 
start_time <- proc.time()
geepack_fit_zcor <- geeglm(log_ra ~ use_antibiotic_past_year, family = gaussian,
                           data = log_ra, id = sample_id, waves = OTU_name,
                           corstr = "userdefined", zcor = zcor_reduced)

diff <- proc.time() - start_time
write_rds(diff, here::here("R","american_gut_analysis","normal_antibiotic_time2.rds"))
write_rds(geepack_fit_zcor,here::here("R","american_gut_analysis","normal_antibiotic_mod_wav.rds") )
# Graphical display 



############################################
# Previous attempts with other covariates: 
# try with plants as covariates 
# For now, collapse into two groups so the model will run overnight. 
# 0 if < 10 plants 1 if > 10 plants
# Filter and then adjust zcor to only include samples with plant data
# Now there are only 91 samples (vs 94)

# log_ra_plants <- log_ra %>% 
#     filter(TYPES_OF_PLANTS != "unknown" & TYPES_OF_PLANTS != "no_data") %>% 
#     mutate(plant_ind = ifelse(TYPES_OF_PLANTS %in% c("Less than 5","6 to 10"), 0 , 1) )
# 
# # adjust zcor: 
# zcor_plants <- adjust_zcor(log_ra_plants, 162, 91, r_mat, zcor1, "ra")
# 
# 
# # run model (overnight) takes ~10 hours to complete
# start_time <- proc.time()
# geepack_fit_zcor_plant2cat <- geeglm(log_ra ~ age_present + plant_ind, family = gaussian,
#                                      data = log_ra_plants, id = sample_id2,
#                                      corstr = "userdefined", zcor = zcor_plants)
# end_time <- proc.time()
# diff <- end_time - start_time
# write_rds(diff, here::here("R","american_gut_analysis","normal_age_time_plant2cat.rds"))
# write_rds(geepack_fit_zcor_plant2cat,here::here("R","american_gut_analysis","normal_age_gee_mod_plant2cat.rds") )
# 
# summary(geepack_fit_zcor_plant2cat)














# Check correct dimension of zcor
clusz <- log_ra %>% 
     group_by(sample_id) %>%
     summarize(n = sum(!is.na(ra)))  %>% pull(n)
sum(clusz * (clusz - 1) / 2)
nrow(zcor_reduced)

