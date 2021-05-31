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
    mutate(sample_id2 = factor(sample_id),
           age_present = as.numeric(age_present))

# Model fit using geepack 
library(geepack)
geepack_fit_1 <- geeglm(log_ra ~ age_present, family = gaussian, data = log_ra, id = sample_id2, 
                        corstr = "independence")


# If we have a 'userdefined' correlation structure from zcor
# takes SEVERAL HOURS 
start_time <- proc.time()
geepack_fit_zcor <- geeglm(log_ra ~ age_present, family = gaussian, data = log_ra, id = sample_id2,
                           corstr = "userdefined", zcor = zcor_reduced)
end_time <- proc.time()
diff <- end_time - start_time
write_rds(diff, here::here("R","american_gut_analysis","normal_age_time.rds"))
write_rds(geepack_fit_zcor,here::here("R","american_gut_analysis","normal_age_gee_mod.rds") )




############################################
# try with plants as covariates 
# For now, collapse into two groups so the model will run overnight. 
# 0 if < 10 plants 1 if > 10 plants
# Filter and then adjust zcor to only include samples with plant data
# Now there are only 91 samples (vs 94)

log_ra_plants <- log_ra %>% 
    filter(TYPES_OF_PLANTS != "unknown" & TYPES_OF_PLANTS != "no_data") %>% 
    mutate(plant_ind = ifelse(TYPES_OF_PLANTS %in% c("Less than 5","6 to 10"), 0 , 1) )

# adjust zcor: 
zcor_plants <- adjust_zcor(log_ra_plants, 162, 91, r_mat, zcor1, "ra")


# run model (overnight) takes ~10 hours to complete
start_time <- proc.time()
geepack_fit_zcor_plant2cat <- geeglm(log_ra ~ age_present + plant_ind, family = gaussian,
                                     data = log_ra_plants, id = sample_id2,
                                     corstr = "userdefined", zcor = zcor_plants)
end_time <- proc.time()
diff <- end_time - start_time
write_rds(diff, here::here("R","american_gut_analysis","normal_age_time_plant2cat.rds"))
write_rds(geepack_fit_zcor_plant2cat,here::here("R","american_gut_analysis","normal_age_gee_mod_plant2cat.rds") )

summary(geepack_fit_zcor_plant2cat)














# Check correct dimension of zcor
clusz <- log_ra %>% 
    group_by(sample_id) %>%
    summarize(n = sum(!is.na(ra)))  %>% pull(n)
sum(clusz * (clusz - 1) / 2)
nrow(zcor_reduced)


