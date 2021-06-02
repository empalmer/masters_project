#filtered_data <- read_rds(here::here("Data","American_Gut","genus_threshold_10.rds"))
#zcor_full <- read_rds(here::here("Data","American_Gut","zcor_full.rds"))

library(tidyverse)

## Load in the "small" data 
## There are 100 samples and 164 different OTUs/taxa
## Will have 100*164 rows 
#filtered_data <- read_rds(here::here("Data","temp_data","prevalence100.rds"))
filtered_data <- read_rds(here::here("Data","temp_data","antibiotic_data_100.rds"))
dim(filtered_data)

# load in the R and zcor matrices for 1 sample 
r_zcor <- read_rds(here::here("Data","temp_data","r_zcor_antibiotic_100.rds"))

R1 <- r_zcor[[1]]
zcor1 <- r_zcor[[2]]

# Calculate the zcor for all samples 
# has dimension 100*choose(164,2) x 55
# number of correlations
(max_R <- max(R1) - 1)
#number of samples
n <- 100
# number of OTUs
N <- nrow(R1)
#create zcor matrix for all samples. 
zcor = matrix(nrow = n*choose(N , 2), ncol = max_R)
for (i in 1:max_R ) {
    zcor[, i] = rep(as.numeric(R1[lower.tri(R1)] == i + 1),n)
}
dim(zcor)




library(geepack)
# If we have a 'userdefined' correlation structure from zcor
# Time and save model output 
start_time <- proc.time()
# Model fit 
geepack_fit_zcor <- geeglm(presence ~ use_antibiotic_past_year, family = binomial, data = filtered_data, 
                           corstr = "userdefined", zcor = zcor, id = sample_id)

diff <- proc.time() - start_time

# Save outputs
write_rds(diff, here::here("R","american_gut_analysis","logistic_antibiotic_time.rds"))
write_rds(geepack_fit_zcor,here::here("R","american_gut_analysis","logistic_antibiotic_gee_mod.rds") )










