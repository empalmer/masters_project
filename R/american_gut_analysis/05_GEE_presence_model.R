#filtered_data <- read_rds(here::here("Data","American_Gut","genus_threshold_10.rds"))
#zcor_full <- read_rds(here::here("Data","American_Gut","zcor_full.rds"))

library(tidyverse)

## Load in the "small" data 
## There are 94 samples and 162 different OTUs/taxa
## Will have 94*162 rows 
filtered_data <- read_rds(here::here("Data","temp_data","prevalence100.rds"))

# load in correlation matrix for geeM package
# has dimension 162x162
r_mat <- read_rds(here::here("Data","temp_data","R_100.rds"))
# Load in the zcor for geepack package
# has dimension 94*choose(162,2) x 54
# too large to store, calculate below
#zcor <- read_rds(here::here("Data","temp_data","zcor_100_samples.rds"))


max_R <- max(r_mat) - 1
n <- 94
N <- nrow(r_mat)
zcor = matrix(nrow = n*choose(N , 2), ncol = max_R )
for (i in 1:max_R ) {
    zcor[, i] = rep(as.numeric(r_mat[lower.tri(r_mat)] == i + 1),n)
}
dim(zcor)





# Some transformations to make the code work 
# geepack needs sample_id as a factor, otherwise it tries to convert to 
# numeric and runs into problems. 
# Age should be treated as numeric, as it is originally a character
filtered_data <- filtered_data %>% 
    mutate(sample_id2 = factor(sample_id), 
           AGE = as.numeric(AGE), 
           presence0 = ifelse(presence,1,0))


library(geepack)
# Try the other "simpler" correlation structures built in to geepack
# These all will run in a reasonable amount of time. (< 1m) 
#geepack_fit_1 <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, id = sample_id2, 
#                      corstr = "independence")
#geepack_fit_2 <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, id = sample_id2, 
#                        corstr = "exchangeable")
#geepack_fit_3 <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, id = sample_id2, 
                        corstr = "ar1")
#summary(geepack_fit_1)

# Try the unstructured correlation structure, which should be more complicated than userdefined 
#geepack_fit_4 <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, id = sample_id2, 
#                        corstr = "unstructured")
# If we have a 'userdefined' correlation structure from zcor
# At least, it will not run for me
#geepack_fit_zcor <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, 
#                      corstr = "userdefined", zcor = zcor, id = sample_id2)


# If we have a 'userdefined' correlation structure from zcor
# Time and save model output 
start_time <- proc.time()
geepack_fit_zcor <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, 
                           corstr = "userdefined", zcor = zcor, id = sample_id2)
end_time <- proc.time()
diff <- end_time - start_time
write_rds(diff, here::here("R","american_gut_analysis","logistic_age_time.rds"))
write_rds(geepack_fit_zcor,here::here("R","american_gut_analysis","logistic_age_gee_mod.rds") )













# Next try the geeM package, which uses the original correlation (sqaure)
# matrix for userdefined 
# Seems a LOT slower. Won't even run with independence structure...
# (or takes > 5 min to run, I've always stopped it after 5 min)
#library(geeM)

#geem_fit <- geem(presence ~ AGE, data = filtered_data, family = binomial, id = sample_id,     corstr = "independence")

# Since it takes so long for independence, I havent tried the userdefined option
#geem(presence ~ AGE, data = filtered_data, family = binomial, id = sample_id,
#     corstr = "userdefined", corr.mat = r_mat)


