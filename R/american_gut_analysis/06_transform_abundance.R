library(tidyverse)
# Load data with 100 samples
#filtered_data <- read_rds(here::here("Data","American_Gut","prevalence100.rds"))
# load in correlation matrix for geeM package
# has dimension 164x164
#r_mat <- read_rds(here::here("Data","American_Gut","R_100.rds"))
# Load in the zcor for geepack package
# has dimension 94*choose(162,2) x 53
#zcor <- read_rds(here::here("Data","American_Gut","zcor_100_samples.rds"))

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


sampling_depth <- read_rds(here::here("Data","American_Gut","sampling_depth.rds"))



# question here: is abundance OTU count/ sampling depth? 
# Is this before/after removing low occurance OTUs?
# For both predictor (RA) and response (AGE) (for now)
# replace everything with NA if the otu is not present in the sample 
# Then calculate RA 
# -log10 transform RA
joined <- filtered_data %>% 
    left_join(sampling_depth) %>% 
    mutate(nonzero_antibiotic = ifelse(OTU_value > 0, use_antibiotic_past_year, NA),
           ra = ifelse(OTU_value > 0, OTU_value/sampling_depth, NA)) %>% 
    mutate(log_ra = -log10(ra), 
           sample_id = factor(sample_id))
dim(joined)

write_rds(joined, here::here("Data","temp_data","log_ra_antibiotic_100.rds"))


# not really normally distributed. Maybe because I removed low-prevalence OTUs earlier?
hist(joined$log_ra)


# have to adjust the zcor now that there will be "missing" data when there are 
# 0 count OTU values
# 1st, make sure to load in the functions below. 

zcor1 <- zcor[1:choose(162,2),]
write_rds(zcor1, here::here("Data","temp_data","zcor_100_1.rds"))

zcor_reduced <- adjust_zcor(joined,162, 92, r_mat, zcor1, "ra")
dim(zcor)
dim(zcor_reduced)    
    
write_rds(zcor_reduced, here::here("Data","American_Gut","zcor_reduced_100.rds"))

    
