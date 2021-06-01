 library(tidyverse)
 library(geepack)
 
 
# Start with RA model, since is faster 
log_ra <- read_rds(here::here("Data","temp_data","log_ra_50.rds")) %>% 
    mutate(sample_id = factor(sample_id))
zcor_ra <- read_rds(here::here("Data","temp_data","zcor_ra_50.rds"))



start_time <- proc.time()
mod_normal <- geeglm(log_ra ~ antibiotic_year, family = gaussian,
                           data = log_ra, id = sample_id,
                           corstr = "userdefined", zcor = zcor_ra)
end_time <- proc.time()
diff <- end_time - start_time
write_rds(diff, here::here("R","american_gut_analysis","normal_antibiotic50_time.rds"))
write_rds(mod_normal,here::here("R","american_gut_analysis","normal_antibiotic50_mod.rds") )



# Next logsitic model - which is slower. 
presence_data <- read_rds(here::here("Data","temp_data","cleaned_data_50.rds")) 
zcor_logistic <- read_rds(here::here("Data","temp_data","zcor_logistic_50.rds"))

start_time2 <- proc.time()
mod_logistic <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, 
                           corstr = "userdefined", zcor = zcor, id = sample_id2)
end_time2 <- proc.time()
diff2 <- end_time2 - start_time2
write_rds(diff2, here::here("R","american_gut_analysis","logistic_age_time.rds"))
write_rds(mod_logistic,here::here("R","american_gut_analysis","logistic_age_gee_mod.rds") )
