### Examples with the gut data, running throught the pipeline of fxns in the right order 
# Load the data 
library(tidyverse)

#load in data
gut_data <- read_rds(here::here("Data","American_Gut","gut_data_taxonomy.rds"))
taxnames <- c("k","p","c","o","f","g","s")

small <- head(gut_data)


# Filter out low sampling depth 
## First remove low sampling depth values. 
# sampling depth would be the sum of each column
# total OTUS for each sample
# testing this workflow
# small_sample <- small[,10:100] %>% 
#     summarize_all(sum)
# length(small_sample > 0)
# colnames(small_sample)[small_sample>0]
# colnames(small_sample > 0)

# Takes 1 min to run
sampling_depth <- gut_data %>% 
    select(-c("#OTU ID","k","p","c","o","f","g","s")) %>% 
    summarize_all(sum)

cutoff <- 100
low_sampling_depth_sample_names <- colnames(sampling_depth)[sampling_depth < cutoff]

gut_drop_low_sampling_depth <- gut_data %>% 
    select(-all_of(low_sampling_depth_sample_names))
# drops ~400 samples
dim(gut_drop_low_sampling_depth)
dim(gut_data)


# View distribution of sampling depths. 
# ASK ASK what is a good cutoff? 100? 
ggplot(data.frame(x = t(sampling_depth))) + 
    geom_histogram(aes( x= x))+
    xlim(0,10000)
hist(t(sampling_depth))

# write if need to save
# write_rds(gut_drop_low_sampling_depth, here::here("Data","American_Gut","sampling_depth.rds"))    



#### Filter to be only kingdom Bacteria OTUs. 
gut_bacteria <- gut_drop_low_sampling_depth %>% 
    filter(k == "k__Bacteria")
dim(gut_bacteria)
# removes ~ 100 rows 




######## Filter to only be samples from one body site 
# metadata 
selected_metadata <- read_rds(here::here("Data","American_Gut","selected_metadata.rds"))
head(selected_metadata)

# just select one sample site 
selected_sample_site <- selected_metadata %>% 
    filter(BODY_SITE == "UBERON:feces") %>% 
    pull(`#SampleID`)

# Subset the columns to be only the samples at the given sample site 
names(gut_bacteria)
gut_bacteria_site <- gut_bacteria[, names(gut_bacteria) %in% c("#OTU ID", selected_sample_site, taxnames) ]
dim(gut_bacteria_site)
# removes about 1000 samples
write_rds(gut_bacteria_site, here::here("Data","American_Gut","gut_bacteria_site.rds"))





# fill blanks and sum at genus level... - done in job_sum.R file 



# takes about 6 minutes to run 
sum_genus_filtered <- read_rds(here::here("Data","American_Gut","sum_genus_filtered.rds"))
dim(sum_genus_filtered)    
# this filters down a LOT!! 




# genera sparsity - currently done using presence/absense 
# Potential do do it based on abundance, save this option for later 
# Filter out OTU rows where the proportion of 0 counts is <.5 across columns 
colnames(sum_genus_filtered)[ncol(sum_genus_filtered)]


# start 9:01
# Calculate the rowwise mean across columns
# For every row, calculate the proportion of columns that are non zero
genera_sparsity <- sum_genus_filtered  %>% 
    rowwise() %>% 
    mutate(per_nonzero = mean(c_across(`000003350.1131718`:`000004625.1131944`) > 0))


# Only include an OTU if it is present in more than 5% of samples
genus_thresholded <- genera_sparsity %>% 
    filter(per_nonzero > .1) %>% 
    select(-per_nonzero)
dim(genus_thresholded)
# Now only 215 OTUs included! 




write_rds(genus_thresholded, here::here("Data","American_Gut","genus_threshold_10.rds"))    
    
    






    
# test to make sure my dplyr code is working
# s1 <- c(0,0,0,0,0,0,0,0)
# s2 <- c(1,0,0,0,0,0,0,0)
# s3 <- c(1,1,1,1,0,0,0,0)
# 
# test_colsum <- data.frame(s1,s2,s3)
# test_colsum %>% rowwise() %>%
#     mutate(per = mean(c_across(s1:s3 ) > 0))

