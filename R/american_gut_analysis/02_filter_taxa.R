### Examples with the gut data, running throught the pipeline of fxns in the right order 
# Load the data 
library(tidyverse)

#load in data
gut_data <- read_rds(here::here("Data","American_Gut","gut_data_taxonomy.rds"))
taxnames <- c("k","p","c","o","f","g","s")
dim(gut_data)

# Filter out low sampling depth 
## First remove low sampling depth values. 
# sampling depth would be the sum of each column
# total OTUS for each sample

# Takes 1 min to run
sampling_depth <- gut_data %>% 
    select(-c("#OTU ID","k","p","c","o","f","g","s")) %>% 
    summarize_all(sum)

# Decide on cutoff to filter out 
cutoff <- 500
# list of samples to remove
low_sampling_depth_sample_names <- colnames(sampling_depth)[sampling_depth < cutoff]

gut_drop_low_sampling_depth <- gut_data %>% 
    select(-all_of(low_sampling_depth_sample_names))
# drops ~400 samples
dim(gut_drop_low_sampling_depth)



# View distribution of sampling depths. 
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
# load metadata 
selected_metadata <- read_rds(here::here("Data","American_Gut","selected_metadata.rds"))
head(selected_metadata)

# just select one sample site 
# and has a value for antibiotic 
unique(selected_metadata$ANTIBIOTIC_SELECT)
selected_sample_site <- selected_metadata %>% 
    filter(BODY_SITE == "UBERON:feces") %>% 
    filter(ANTIBIOTIC_SELECT != "unknown" & ANTIBIOTIC_SELECT != "no_data") %>% 
    pull(`#SampleID`)
    #filter(ANTIBIOTIC_SELECT == "Not in the last year") %>% 
    

# Subset the columns to be only the samples at the given sample site 
names(gut_bacteria)
filtered_data <- gut_bacteria[, names(gut_bacteria) %in% c("#OTU ID", selected_sample_site, taxnames) ]

# removes about 1000 samples


# save intermediate step 
#write_rds(gut_bacteria_site, here::here("Data","American_Gut","gut_bacteria_site.rds"))
# fill blanks and sum at genus level
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))
# Load data ---- 
#filtered_data <- read_rds(here::here("Data","American_Gut","gut_bacteria_site.rds"))



# Fill in blanks ---- 
gut_full_blank <- blank_taxa_fill(filtered_data)

# sum at genus level  ---- 
gut_full_sum_genus_level <- sum_at_taxa_level(gut_full_blank,c("k","p","c","o","f","g"))
sum_genus_filtered <- gut_full_sum_genus_level

# write rds ---- 
write_rds(gut_full_sum_genus_level , here::here("Data","American_Gut","sum_genus_filtered.rds"))




sum_genus_filtered <- read_rds(here::here("Data","American_Gut","sum_genus_filtered.rds"))
dim(sum_genus_filtered)    
# this filters down a LOT!! 


# genera sparsity - currently done using presence/absense 
# Potential do do it based on abundance, save this option for later 
# Filter out OTU rows where the proportion of 0 counts is <.5 across columns 
start <- colnames(sum_genus_filtered)[8]
end <- colnames(sum_genus_filtered)[ncol(sum_genus_filtered)]


# Calculate the rowwise mean across columns
# For every row, calculate the proportion of columns that are non zero
genera_sparsity <- sum_genus_filtered  %>% 
    rowwise() %>% 
    mutate(per_nonzero = mean(c_across(start:end) > 0))


# Only include an OTU if it is present in more than 5% of samples
genus_thresholded <- genera_sparsity %>% 
    filter(per_nonzero > .1) %>% 
    select(-per_nonzero)
dim(genus_thresholded)

# Now only 165 OTUs included! (for 10% sparcity )


## Save full dataset (ie not 100 randomly sampled )
write_rds(genus_thresholded, here::here("Data","American_Gut","antibiotic_otu_table_all_samples.rds"))  





#### Randomly select 100 samples to use 
head(colnames(genus_thresholded))

sampled_cols <- colnames(genus_thresholded)[c(1:7,
                                            sample(8:length(colnames(genus_thresholded)), size = 100))]
sampled_cols
filtered_data_100_samples <- genus_thresholded[,sampled_cols]
dim(filtered_data_100_samples)


write_rds(filtered_data_100_samples, here::here("Data","American_Gut","antibiotic_otu_table_100.rds"))  

#write_rds(genus_thresholded, here::here("Data","American_Gut","genus_threshold_10.rds"))    
    
    




