
# Load the data 
library(tidyverse)

#load in data
gut_data <- read_rds(here::here("Data","American_Gut","gut_data_taxonomy.rds"))
taxnames <- c("k","p","c","o","f","g","s")
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
#### Filter to be only kingdom Bacteria OTUs. 
gut_bacteria <- gut_drop_low_sampling_depth %>% 
    filter(k == "k__Bacteria")
dim(gut_bacteria)
# removes ~ 100 rows 

######## Filter to only be samples from one body site 
#### AND no antibiotics in the past year. 
# metadata 
selected_metadata <- read_rds(here::here("Data","American_Gut","selected_metadata.rds"))
head(selected_metadata)

# just select one sample site 
selected_sample_site <- selected_metadata %>% 
    filter(BODY_SITE == "UBERON:feces") %>% 
    filter(ANTIBIOTIC_SELECT == "Not in the last year") %>% 
    pull(`#SampleID`)

selected_sample_site <- selected_sample_site[1:100]

# Subset the columns to be only the samples at the given sample site 
names(gut_bacteria)
gut_bacteria_site <- gut_bacteria[, names(gut_bacteria) %in% c("#OTU ID", selected_sample_site, taxnames) ]
dim(gut_bacteria_site)
colnames(gut_bacteria_site)
# removes about 1000 samples


filtered_data <- gut_bacteria_site

# fill blanks and sum at genus level.
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))


# Fill in blanks ---- 
gut_full_blank <- blank_taxa_fill(filtered_data)

# sum at genus level  ---- 
sum_genus_filtered <- sum_at_taxa_level(gut_full_blank,c("k","p","c","o","f","g"))
dim(sum_genus_filtered)




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


write_rds(genus_thresholded, here::here("Data","American_Gut","samples100.rds"))  








