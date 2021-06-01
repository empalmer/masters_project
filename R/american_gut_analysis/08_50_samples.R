### Examples with the gut data, running throught the pipeline of fxns in the right order 
# Load the data 
library(tidyverse)

#load in data
gut_data <- read_rds(here::here("Data","American_Gut","gut_data_taxonomy.rds"))
taxnames <- c("k","p","c","o","f","g","s")


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

rownames_sampling_depth <- rownames(t(sampling_depth))
sampling_depth50 <- as_tibble(t(sampling_depth)) %>% 
    mutate(sample_id = rownames_sampling_depth, 
           sampling_depth = V1) %>% 
    select(-V1)

write_rds(sampling_depth50, here::here("Data","American_Gut","sampling_depth50.rds")) 


cutoff <- 500
low_sampling_depth_sample_names <- colnames(sampling_depth)[sampling_depth < cutoff]

gut_drop_low_sampling_depth <- gut_data %>% 
    select(-all_of(low_sampling_depth_sample_names))
# drops ~400 samples
dim(gut_drop_low_sampling_depth)



# View distribution of sampling depths. 
# What is a good cutoff? 100? 
ggplot(data.frame(x = t(sampling_depth))) + 
    geom_histogram(aes( x= x))+
    xlim(0,10000) + 
    geom_vline(xintercept = 500)
hist(t(sampling_depth))

# write if need to save
#write_rds(gut_drop_low_sampling_depth, here::here("Data","American_Gut","sampling_depth50.rds"))    



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
selected_sample_site_ids <- selected_metadata %>% 
    filter(BODY_SITE == "UBERON:feces") %>% 
    pull(`#SampleID`)
gut_bacteria_site <- gut_bacteria[, names(gut_bacteria) %in% c("#OTU ID", selected_sample_site_ids, taxnames) ]
# This is the dimension of the original dataset before other filtering - 
unfiltered_dim <- dim(gut_bacteria_site)
################################
# Now randomly sample only 50 samples, even 100 wont run 
sample_ids50 <- selected_sample_site_ids[sample(1:length(selected_sample_site_ids), 50, replace = F)]

# Subset the columns to be only the samples at the given sample site 
gut_50 <- gut_bacteria[, names(gut_bacteria) %in% c("#OTU ID", sample_ids50, taxnames) ]



# fill blanks and sum at genus level for 50 samples
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))

# Fill in blanks ---- 
gut_full_blank <- blank_taxa_fill(gut_50)

# sum at genus level  ---- 

#sum_at_taxa_level(gut_full_blank, c("k","p","c","o","f","g"))

sum_genus_filtered <-  gut_full_blank %>% 
    group_by(k,p,c,o,f,g) %>% 
    summarize_if(is.numeric, sum)


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

names(genus_thresholded)
# Now only 160 OTUs included! (for 10% sparcity )

write_rds(genus_thresholded, here::here("Data","temp_data","samples_50.rds"))    
samples_50 <- genus_thresholded
dim(samples_50)


### Create R matrix 
taxa_count_list <- tax2cor(samples_50, taxa_levels = c("k","p","c","o","f","g"))
taxa_count_list

#extract just the counts - the dimension list to give to the cor2zcor fxn. 
dim <- map(taxa_count_list,2)
dim


taxa_list <- list("k",c("k","p"),c("k","p","c"),c("k","p","c","o"),
                  c("k","p","c","o","f"),c("k","p","c","o","f","g"))

all <- map(taxa_list, ~taxa_group_tally(samples_50,.x))
dim <- map(all, "n")
dim

source(here::here("GEE_Chen2020","cor2zcor.r"))

R_zcor <- cor2zcor(n = 1, dim = dim, par = 1)

write_rds(R_zcor, here::here("Data","temp_data","r_zcor_50.rds")) 


R1 <- R_zcor[[1]]
zcor1 <- R_zcor[[2]]


(max_cor <- max(R1) -1)



## Transform data to correct structure: 
prep_step <- samples_50 %>% 
    ungroup() %>% 
    select(-c("k","p","c","o","f","g"))

otu_lookup <- samples_50  %>% 
    select(`#OTU ID`,"k","p","c","o","f","g")

rownames(prep_step) <- prep_step$`#OTU ID`
prep_step <- prep_step %>% select(-`#OTU ID`)

test_t <- t(prep_step)
colnames(test_t) <- otu_lookup$`#OTU ID`
colnames(test_t) <- otu_lookup$g
colnames(test_t)
sample_ids <- rownames(test_t)
test_t2 <- as_tibble(test_t) %>% 
    mutate(sample_id = sample_ids)
names(test_t2)
dim(test_t2)
dim(samples_50)


## Load in the metadata to join it with the OTU data
meta <- read_rds(here::here("Data","American_Gut","selected_metadata.rds"))

# Join, remove unnecessary columns 
combined <- test_t2 %>% 
    left_join(meta, by = c("sample_id" = "#SampleID")) %>% 
    select(-BODY_SITE)

filter_antibiotic  %>% 
    group_by(ANTIBIOTIC_SELECT) %>% 
    tally()
filter_antibiotic <- combined %>% 
    filter(ANTIBIOTIC_SELECT != "no_data" & ANTIBIOTIC_SELECT != "unknown")

dim(filter_antibiotic)
# 43 samples to be used 
names(filter_antibiotic)

gut_long50 <- filter_antibiotic %>% 
    pivot_longer(cols = !c("sample_id","TYPES_OF_PLANTS","AGE","ANTIBIOTIC_SELECT"), 
                 names_to = "OTU_name", 
                 values_to = "OTU_value")


## convert into presence/absence 
cleaned_data <- gut_long50 %>% 
    mutate(presence = OTU_value > 0, 
           AGE = as.numeric(AGE),
           antibiotic_year = ifelse(ANTIBIOTIC_SELECT == "Not in the last year", 1, 0))


write_rds(cleaned_data, here::here("Data","temp_data","cleaned_data_50.rds")) 



log_ra_50 <- cleaned_data %>% 
    left_join(sampling_depth50) %>% 
    mutate(antibiotic_year = ifelse(OTU_value > 0, antibiotic_year, NA),
           ra = ifelse(OTU_value > 0, OTU_value/sampling_depth, NA)) %>% 
    mutate(log_ra = -log10(ra),
           sample_id = factor(sample_id))

write_rds(log_ra_50, here::here("Data","temp_data","log_ra_50.rds"))

# full zcor for logistic model 
max_R <- max(R1) - 1
# sample size 
n <- 43
N <- nrow(R1)
zcor_logistic = matrix(nrow = n*choose(N , 2), ncol = max_R )
for (i in 1:max_R ) {
    zcor_logistic[, i] = rep(as.numeric(R1[lower.tri(R1)] == i + 1),n)
}
dim(zcor_logistic)

write_rds(zcor_logistic, here::here("Data","temp_data","zcor_logistic_50.rds"))



# full zcor for logRA model 
rzcor <- read_rds(here::here("Data","temp_data","r_zcor_50.rds")) 
zcor1 <- rzcor[[2]]
R1 <- rzcor[[1]]
dim(zcor1)
nsamples <- 43
nOTU <- nrow(R1)
zcor_reduced <- adjust_zcor(log_ra_50, nOTU, nsamples, R1, zcor1, "ra")
dim(zcor_reduced)
write_rds(zcor_reduced, here::here("Data","temp_data","zcor_ra_50.rds"))
