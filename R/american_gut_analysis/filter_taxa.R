### Examples with the gut data, running throught the pipeline of fxns in the right order 
# Load the data 
library(tidyverse)


gut_data <- read_rds(here::here("Data","American_Gut","gut_data_taxonomy.rds"))

taxnames <- c("k","p","c","o","f","g","s")


gut_bacteria <- gut_data %>% 
    filter(k == "k__Bacteria")
colnames(gut_bacteria)

# metadata 
selected_metadata <- read_rds(here::here("Data","American_Gut","selected_metadata.rds"))
head(selected_metadata)

# just select one sample site 
selected_sample_site <- selected_metadata %>% 
    filter(BODY_SITE == "UBERON:feces") %>% 
    pull(`#SampleID`)


# Subset the columns to be only the samples at the given sample site 
gut_bacteria_site <- gut_bacteria %>% 
    select_at(c("#OTU ID", selected_sample_site, taxnames))


write_rds(gut_bacteria_site, here::here("Data","American_Gut","gut_bacteria_site.rds"))

# fill blanks sum at order level... - done in job_sum.R file 
# takes about 6 minutes to run 
sum_genus_filtered <- read_rds(here::here("Data","American_Gut","sum_genus_filtered.rds"))
    
