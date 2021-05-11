
library(tidyverse)


file <- here::here("Data","American_Gut","AG_100nt_taxonomy.txt")
colnames = c("OTU_id", 1:4828,"k","p","c","o","f","g","s")
taxnames <- c("k","p","c","o","f","g","s")

# first 50 rows to make testing and computation quicker
gut_subset_tax <- read_tsv(file, skip = 1, n_max = 50)

# full dataset
gut_full_tax <- read_tsv(file, skip = 1)

# Both the above have a "taxonomy" column, separate into multiple 
gut_subset <- gut_subset_tax %>% 
    separate(col = "taxonomy",into = taxnames, sep = "; ")

gut_full <- gut_full_tax %>% 
    separate(col = "taxonomy",into = taxnames, sep = "; ")
# as factor OTU id col 



write_rds(gut_full, here::here("Data","American_Gut","gut_data_taxonomy.rds"))





###### Metadata stuff 
file_meta <- here::here("Data","American_Gut","AG_100nt.txt")
meta <- read_tsv(file_meta)


# change the selected once we have good covariates to use 
meta_select <- meta %>% 
    select(`#SampleID`, BODY_SITE, CARBOHYDRATE_PER, AGE)


write_rds(meta_select,here::here("Data","American_Gut","selected_metadata.rds"))
