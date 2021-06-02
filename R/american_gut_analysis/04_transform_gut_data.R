
## Transform American Gut data to be correct format for GEE
# and connect with metadata 
## Currently it is rows = OTU, sample = Col 

# Load data 
#filtered_data <- read_rds(here::here("Data","American_Gut","genus_threshold_10.rds"))
# try using only 100 samples to make it run... 
#filtered_data <- read_rds(here::here("Data","American_Gut","samples100.rds"))
filtered_data <- read_rds(here::here("Data","American_Gut","antibiotic_otu_table_100.rds"))

## Step 1 - transpose data 
# Create lookup table 
otu_lookup <- filtered_data %>% 
     select(`#OTU ID`,"k","p","c","o","f","g")



# Remove OTU id col 
# Remove k-s cols - not needed now. 
# transposing and keeping row and colnames is annoying...
prep_step <- filtered_data %>% 
     ungroup() %>% 
     select(-c("k","p","c","o","f","g"))
colnames(prep_step)
head(prep_step)

rownames(prep_step) <- prep_step$`#OTU ID`
prep_step <- prep_step %>% select(-`#OTU ID`)

test_t <- t(prep_step)
colnames(test_t) <- otu_lookup$`#OTU ID`
colnames(test_t)
sample_ids <- rownames(test_t)
test_t2 <- as_tibble(test_t) %>% 
    mutate(sample_id = sample_ids)
dim(test_t2)
dim(filtered_data)

# try using dplyr - clean this step up later. 
# test_t <- prep_step %>% 
#     mutate(group = 1) %>%
#     spread(HEADER, price)
#install.packages("janitor")
# library(janitor)
# gut_t_col <- gut_transposed %>% 
#     row_to_names(row_number = 1)
# 
# gut_t_colnames <- gut_t_col %>% 
#     as_tibble() %>% 
#     rownames_to_column()
# rownames(gut_t_col)
# colnames(gut_t_colnames)
# head(gut_t_colnames)



## Load in the metadata to join it with the OTU data
meta <- read_rds(here::here("Data","American_Gut","selected_metadata.rds"))

# Join, remove unnecessary columns 
combined <- test_t2 %>% 
    left_join(meta, by = c("sample_id" = "#SampleID")) %>% 
    select(-BODY_SITE)

#write_rds(combined, here::here("Data","American_Gut","anti.rds"))


## GEE model has 1 row per observation. 1 row per OTU per sample.
## pivot_longer 
gut_otu_response <- combined %>% 
    pivot_longer(cols = !c("sample_id","TYPES_OF_PLANTS","AGE","ANTIBIOTIC_SELECT"), 
                 names_to = "OTU_name", 
                 values_to = "OTU_value")


## convert into presence/absence 
gut_otu_repsonse_presence <- gut_otu_response %>% 
    mutate(presence = OTU_value > 0) %>% 
    mutate(sample_id = factor(sample_id), 
           AGE = as.numeric(AGE), 
           use_antibiotic_past_year = ifelse(ANTIBIOTIC_SELECT == "Not in the last year", 0, 1))



write_rds(gut_otu_repsonse_presence, here::here("Data","temp_data","antibiotic_data_100.rds"))



### create zcor for full matrix 
(sample_size <- nrow(combined))

as.matrix(1:10, nrow = 2)
rep(matrix(c(1,2,3, 11,12,13), nrow = 2, ncol = 3),2)


max_structure_OTU <- max(structure.OTU) - 1
n <- sample_size

zcor_presence = matrix(nrow = n*choose(N , 2), ncol = max_structure_OTU )
for (i in 1:max_structure_OTU ) {
    zcor_presence[, i] = rep(as.numeric(structure.OTU[lower.tri(structure.OTU)] == i + 1),n)
}
dim(zcor_presence)