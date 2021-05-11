

## Transform American Gut data to be correct format for GEE 
## Currently it is rows = OTU, sample = Col 
## Step 1 - transpose data 
# Remove k-s cols 
# Create new lookup table 
otu_lookup <- gut_full_sum_species_level %>% 
    select(`#OTU ID`,"k","p","c","o","f","g","s")

# Remove OTU id col 
prep_step <- gut_full_sum_species_level %>% 
    ungroup() %>% 
    select(-c("k","p","c","o","f","g","s"))
colnames(prep_step)


gut_transposed <- t(prep_step)
install.packages("janitor")
library(janitor)
gut_t_col <- gut_transposed %>% 
    row_to_names(row_number = 1)

gut_t_colnames <- gut_t_col %>% 
    as_tibble() %>% 
    rownames_to_column(var = "subject_id")
colnames(gut_t_colnames)


## insert metadata files here.... 

file_meta <- here::here("Data","American_Gut","AG_100nt.txt")
meta <- read_tsv(file_meta, n_max = 50)


# change the selected once we have good covariates to use 
meta_select <- meta %>% 
    select(`#SampleID`, BODY_SITE, CARBOHYDRATE_PER, AGE)



## pivot_longer 
gut_otu_response <- gut_t_colnames %>% 
    pivot_longer(cols = !subject_id , 
                 names_to = "OTU_name", 
                 values_to = "OTU_value")
