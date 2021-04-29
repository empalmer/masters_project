 
# This is mostly done in the models.R file.... 

## Get twin data

source(here::here("R","GEE","clean_twindata.R"))
clean_twindata <- clean_twinexample_data()

### Option A: 
# Calculate RA as OTU count / total_read column


ra_A <- clean_twindata %>% 
    mutate(ra = value/total_read) %>% 
    mutate(logt = -log10(ra))

# for two part models - replace 0 counts with NA so wont be included in model
remove_0abund <- ra_A %>% 
    mutate(obesity = ifelse(ra == 0, NA, obesity), 
           logt = ifelse(ra == 0, NA, logt))

# for one part models - replace with 5 so -log10(0) isnt infinity
# unclear from paper which value was used here... 
replace_0abund <- ra_A %>% 
    mutate(logt = ifelse(ra == 0, 5, logt))


clusz_abund <- get_clusz(remove_0abund, family)
zcor_abund <- adjust_zcor(remove_0abund, max_size = 36, dim = list(9,c(4,1,4)), par = c(2,2), n = 54)
nrow(zcor_abund)
sum(clusz_abund * (clusz_abund - 1) / 2)









### Option B: 
# Calculate RA as OTU count / total OTUs included (14)

### Option C: 
# Calculate RA as OTU count / sum OTUs used in paper (9)