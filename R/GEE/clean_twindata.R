library(tidyverse)
library(readxl)
twindata <- read_excel("Data/twindata.xls")


# returns cleaned and formated twin data for this specific example
clean_twinexample_data <- function(){

## Data dimension check
# Y is OTU 
# Y should have dim K x J where 
# K <- 54 - samples/clusters 108 individuals, 54 twin pairs
# J <- 36 - measurements - 2 time points, 2 twins 9 otus 
# 2 * 2 * 9 = 36
#  K*J = 54*36 = 1944 rows in data 



# Subset the data to only include OTUs used in paper 
# In order Clostridiales,
# Create vector of the 9 OTUs 
otus_used <- c("Blautia","Coprococcus","Lachnospira",
               "Roseburia","Eubacterium","Ruminococcus",
               "Faecalibacterium","Oscillospira","Ruminococcus.1")
meta_vars <- c("HOST_SUBJECT_ID", "SAMPLEDATE", "ZYGOSITY",
               "FAMILY", "OBESITYCAT")


data_filter_otus <- twindata %>% 
    select(all_of(c(meta_vars,otus_used))) 
colnames(data_filter_otus)
dim(data_filter_otus)



# Also remove 1 observation where there is not a measured obesity value 
twins_no_na <- data_filter_otus %>% 
    filter(!is.na(OBESITYCAT))


# Mothers were also sampled in the study - remove for this analysis. 
# remove zygosity column as not important for later work 
# Create an ID column that includes the subject ID and the family
# arranging by family is VERY IMPORTANT
data_twins_only <- twins_no_na %>% 
    dplyr::filter(ZYGOSITY != "NA") %>% 
    mutate(twin_family = paste(HOST_SUBJECT_ID, FAMILY)) %>% 
    select(-c("ZYGOSITY")) %>% 
    arrange(FAMILY)



# Missing twin data 
# Figure out where missing data occurs 
data_twins_only  %>% 
    group_by(HOST_SUBJECT_ID) %>% summarize(n = n()) %>% 
    filter(n != 2)
data_twins_only  %>% 
    group_by(FAMILY) %>% summarize(n = n()) 




# Introduce NA rows for missing data 
# This is done so we have the same cluster size and rows for every observation
# even if it is missing
# This is done to make the final data have the dimensinos the paper says it will

expected_combos <- expand_grid(
    twin_family = unique(data_twins_only$twin_family),
    time = c("TimePoint1","TimePoint2")) %>% 
    mutate(twin_id = rep(c("A","A","B","B"),54))
data_include_missing <- left_join(expected_combos, data_twins_only,
                                  by = c("twin_family" = "twin_family",
                                         "time" = "SAMPLEDATE")) %>% 
    separate(twin_family, c("id","family"))




# pivot to be correct data structure, each row is a 
# different observation on timepoints or OTU value 
# pivot where we dont introduce NA values
long_data_no_na <- pivot_longer(data_twins_only, 
                                cols = `Blautia`:`Ruminococcus.1`, 
                                names_to = "OTU") %>% 
    separate(twin_family, c("id","family"))%>% 
    mutate(OTU = factor(OTU, levels = otus_used))
dim(long_data_no_na)
## get a vector representing how many observations in each cluster 
long_data_no_na %>% 
    group_by(family) %>% tally() %>% 
    pull(n)


# pivot to be correct data structure, each row is a 
# different observation on timepoints or OTU value 
# and where we do so it is the "right dimension 
long_data<- pivot_longer(data_include_missing, 
                         cols = `Blautia`:`Ruminococcus.1`, 
                         names_to = "OTU") %>% 
    mutate(OTU = factor(OTU, levels = otus_used))
dim(long_data)
# This is the same dimension as in the paper
# see we have 36 observations per cluster 4 timepoints and 9 otus
long_data %>% 
    group_by(family) %>% tally() %>% 
    pull(n)


# Change obesity into a Yes/ No status 

# Change obesity column to yes/no 
long_data_neat <- long_data %>% 
    mutate(obesity = ifelse(OBESITYCAT == "Lean",0,1),
           family = as.factor(family)) %>% 
    rename(subject_id = HOST_SUBJECT_ID) %>% 
    select(-c(OBESITYCAT,FAMILY))

long_data_neat %>% 
    group_by(family) %>% tally() %>% 
    pull(n)


## sort to be by OTU, 
long_data_neat <- long_data_neat %>% 
    arrange(family, OTU)

return(long_data_neat)
}