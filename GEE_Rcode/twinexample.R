library(tidyverse)
library(readxl)
twindata <- read_excel("Data/twindata.xls")



# OTU correlation structure - KNOWN, shown in paper
dim <- list(9,c(4,1,4))


# Gamma matrix - same as in paper 
gamma <- cor2zcor(1,dim,1,corstr="exchangeable",otustr="exchangeable")[[1]] -1
gamma
dim(gamma)
# Omega matrix - for 2 twins and 2 timepoints
omega <- cor2zcor(1,dim =1,c(2,2),corstr="exchangeable",otustr="exchangeable")[[1]]-1
omega
dim(omega)


# R matrix 9x4 is dim 36 - same as paper 
R <- cor2zcor(1,dim,c(2,2),corstr="exchangeable",otustr="exchangeable")
dim(R[[2]])
dim(R[[1]])

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


 
# Mothers were also sampled in the study - remove for this analysis. 
data_twins_only <- data_filter_otus%>% 
    dplyr::filter(ZYGOSITY != "NA") %>% 
    mutate(twin_family = paste(HOST_SUBJECT_ID, FAMILY)) %>% 
    select(-c("HOST_SUBJECT_ID", "FAMILY","ZYGOSITY"))


# Missing twin data 
# Figure out where missing data occurs 
data_twins_only  %>% 
    group_by(HOST_SUBJECT_ID) %>% summarize(n = n()) %>% 
    filter(n != 2)
data_twins_only  %>% 
    group_by(FAMILY) %>% summarize(n = n())


# Introduce NA rows for missing data 
expected_combos <- expand_grid(
    twin_family = unique(data_twins_only$twin_family),
    time = c("TimePoint1","TimePoint2"))

data_include_missing <- left_join(expected_combos, data_twins_only,
                                  by = c("twin_family" = "twin_family",
                                         "time" = "SAMPLEDATE")) %>% 
    separate(twin_family, c("id","family"))



# pivot to be correct data structure, each row is a 
# different observation on timepoints or OTU value 
long_data <- pivot_longer(data_include_missing, 
                          cols = `Blautia`:`Ruminococcus.1`, 
                          names_to = "OTU") 
dim(long_data)
# We have the right dimension! Yay. 


# Need to check if the order this data is in will match 
# the given correlation structure matrix R 





# Presence/absence 
data0 <- long_data %>% 
    mutate(presence = ifelse(value == 0, 0, 1)) %>% 
    mutate(obesity = as.factor(OBESITYCAT),
           family = as.factor(family))
dim(data0)


library(geepack)
# need to also specify id 
# I think 'family' might work as a 'cluster?'
zcor <- R[[2]]
dim(zcor)
geeglm(
    presence ~ obesity,
    data = data0,
    family = "binomial",
    id = factor(family),
    corstr = "userdefined",
    zcor = zcor
)

# How to calculate RA? is it only of these OTUs? do i use the total column? 

# Fit model 