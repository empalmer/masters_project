library(tidyverse)
library(readxl)
twindata <- read_excel("Data/twindata.xls")
# Load cor2zcor function 
source(here::here("GEE_Chen2020","cor2zcor.r"))


# What is total read column? 
colnames(twindata)
rowSums(twindata[,10:23])


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
# n is 54 since we have 54 families or clusters 
R <- cor2zcor(54,dim,c(2,2),corstr="exchangeable",otustr="exchangeable")
dim(R[[2]])
#34020
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


# Also remove 1 observation where there is not a measured obesity value 
twins_no_na <- data_filter_otus %>% 
    filter(!is.na(OBESITYCAT))


# Mothers were also sampled in the study - remove for this analysis. 
data_twins_only <- twins_no_na%>% 
    dplyr::filter(ZYGOSITY != "NA") %>% 
    mutate(twin_family = paste(HOST_SUBJECT_ID, FAMILY)) %>% 
    select(-c("ZYGOSITY"))




# Missing twin data 
# Figure out where missing data occurs 
data_twins_only  %>% 
    group_by(HOST_SUBJECT_ID) %>% summarize(n = n()) %>% 
    filter(n != 2)
data_twins_only  %>% 
    group_by(FAMILY) %>% summarize(n = n()) 




# Introduce NA rows for missing data 
# This is done to make the final data have the dimensinos the paper says it will
# but I am not sure it will work for the geeglm() function...
expected_combos <- expand_grid(
    twin_family = unique(data_twins_only$twin_family),
    time = c("TimePoint1","TimePoint2"))
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
    separate(twin_family, c("id","family"))
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
                                names_to = "OTU") 
dim(long_data)
# This is the same dimension as in the paper



# Presence/absence  data 
data0 <- long_data %>% 
    mutate(presence = ifelse(value == 0, 0, 1)) %>% 
    mutate(obesity = as.factor(OBESITYCAT),
           family = as.factor(family))
dim(data0)

# Change obesity column to yes/no 
data0 <- data0 %>% 
    mutate(obesity = ifelse(obesity == "Lean",0,1))




library(geepack)
# need to also specify id 
# I think 'family' might work as a 'cluster?'
zcor <- R[[2]]
dim(zcor)


# Some of the errors that occur
#Error in geese.fit(xx, yy, id, offset, soffset, w, waves = waves, zsca,  : 
#nrow(zcor) need to be equal sum(clusz * (clusz - 1) / 2) for unstructured or userdefined corstr.
#clusz	- integer vector giving the number of observations in each cluster

# some checks to see we have right dimensions... 
clusz <- data0 %>% 
    group_by(family) %>% 
    summarize(n = n())
clusz2 <- na.omit(data0) %>% 
    group_by(family) %>% 
    summarize(n = n())
clusz <- clusz$n
clusz2 <- clusz2$n
sum(clusz * (clusz - 1) / 2)
sum(clusz2 * (clusz2 - 1) / 2)
nrow(zcor)
nrow(na.omit(data0))


### Make our own zcor matrix to account for missing data???



# Sort so observations of the same cluster are together
# Necessary for geeglm() function 
data0_sorted <- data0 %>% 
    arrange(family)


# fit the model... or try to 
geeglm(
    presence ~ obesity,
    data = data0,
    family = "binomial",
    id = family,
    corstr = "userdefined",
    zcor = zcor, 
    waves = paste(time,OTU)
)




# About zcor: 
# Given a n_i x m matrix X_i of covariates, 
# the upper diagonal correlations parameter r_i of the working 
# correlation matrix R_i(alpha) can be written as r_i = X_i alpha 
# 
# The zcor argument takes the concatenated matrices (X1, ... , X_k) as the 
# design matrix for the working correlation 


# About missing data: 
# In case of missing values, the GEE estimates are consistent if the values are missing com- pletely at random (Rubin 1976). The geeglm function assumes by default that observations are equally separated in time. Therefore, one has to inform the function about different sep- arations if there are missing values and other correlation structures than the independence or exchangeable structures are used.
# The waves arguments takes an integer vector that indicates that two observations of the same cluster with the values of the vector of k respectively l have a correlation of rkl.

# How to calculate RA? is it only of these OTUs? do i use the total column? 

