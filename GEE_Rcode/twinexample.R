library(tidyverse)
library(readxl)
twindata <- read_excel("twindata.xls")



# OTU correlation structure - KNOWN, shown in paper
dim <- list(9,c(4,1,4))


# Gamma matrix - same as in paper 
cor2zcor(1,dim,1,corstr="exchangeable",otustr="exchangeable")[[1]] -1
# Omega matrix
cor2zcor(1,dim =1,c(2,2),corstr="exchangeable",otustr="exchangeable")[[1]]-1


# R matrix 9x4 is dim 36 - same as paper 
R <- cor2zcor(1,dim,c(2,2),corstr="exchangeable",otustr="exchangeable")[[1]]
dim(R)


# Y is OTU 
# Y should have dim K x J where 
# K <- 54
# J <- 36
54*36
#  = 1944








## Need to remove OTUs not included in the correlation tree
twin_select <- twindata[,-c(1,2,3,4,8)] %>% 
    select(-total_read) %>% 
    select(-Bacteroides) %>% 
    select(-Prevotella) %>% 
    select(-Alistipes) %>% 
    select(-`Erysipelotrichaceae..g__Clostridium`) %>% 
    select(-Lachnospiraceae..g__Clostridium)
colnames(twin_select)
dim(twin_select)    

# seems like there are parents in the model??? ie 6 observations in 1 family 
# Remove parents from data 
twin_noparents <- twin_select %>% 
    dplyr::filter(ZYGOSITY != "NA")


# Missing twin data 
twin_select %>% 
    group_by(HOST_SUBJECT_ID) %>% summarize(n = n())
twin_select %>% 
    group_by(FAMILY) %>% summarize(n = n())



nrow(twindata)
twin_filter <- twindata %>% 
    dplyr::filter(!is.na(ZYGOSITY))
nrow(twin_filter)

# pivot to be correct data structure, each row is a 
long_twin <- pivot_longer(twin_noparents, 
                          cols = `Blautia`:`Ruminococcus.1`, 
                          names_to = "OTU")
dim(long_twin)



# Separate into 0 and + parts 


# Fit model 