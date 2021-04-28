library(tidyverse)
library(readxl)
twindata <- read_excel("Data/twindata.xls")
# Load cor2zcor function 
source(here::here("GEE_Chen2020","cor2zcor.r"))


# What is total read column? 
colnames(twindata)

# Looks like there is not an even sequencing depth? 
rowSums(twindata[,10:23])



########## Exploring the correlation structure 
# OTU correlation structure - KNOWN, shown in paper
dim <- list(9,c(4,1,4))

# Gamma matrix - same as in paper 
# Shows correlations between OTUs based on taxonomic tree
gamma <- cor2zcor(1,dim = list(9,c(4,1,4)),1,corstr="exchangeable",otustr="exchangeable")[[1]] -1
gamma
dim(gamma)
# Omega matrix - for 2 twins and 2 timepoints
# Shows 
omega <- cor2zcor(1,dim =1,c(2,2),corstr="exchangeable",otustr="exchangeable")[[1]]-1
omega
dim(omega)

# Integrative correlation matrix 
# R matrix 9x4 is dim 36 - same as paper 
# n is 54 since we have 54 families or clusters 
# This R matrix will have to be adjusted for missing data 
# This is just for initial exploration 
# so n = 1 to show the structure
R_oneobs <- cor2zcor(1,dim = list(9,c(4,1,4)),c(2,2),corstr="exchangeable",otustr="exchangeable")
dim(R_oneobs[[2]])
#34020
dim(R_oneobs[[1]])



## source cleaned dataset from clean_twindata.R. 
# includes a function clean_twinexample_data() that returns a dataframe 
# with correct structure and cleaned
# Has columns: id, family, time, twin_id, subject_id, OTU, value, obesity 
# This column will be sorted by Family 
# within family it will be sorted .... 
source(here::here("R","GEE","clean_twindata.R"))

clean_twindata <- clean_twinexample_data()
dim(clean_twindata)




### Convert data into presence/absence
# Presence/absence  data 
data0 <- clean_twindata %>% 
    mutate(presence = ifelse(value == 0, 0, 1)) %>% 
    arrange(family, OTU, time)
dim(data0)





# Fit model using GEEPACK
library(geepack)
# Id: 
# cluster - family we have 54 families and they "should" have 36 observations
# some missing 
# TODO - Wave indication of correlation parameter for twin/time/otu combo
# mutate(cor_id = paste0(twin_id, time, OTU))
#Initial idea for waves? Figure out later 
# if it was the unique combo we woudl have too many correlations to estimate

### Reorder the clusters in the data to have all the same OTUs together (done in twin function)
# order by family (cluster), otu, sample (twin), time 
#data0_sorted <- data0 %>%
#    mutate(OTU = factor(OTU, levels = c("Blautia",Coprococcus)))
#    fct_relevel(OTU, "Blautia")
    




## Load in corrected Integrative correlation matrix "R"
# That accounts for missing data 
source(here::here("R","GEE","create_zcor.R"))

clusz <- get_clusz(data0, family)
zcor <- adjust_zcor(data0, max_size = 36, dim = list(9,c(4,1,4)), par = c(2,2), n = 54)

nrow(zcor)
sum(clusz * (clusz - 1) / 2)


# fit the model...
mod1 <- geeglm(
    presence ~ obesity,
    data = data0,
    family = "binomial",
    id = family,
    corstr = "userdefined",
    zcor = zcor
)
summary(mod1)


R <- cor2zcor(1,dim = list(9,c(4,1,4)),c(2,2),corstr="exchangeable",otustr="exchangeable")[[1]]
View(R)

# How to calculate RA? is it only of these OTUs? do i use the total column? 




