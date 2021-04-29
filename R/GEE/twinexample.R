library(tidyverse)
library(readxl)
twindata <- read_excel("Data/twindata.xls")
# Load cor2zcor function 
source(here::here("GEE_Chen2020","cor2zcor.r"))



###### Understanding the data
# What is total read column? 
colnames(twindata)

# The sum of the counts of all the OTUs in the paper does not equal the total read column
rowSums(twindata[,10:23])
twindata$total_read


########## Exploring the correlation structure 
# OTU correlation structure - KNOWN, shown in paper
dim <- list(9,c(4,1,4))

# Gamma matrix - same as in paper 
# Shows correlations between OTUs based on taxonomic tree
gamma <- cor2zcor(1,dim = list(9,c(4,1,4)),1,corstr="exchangeable",otustr="exchangeable")[[1]] -1
gamma
dim(gamma)
# Omega matrix - for 2 twins and 2 timepoints
# Same as paper
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
R_oneobs[[1]]
dim(R_oneobs[[2]])
#34020
dim(R_oneobs[[1]])



## source cleaned dataset from clean_twindata.R. 
# includes a function clean_twinexample_data() that returns a dataframe 
# with correct structure and cleaned
# Has columns: id, family, time, twin_id, subject_id, OTU, value, obesity 
# This column will be sorted by Family 
# within family it will be sorted by OTU 
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
# Id: cluster - family we have 54 families and they "should" have 36 observations - but missing obs

## Load in corrected Integrative correlation matrix "R"
# That accounts for missing data 
source(here::here("R","GEE","create_zcor.R"))
zcor <- adjust_zcor(data0, max_size = 36, dim = list(9,c(4,1,4)), par = c(2,2), n = 54)

# check that we have the correct dimension of zcor to account for missing structure
clusz <- get_clusz(data0, family)
nrow(zcor)
sum(clusz * (clusz - 1) / 2)



# fit the presence/absence part of the model 
mod1 <- geeglm(
    presence ~ obesity,
    data = data0,
    family = "binomial",
    id = family,
    corstr = "userdefined",
    zcor = zcor,
    scale.fix = TRUE
)
summary(mod1)

### ... not the same result as the paper (beta0 = -.511)



