library(tidyverse)
library(readxl)
twindata <- read_excel("Data/twindata.xls")
# Load cor2zcor function 
source(here::here("GEE_Chen2020","cor2zcor.r"))

## Test to see if geeglm() and correlation structure will work 
# if we DO have equal observations in each cluster. Remove families if they do not include complete observations


otus_used <- c("Blautia","Coprococcus","Lachnospira",
               "Roseburia","Eubacterium","Ruminococcus",
               "Faecalibacterium","Oscillospira","Ruminococcus.1")
meta_vars <- c("HOST_SUBJECT_ID", "SAMPLEDATE", "ZYGOSITY",
               "FAMILY", "OBESITYCAT")


data_filter_otus <- twindata %>% 
    select(all_of(c(meta_vars,otus_used)))


data_twins_only <- data_filter_otus%>% 
    dplyr::filter(ZYGOSITY != "NA") %>% 
    mutate(twin_family = paste(HOST_SUBJECT_ID, FAMILY)) %>% 
    select(-c("ZYGOSITY"))



remove_missing_timepoints <- data_twins_only  %>% 
    group_by(FAMILY) %>% filter(n() == 4 )

remove_missing_timepoints <- remove_missing_timepoints %>% 
    filter(FAMILY != 37)

length(unique(remove_missing_timepoints$FAMILY))

# now only 39 clusters 
# dimension should be 39 * 36 = 1404
long_data <- pivot_longer(remove_missing_timepoints, 
    cols = `Blautia`:`Ruminococcus.1`, 
    names_to = "OTU") 
dim(long_data)




# Presence/absence 
data0 <- long_data %>% 
    mutate(presence = ifelse(value == 0, 0, 1)) %>% 
    mutate(obesity = as.factor(OBESITYCAT),
           family = as.factor(FAMILY))
dim(data0)
# Change obesity column to yes/no 
data0 <- data0 %>% 
    mutate(obesity = ifelse(OBESITYCAT == "Lean","No","Yes"))

data0 %>% 
    group_by(OBESITYCAT) %>% tally()
data0 %>% 
    group_by(obesity) %>% tally()


# Removing all missing data info, so 40 clusters now
dim <- list(9,c(4,1,4))
R <- cor2zcor(39,dim,c(2,2),corstr="exchangeable",otustr="exchangeable")

zcor <- R[[2]]
dim(zcor)
nrow(zcor)

# check right dimensions for zcor 
clusz <- data0 %>% 
    group_by(family) %>% 
    summarize(n = n())
clusz <- clusz$n
sum(clusz * (clusz - 1) / 2)


######## THIS IS VERY IMPORTANT! ######### 
data0_sorted <- data0 %>% 
    arrange(family)


## Model 
mod_missing_removed <- geeglm(
    presence ~ obesity,
    data = data0_sorted,
    family = "binomial",
    id = family,
    corstr = "userdefined",
    zcor = zcor
)

mod_missing_removed
summary(mod_missing_removed)



