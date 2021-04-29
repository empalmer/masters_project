# Repeat the process of tidying/formatting the data 
# Not so many steps needed for non GEE models 

# use this one for presence/absense
cleaned_twin_data <- twindata %>% 
    mutate(all_otu_sum = rowSums(twindata[,10:23])) %>% # column used for alternative RA definition
    filter(ZYGOSITY != "NA") %>% # remove mothers 
    select(all_of(c(meta_vars <- c("HOST_SUBJECT_ID", "SAMPLEDATE", "ZYGOSITY",
                                   "FAMILY", "OBESITYCAT","total_read", "all_otu_sum","Blautia","Coprococcus","Lachnospira",
                                   "Roseburia","Eubacterium","Ruminococcus",
                                   "Faecalibacterium","Oscillospira","Ruminococcus.1")))) %>% # only use 9 OTUs
    pivot_longer(cols = `Blautia`:`Ruminococcus.1`, # convert data to be in a long format
                 names_to = "OTU")  %>% 
    mutate(presence = ifelse(value == 0, 0, 1), # define presence/absence data 0 if zero count 1 otherwise
           obesity = ifelse(OBESITYCAT == "Lean",0,1), # convert obesity (lean/overweight/obese) to binary lean vs overweight/obese
           ra = value/total_read, # one definition of RA count/total read column. 
           ra2 = value/all_otu_sum # another definition of RA count/sum of OTUs in data 
           ) %>% 
    mutate(logt = -log10(ra)) # calculate transformed -log10 of RA - normally distributed


# use for 1 part models (W)
# replace with 5 so -log10(0) isnt infinity
# unclear from paper which value was used here... 
replace_0 <- cleaned_twin_data %>% 
    mutate(logt = ifelse(value == 0,5,logt )) # replace -log10 (0) with 5 

# use for 2 part models (Y)
# replace 0 counts with NA so wont be included in model
remove_0 <- cleaned_twin_data %>% 
    mutate(logt = ifelse(value == 0, NA, logt),
           obesity = ifelse(value == 0, NA, logt))


# Commented model code is from the simulation in the paper
# X is 0/1 obesity 
# Y is 0/1 presence absence 
# Z is -log10 RA when 0s are removed
# W is -log10 RA with 0 replaced by some value
# data: 
#data0: presence/absence ~ obesity 
#remove_0abund: counts of zero are removed from the data *y = logt
#replace_0abund: counts of zero have log transformed value of 6 (y = logt)


## Model 1: GEE0
# GEE model on presence/absense data 
# result beta = -.511
# fit1 = geeglm(
#     Y ~ X,
#     family = "binomial",
#     id = id,
#     corstr = "userdefined",
#     zcor = zcor
# )
# included in twinexample file

## Model 2: GEE+
# GEE model on non-zero RA (linear regression)
# result beta = -0.017
# fit2 = geeglm(
#     Z ~ X,
#     family = "gaussian",
#     id = id,
#     corstr = "userdefined",
#     zcor = zcor.reduced
# )
# included in twinexample file

## Model 3 is combined GEE0 and GEE+




# This should be the easiest model to compare to ... (mod4_B)
# mod4_B should be the same - as it doesn't depend on how RA is calculated 
# But stil different results? 
## Model 4: 2P_ind 
# two part independence, assuming no correlation, logistic 
# result beta0: -.496, beta+: 0.14
# my results: beta0: -0.884, beta+: 0.000753
#fit4.B = glm(Y ~ X, family = "binomial")
#fit4.N = lm(Z ~ X)

mod4_B <- glm(presence ~ obesity , data = cleaned_twin_data, family = "binomial")
mod4_N <- lm(logt ~ obesity, data = remove_0abund)
summary(mod4_B)
summary(mod4_N)




## Model 5: 1P_GEE 
# same correlation structure, only one GEE linear model on all data 
# result beta = -0.041
# fit5 = geeglm(
#     W ~ X,
#     family = "gaussian",
#     id = id,
#     corstr = "userdefined",
#     zcor = zcor
# )



## Model 6: 1P_ind one part independence, 
#assuming no correlation and only one simple linear model for all 0 and nonzero RAs
# result beta = -0.24
# fit6 = lm(W ~ X)
# this will depend on the value of what 0 counts are replace with. Currently using 5
# I tried 6 and 4 as well, no change 

mod6 <- lm(logt~obesity, data = replace_0abund)
summary(mod6)

# my result beta = 0.08


