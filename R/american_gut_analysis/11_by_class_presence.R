 ## Import data: 
library(tidyverse)
# Load in my functions 
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))
# and paper functions 
source(here::here("GEE_Chen2020","cor2zcor.r"))


# load in the pre-filtered data: 
# Sampling depth above 500
# Bacteria Kingdom only 
# Sum to genus level 
# Filter genera sparsity 
# Fecal site only 
# Has data for antibiotic use 
# randomly select 100 samples

# Try both for 100 samples and full samples (3354)
#filtered_data <- read_rds(here::here("Data","American_Gut","antibiotic_otu_table_100.rds"))
filtered_data <- read_rds(here::here("Data","American_Gut","antibiotic_otu_table_all_samples.rds"))

# Full data - transformed with metadata info 
presence_full <- read_rds(here::here("Data","American_Gut","antibiotic_data_all_samples.rds"))



## Calculate taxa tree and dimensions 
taxa_count_list <- tax2cor(filtered_data, taxa_levels = c("k","p","c","o","f","g"))
taxa_count_list

(class_taxa_counts <- taxa_count_list[[3]])


classes <- class_taxa_counts$c
classes
## Lets run this model on only class. 
# Run through on 1 class. do class c__Bacilli - has 21 members. 

bacilli_class <- filtered_data %>% 
    ungroup() %>% 
    filter(c == "c__Bacilli") %>% 
    select(-k,-p)


# recalculate taxa tree and counts 
bacilli_taxa_tree <- tax2cor(bacilli_class, taxa_levels = c("c","o","f","g"))
bacilli_taxa_tree

dim_bacilli <- map(bacilli_taxa_tree,2)
dim_bacilli

# calculate R and zcor from this (for 1 sample)
bacilli_R_zcor <- cor2zcor(n = 1, dim = dim_bacilli, par = 1)
bacilli_R_zcor[[1]]


# load in full transformed data, and then filter only OTUs by OTUid in selected class 
#presence_full <- read_rds(here::here("Data","temp_data","antibiotic_data_100.rds"))




##### Test i
bacilli_names <- as.character(bacilli_class$`#OTU ID`)

presence_bacilli <- presence_full %>% 
    filter(OTU_name %in% bacilli_names)


#calculate full zcor. 
# for logistic model we can just calculate it using cor2zcor since no "missing"
zcor_full_bacilli <- cor2zcor(n = 100, dim = dim_bacilli, par = 1)[[2]]



library(geepack)
# Model fit 
geepack_fit_zcor <- geeglm(presence ~ use_antibiotic_past_year, family = binomial, data = presence_bacilli, 
                           corstr = "userdefined", zcor = zcor_full_bacilli, id = sample_id)
# Runs in less than a minute!! omg 
summary(geepack_fit_zcor)



#################################################
#################################################

class_filtered_data <- filtered_data %>% 
    ungroup() %>% 
    select(-k,-p)



fit_on_class <- function(class){
    class_data <- class_filtered_data %>% 
        filter(c == class)
    
    class_taxa_tree <- tax2cor(class_data, taxa_levels = c("c","o","f","g"))
    dim_class <- map(class_taxa_tree,2)
    
    class_names <- as.character(class_data$`#OTU ID`)
    presence_class <- presence_full %>% 
        filter(OTU_name %in% class_names)

    #calculate full zcor. 
    # for logistic model we can just calculate it using cor2zcor since no "missing"
    # 100 samples
    # now do for full samples 
    
    zcor_full_class <- cor2zcor(n = 3347, dim = dim_class, par = 1)[[2]]

    # Model fit 
    fit_class <- geeglm(presence ~ use_antibiotic_past_year, family = binomial, data = presence_class, 
                               corstr = "userdefined", zcor = zcor_full_class, id = sample_id)
    return(summary(fit_class))
}

fit_on_one <- function(class){
    class_data <- class_filtered_data %>% 
        filter(c == class)
    
    class_names <- as.character(class_data$`#OTU ID`)
    presence_class <- presence_full %>% 
        filter(OTU_name %in% class_names)
    
    # Model fit 
    fit_class <- geeglm(presence ~ use_antibiotic_past_year, family = binomial, data = presence_class, 
                        corstr = "independence", id = sample_id)
    return(summary(fit_class))
}



classes_over_1 <- class_taxa_counts %>% 
    filter(n > 1) %>% 
    pull(c)
classes_over_1
classes_equal_1 <-class_taxa_counts %>% 
    filter(n == 1) %>% 
    pull(c)
classes_under_1

loop_through_classes_over_1 <- map(classes_over_1, fit_on_class)
# The class with c__Clostridia with 45 obs took much longer than the others. Still under 5 minutes though. 
loop_through_classes_over_1

loop_through_classes_equal_1 <- map(classes_equal_1, fit_on_one)


trial_res <- loop_through_classes_equal_1[[1]]
coef <- trial_res$coefficients
coef$`Pr(>|W|)`[2]


get_p_value <- function(summary){
    coefs <- summary$coefficients 
    p_val <- coefs$`Pr(>|W|)`[2]
    return(p_val)
}

all_results <- c(loop_through_classes_over_1,loop_through_classes_equal_1)

write_rds(all_results, here::here("R","american_gut_analysis","class_loop_presence.rds"))

p_vals <- map_dbl(all_results, get_p_value)
p_vals

# Significant results!!! 
sum(p.adjust(p_vals, "fdr") < .05)

p_vals_adjusted <- p.adjust(p_vals, "fdr")


class_labels <- c(classes_over_1,classes_equal_1)
class_labels[p_vals_adjusted < .05]

