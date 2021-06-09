## Import data: 
library(tidyverse)
library(geepack)
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
filtered_data <- read_rds(here::here("Data","American_Gut","antibiotic_otu_table_all_samples.rds"))
#filtered_data <- read_rds(here::here("Data","American_Gut","antibiotic_otu_table_100.rds"))


class_filtered_data <- filtered_data %>% 
    ungroup() %>% 
    select(-k,-p)


log_ra_full <- read_rds(here::here("Data","American_Gut","log_ra_antibiotic_all_samples.rds"))
#log_ra_full <- read_rds(here::here("Data","temp_data","log_ra_antibiotic_100.rds"))


log_ra_full <- log_ra_full %>% 
    mutate(OTU_name = factor(OTU_name))


## Calculate taxa tree and dimensions 
taxa_count_list <- tax2cor(filtered_data, taxa_levels = c("k","p","c","o","f","g"))
taxa_count_list

(class_taxa_counts <- taxa_count_list[[3]])

fit_ra_on_class <- function(class){
    class_data <- class_filtered_data %>% 
        filter(c == class)
    
    class_depth <- class_data %>% 
        summarise_if(is.numeric,sum)
    class_depth
    non_blank_sample <- names(class_depth[,(class_depth>0)[1,]])[-1]
    
    class_taxa_tree <- tax2cor(class_data, taxa_levels = c("c","o","f","g"))
    dim_class <- map(class_taxa_tree,2)
    
    class_names <- as.character(class_data$`#OTU ID`)
    ra_class <- log_ra_full %>% 
        filter(OTU_name %in% class_names) %>% 
        filter(sample_id %in% non_blank_sample)
    
    ra_class %>% 
        group_by(sample_id) %>% 
        summarize(sum(is.na(ra)))
    
    #calculate R and zcor for one obs 
    R_zcor_1 <- cor2zcor(n = 1, dim = dim_class, par = 1)
    R1 <- R_zcor_1[[1]]
    zcor1 <- R_zcor_1[[2]]
    
    # if we have exchangeable structure 
    if(max(R1) - 1 > 1){
        #calculate full zcor. 
        # need to adjust based on 0 otu values 
        # for logistic model we can just calculate it using cor2zcor since no "missing"
        # 100 samples
        n_samples <- length(non_blank_sample)
        n_OTUs <- nrow(R1)
        adjusted_zcor <- adjust_zcor(ra_class, n_OTUs, n_samples, R1, zcor1, "ra")
        
        # Model fit 
        
        fit_class <- geeglm(log_ra ~ nonzero_antibiotic, family = gaussian, data = ra_class, 
                            corstr = "userdefined", zcor = adjusted_zcor, id = sample_id)
    }else{
        fit_class <- geeglm(log_ra ~ nonzero_antibiotic, family = gaussian, data = ra_class, 
                            corstr = "exchangeable", id = sample_id, waves = OTU_name)
    }
    
    return(summary(fit_class))
}

fit_ra_on_one <- function(class){
    class_data <- class_filtered_data %>% 
        filter(c == class)
    
    class_names <- as.character(class_data$`#OTU ID`)
    ra_class <- log_ra_full %>% 
        filter(OTU_name %in% class_names)
    
    # Model fit 
    fit_class <- geeglm(presence ~ use_antibiotic_past_year, family = gaussian, data = ra_class, 
                        corstr = "independence", id = sample_id)
    return(summary(fit_class))
}




classes_over_1_ra <- class_taxa_counts %>% 
    filter(n > 1) %>% 
    pull(c)

classes_equal_1_ra <-class_taxa_counts %>% 
    filter(n == 1) %>% 
    pull(c)


fit_ra_on_class("c__Coriobacteriia")
class <- "c__Coriobacteriia"
fit_ra_on_class("c__Actinobacteria")

ra_loop_through_classes_over_1 <- map(classes_over_1_ra, fit_ra_on_class)
ra_loop_through_classes_equal_1 <- map(classes_equal_1_ra, fit_ra_on_one)



get_p_value <- function(summary){
    coefs <- summary$coefficients 
    p_val <- coefs$`Pr(>|W|)`[2]
    return(p_val)
}

all_results_ra <- c(ra_loop_through_classes_over_1,ra_loop_through_classes_equal_1)

write_rds(all_results_ra, here::here("R","american_gut_analysis","class_loop_ra.rds"))


p_vals_ra <- map_dbl(all_results_ra, get_p_value)
p_vals_ra

# No significant results :( 
sum(p.adjust(p_vals_ra, "fdr") < .05)

p_vals_adjusted_ra <- p.adjust(p_vals_ra, "fdr")


class_labels_ra <- c(classes_over_1_ra,classes_equal_1_ra)
class_labels_ra[p_vals_adjusted_ra < .05]


cauchy.test <- 0.5 * tan((0.5 - p_vals_adjusted_ra) * pi) + 0.5 * tan((0.5 - p_vals_adjusted) * pi)
cauchy_adjust <- 1 - pcauchy(cauchy.test)


class_labels_ra[cauchy_adjust< .05]
