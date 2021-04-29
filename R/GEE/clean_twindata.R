library(tidyverse)
library(readxl)



# returns cleaned and formatted twin data for this specific example
# made into a function so it is easier to load elsewhere
clean_twinexample_data <- function() {
    twindata <- read_excel("Data/twindata.xls")
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
    otus_used <- c(
        "Blautia",
        "Coprococcus",
        "Lachnospira",
        "Roseburia",
        "Eubacterium",
        "Ruminococcus",
        "Faecalibacterium",
        "Oscillospira",
        "Ruminococcus.1"
    )
    # which of the other columns are needed 
    meta_vars <- c(
        "HOST_SUBJECT_ID",
        "SAMPLEDATE",
        "ZYGOSITY",
        "FAMILY",
        "OBESITYCAT",
        "total_read"
    )
    
    # Select only the neede OTUs and other columns 
    data_filter_otus <- twindata %>%
        select(all_of(c(meta_vars, otus_used)))


    
    # Also remove 1 observation where there is not a measured obesity value
    twins_no_na <- data_filter_otus %>%
        filter(!is.na(OBESITYCAT))
    
    # Or impute it since the obesity category is likely the same as the other timepoint?
    # just trying to make results match - but this doesnt change anything
    #twins_no_na <- data_filter_otus %>%
    #    mutate(OBESITYCAT = ifelse(HOST_SUBJECT_ID == "TS131","Obese", OBESITYCAT))
    #twins_no_na %>%
    #    filter(FAMILY == 37)
    
    
    
    # Mothers were also sampled in the study - remove for this analysis.
    # remove zygosity column as not important for later work
    # Create an ID column that includes the subject ID and the family
    # arranging by family  is necessary since the order of the data is important 
    # for the GEE function 
    data_twins_only <- twins_no_na %>%
        dplyr::filter(ZYGOSITY != "NA") %>%
        mutate(twin_family = paste(HOST_SUBJECT_ID, FAMILY)) %>% #id column
        select(-c("ZYGOSITY")) %>%
        arrange(FAMILY)
    
    
    
    # Introduce NA rows for missing data
    # This is done so we have the same cluster size and rows for every observation
    # even if it is missing
    # Missing NA values help later for correcting the zcor matrix
    # This is done to make the final data have the dimensinos the paper says it will
    
    expected_combos <- expand_grid(
        twin_family = unique(data_twins_only$twin_family),
        time = c("TimePoint1", "TimePoint2")
    ) %>%
        mutate(twin_id = rep(c("A", "A", "B", "B"), 54))
    
    data_include_missing <- left_join(
        expected_combos,
        data_twins_only,
        by = c("twin_family" = "twin_family",
               "time" = "SAMPLEDATE")
    ) %>%
        separate(twin_family, c("id", "family"))
    
    
    # pivot to be correct data structure, each row is a
    # different observation on timepoints or OTU value
    # and where we do so it is the "right dimension
    long_data <- pivot_longer(data_include_missing,
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
        mutate(obesity = ifelse(OBESITYCAT == "Lean", 0, 1),
               family = as.factor(family)) %>%
        rename(subject_id = HOST_SUBJECT_ID) %>%
        select(-c(OBESITYCAT, FAMILY))
    


    ## sort to be by OTU - since this is how the correlation matrix is specified
    long_data_neat <- long_data_neat %>%
        arrange(family, OTU)
    
    return(long_data_neat)
}