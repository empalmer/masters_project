#' Helper function to return taxa counts by 
#'
#' @param data 
#' @param group_col 
#'
#' @return vector of taxa counts 
#' @export
#'
#' @examples
taxa_group_tally <- function(data, taxa_col) {
    # Done this way to preserve order of the data 
    taxa_counts <- data %>% 
        group_by_at(taxa_col) %>% 
        add_tally() %>% 
        select(n) %>% 
        distinct()
    return(taxa_counts)
}


#' Helper function to treat blank taxa as distinct
#' 
#' Taxa with black levels have the form "x__;" These should all count as separate 
#' level ids when making the tree. Replace these values with the Unique OTU ids 
#' so counts are correct. 
#' 
#' This should be the first step. 
#'
#' @param ranks 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
blank_taxa_fill <- function(data, taxa_levels = c("k","p","c","o","f","g","s"), otu_id_col) {
    #fix later to be more general..., with a general OTU id col
    # Blank cols will be of the form x__
    # blank_indeces <- data %>%
    #     mutate_at(taxa_levels, function(x){ifelse(nchar(x) == 3, .$`#OTU ID`, x)}) 
    #data <- head_gut %>% relocate(c(all_of(otu_id_col),"k","p","c","o","f","g","s"))
    
    ## FIX THIS LATER to be efficient... 
    # if the taxa at the given level is "blank" fill with the higher level value. 
    blank_indeces <- data %>% 
        mutate(k = ifelse(nchar(k) == 3, `#OTU ID`, k)) %>% 
        mutate(p = ifelse(nchar(p) == 3, k, p)) %>% 
        mutate(c = ifelse(nchar(c) == 3, p, c)) %>%
        mutate(o = ifelse(nchar(o) == 3, c, o)) %>%
        mutate(f = ifelse(nchar(f) == 3, o, f)) %>%
        mutate(g = ifelse(nchar(g) == 3, f, g)) %>%
        mutate(s = ifelse(nchar(s) == 3, g, s)) %>%
    
    return(blank_indeces)
}




#' Tax2cor: given data consisting of separate taxa columns, return a list appropriate for cor2zcor structure
#'
#' @param data with columns "k","p","c","o","f","g","s" Data shoudl have blank taxa names filled 
#' @param ranks A vector of the ranks to extract structure from
#'
#' @return a list that gives the tree structure 
#' @export
#'
#' @examples
tax2cor <- function(data, taxa_levels = c("k","p","c","o","f","g","s")){
    all <- map(taxa_levels, ~taxa_group_tally(data,.x))
    all
    return(all)
} 





#' Sort by the levels, alphabetically 
#' This should be the final order of the data 
#'
#' @param taxa_levels 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
sort_by_taxa <- function(data,taxa_levels = c("k","p","c","o","f","g","s")) {
    sorted <- data %>% 
        arrange_at(taxa_levels)
    return(sorted)
    
}





#' Add up observations at a taxa level 
#'
#' @param data Data, where blank values the specified level are filled with unique placeholders
#' @param taxa_level at what level should we sum? Default to species 
#' @param count_cols vector of the sample columns to sum over
#'
#' @return
#' @export
#'
#' @examples
sum_at_taxa_level <- function(data, taxa_level = "s", ...) {
    summed_data <- data %>% 
        group_by_at(taxa_level) %>% 
        summarize_if(is.numeric,sum)
    return(summed_data)
}


#test_sum <- sum_at_taxa_level(gut_subset[1:30,], taxa_level = c("k","p","c","o","f","g","s"))




#### Functions for converting zcor for missing/0 data 


#' Filter R missing
#' 
#' For a given cluster, filter the corresponding integrative correlation matrix based on the 
#' missing observations in the data. 
#' This function will be looped over all of the clusters in the data 
#'
#' @param clust_index an index from 0 to number of clusters 
#' @param data data that contins missing values 
#' @param max_size size of a full cluster
#' @param zcor_1_full 
#' @param R_full 
#' @param value 
#'
#' @return A data frame containing the adjusted zcor for one cluster
#' @export
#'
#' @examples
filter_R_missing <- function(clust_index, data, max_size, zcor_1_full, R_full, value){
    # Filter the data to be only the rows on the focused on cluster 
    # TODO later fix to be general value column 
    row_start <- clust_index * max_size + 1
    row_end <- row_start + max_size - 1
    #data_cluster <- data[row_start:row_end,]$value
    data_cluster <- data[row_start:row_end,value]
    
    #Figure out which rows of the family are missing (NA values)
    #Of the corresponding rows and columns of R, change to be a 
    # dummy value (-2) chosen, but really any value is ok
    missing_values <- is.na(data_cluster)
    R_full[missing_values,] <- -2
    R_full[,missing_values] <- -2
    
    #convert R matrix into a vector, taking only the lower triangular part of 
    # the matrix. This will be used for indexing
    matrix_as_vector <- as.numeric(R_full[lower.tri(R_full)])
    
    # filter out the rows of the full zcor matrix where there is an NA value (-2)
    zcor_reduced <- zcor_1_full[matrix_as_vector != -2, ]
    
    return(zcor_reduced)
}



library(Matrix)
#' Adjust zcor
#'
#' Function that calls the above function for every family value 
#'
#' @param max_size Size of a "full" cluster
#' @param dim Correlation structure for OTUs
#' @param par Repeated measure correlation structure
#' @param n Number of clusters
#' @param data Data frame including missing observations - must be ordered in terms of cluster 
#' 
#'
#' @return Adjusted Integrative Correlation Matrix R that accounts for missing data
#' 
adjust_zcor <- function(data, max_size, n, R, zcor, value){
    
    clust_index <- 0:(n-1)
    # For every cluster in the data, call the above function to filter
    zcor_filter_missing <- map(clust_index, ~filter_R_missing(.x , data, max_size, zcor, R, value))
    #rbind results of all clusters to create 1 adjusted corrected zcor. 
    zcor_combined_filtered_missing <- reduce(zcor_filter_missing, rbind)
    return(zcor_combined_filtered_missing)
}    
