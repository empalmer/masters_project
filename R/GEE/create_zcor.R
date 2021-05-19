source(here::here("GEE_Chen2020","cor2zcor.r"))


# A function to make it easier to check that we have the correct dimension of the R matrix 
#' Get cluster sizes
#'
#' @param data Data in long format. Can contain missing observations
#' @param cluster_col The column in your data that represents the cluster
#'
#' @return A vector indicating how many observations are in each cluster 
#' @export
#'
#' @examples get_clusz(clean_twindata, family)
get_clusz <- function(data, cluster_col){
    clusz <- data %>% 
        na.omit() %>%
        group_by({{cluster_col}}) %>% 
        tally() %>% 
        pull(n)
    return(clusz)
}


########## ##############################
########## ##############################
########## Reduce the zcor matrix #######
#########################################
# Currently a slightly messy function that is specific to this twin data 
# Later will change to be general to any dataset. 


# parameters for twin data: 
# for testing the function 
# source(here::here("R","GEE","clean_twindata.R"))
# clean_twindata <- clean_twinexample_data()
# clusz <- get_clusz(clean_twindata, family)
# n <- length(clusz)
# dim <- list(9,c(4,1,4))
# par <- c(2,2)
# max_size <- 36
# 
# data <- clean_twindata
# cluster_column <- "family"





#' Filter R missing
#' 
#' For a given cluster, filter the corresponding integrative correlation matrix based on the 
#' missing observations in the data. 
#' This function will be looped over all of the clusters in the data 
#'
#' @param clust_index an index from 0 to number of clusters 
#' @param dim OTU dimension used to build full R and zcor
#' @param par repeated measure used to build full R and zcor
#' @param data data that contins missing values 
#' @param max_size size of a full cluster
#' @param n_cor number of correlations to extimate 
#'
#' @return A data frame containing the adjusted zcor for one cluster
#' @export
#'
#' @examples
filter_R_missing <- function(clust_index, dim, par, data, max_size, n_cor){
    # full R matrix when there are no missing observation
    # move this out of this function loop later for computation time 
    zcor_full <- cor2zcor(1, dim, par,corstr = "exchangeable",otustr = "exchangeable")[[2]]
    R <- cor2zcor(1, dim, par, corstr = "exchangeable", otustr = "exchangeable")[[1]] - 1
    
    
    # Filter the data to be only the rows on the focused on cluster 
    # TODO later fix to be general value column 
    row_start <- clust_index * max_size + 1
    row_end <- row_start + max_size - 1
    data_cluster <- data[row_start:row_end,]$value
    
    
    # Dplyr version of the above, change later 
    # Filter the data for the cluster analyzed in this "loop"
    # TODO fix later to be general cluster column name 
    #data_for_cluster <- data %>% 
    #    filter(family== i)
    
    
    #Figure out which rows of the family are missing (NA values)
    #Of the corresponding rows and columns of R, change to be a 
    # dummy value (-2) chosen, but really any value is ok
    missing_values <- is.na(data_cluster)
    R[missing_values,] <- -2
    R[,missing_values] <- -2
    
    
    #convert R matrix into a vector, taking only the lower triangular part of 
    # the matrix. This will be used for indexing
    matrix_as_vector <- as.numeric(R[lower.tri(R)])
    
    # filter out the rows of the full zcor matrix where there is an NA value (-2)
    zcor_reduced <- zcor_full[matrix_as_vector != -2, ]
    
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
adjust_zcor <- function(data, max_size, dim, par, n){
    #create vector of unique cluster_ids 
    #clust_names <- data %>% 
    #    pull({{cluster_column}}) %>% 
    #    unique()
    # Generl

    # create the unadjusted zcor, if there were no missing
    # rows corresponding to missing obs will be removed 
    cors <- cor2zcor(n, dim, par, corstr = "exchangeable", otustr = "exchangeable")
    zcor_unadjusted <- cors[[2]]
    R_unadjusted <- cors[[1]]
    
    # number of correlation parameters they are estimating 
    # corresponds to the number of columns in zcor
    # 15 for twin example 
    n_cor <- ncol(zcor_unadjusted)
    n_cor

    # Fix this later to be cleaner and more general - here 53 is the number of families in the twin data 
    clust_index <- 0:53
    # For every cluster in the data, call the above function to filter
    zcor_filter_missing <- map(clust_index, ~filter_R_missing(.x ,dim, par, data, max_size, n_cor))
    #rbind results of all clusters to create 1 adjusted corrected zcor. 
    zcor_combined_filtered_missing <- reduce(zcor_filter_missing, rbind)
    return(zcor_combined_filtered_missing)
}

















