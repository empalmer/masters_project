library(lazyeval)

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
    blank_indeces <- data %>%
        mutate_at(taxa_levels, function(x){ifelse(nchar(x) == 3, .$`#OTU ID`, x)}) 
    
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
