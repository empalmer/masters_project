

# examples for testing 
file <- here::here("Data","American_Gut","AG_100nt_taxonomy.txt")
colnames = c("OTU_id", 1:4827,"k","p","c","o","f","g","s")
gut <- read_table2(file, skip = 2, n_max = 50, col_names = colnames)

gut_taxa <- gut %>% 
    select(c("k","p","c","o","f","g","s"))

data <- gut_taxa


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
    # TODO preserve order 
    taxa_counts <- data %>% 
        group_by_at(taxa_col) %>% 
        tally() 
    return(taxa_counts)
}


#' Helper function to treat blank taxa as distinct
#'
#' @param data 
#' @param taxa_col 
#'
#' @return
#' @export
#'
#' @examples
blank_taxa_fill <- function(data,taxa_col) {
    
    blank_indecis <- data %>%
        mutate(remove_1 = )
    
    
}


#' Tax2cor: given data consisting of separate taxa columns, return a list appropriate for cor2zcor structure
#'
#' @param data with columns "k","p","c","o","f","g","s" 
#' @param ranks A vector of the ranks to extract structure from
#'
#' @return a list that gives the tree structure 
#' @export
#'
#' @examples
tax2cor <- function(data, ranks = c("k","p","c")){
    
    #col_ranks <- c("k","p","c","o","f","g","s")
    all <- map(col_ranks, ~taxa_group_tally(data,.x))
    all
    
    
    ## need to introduce noise/fill in  when there is empty values. 
} 