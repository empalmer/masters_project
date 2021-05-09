


# testing the tax2cor function
# do the one in the paper 

o <- rep("o__Clostridiales",9)
f <- c(rep("f__Lachno",4),"f__Eubacter",rep("f__Ruminococcaceae",4))
g <- c(
    "g__Blautia",
    "g__Coprococcus",
    "g__Lachnospira",
    "g__Roseburia",
    "g__Eubacterium",
    "g__Ruminococcus",
    "g__Faecalibacterium",
    "g__Oscillospira",
    "g__Ruminococcus.1"
)


test_data <- data.frame(o,f,g)
(paper_data_sorted <- sort_by_taxa(test_data, taxa_levels = c("o","f","g")))
tax2cor(paper_data_sorted, taxa_levels = c("o","f","g"))

# This works, but it is in a different order than the paper, so the integrative correlation matrix will be different. 





### Examples with the gut data, running throught the pipeline of fxns in the right order 
# Load the data 
file <- here::here("Data","American_Gut","AG_100nt_taxonomy.txt")
colnames = c("OTU_id", 1:4827,"k","p","c","o","f","g","s")
gut <- read_table2(file, skip = 2, n_max = 50, col_names = colnames)

gut_taxa <- gut %>% 
    select(c("OTU_id",))


# try it out 

source(here::here("R","american_gut_analysis","taxonomy_to_correlation_list.R"))

#replace blank levels with the OTU id to count as distinct
test_blank_fill <- blank_taxa_fill(gut)

# should have sorted first... 
test_list <- tax2cor(test_blank_fill, taxa_levels = c("k","p","c","o","f","g","s"))
test_list
#extract just the counts
dim <- map(test_list,2)
dim


## Get the integrative correlation matrix from this dimension 
source(here::here("GEE_Chen2020","cor2zcor.r"))

R <- cor2zcor(n = 1, dim = dim[1:6], par = 1)[[1]]
R
# Number of correlations this will have to estimate
# Would need to do some adjustments to make sure the lowest level is only singletons 
# adding up sublevels - see sum_at_taxa_level fxn. 
max(R) - 1


#maybe make into a function later... find which level is the first to have only singletons 
singletons <- map_lgl(dim, function(x){sum(x > 1) > 0})
singletons
c("k","p","c","o","f","g","s")[which.min(singletons)]


#order the data - this doesnt work yet...
test_sorted <- sort_by_taxa(test_blank_fill)
test_blank_fill[,4825:4835]
test_sorted[,4825:4835]

#This order will not match 


# test ordering idea 
k <- c("M","M","M","M","N","N")
p <- c("Z","Z","Z","Y","X","X")
c <- c("B","D","F","A","C","E")
(df_sort_test <- data.frame(k,p,c)[sample(1:6,6),])

df_sorted <- df_sort_test %>% 
    arrange(k,p,c)
tax2cor(df_sorted, taxa_levels = c("k","p","c"))

df_sorted %>% 
    group_by(p) %>% 
    add_tally() %>% 
    select(n) %>% 
    distinct()


# test summing up over at a level. this function should take in ordered data 
test_sum <- sum_at_taxa_level(test_sorted, "o")


