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
gut_data <- read_rds(here::here("Data","American_Gut","gut_data_taxonomy.rds"))



# try it out 
source(here::here("R","american_gut_analysis","taxonomy_to_correlation_list.R"))

#replace blank levels with the OTU id to count as distinct
gut_full_blank <- blank_taxa_fill(gut_data)
gut_full_blank[,4830:4835]

# sum up at given level default to species 
# This takes FOREVER to run.... 
# TODO make a better function that runs faster to do this step..... 
# Makes sense it takes this long because there are thousands of rows and columsn 
# to group and sum over...
ptm <- proc.time()
gut_full_sum_species_level <- sum_at_taxa_level(gut_full_blank,c("k","p","c","o","f","g","s"))
proc.time() - ptm 

write_rds(gut_full_sum_species_level, 
          here::here("Data","American_Gut","sum_species_level.rds"))

# get the dimensions 
taxa_count_list <- tax2cor(gut_full_sum_species_level, taxa_levels = c("k","p","c","o","f","g","s"))
taxa_count_list


#extract just the counts - the dimension list to give to the cor2zcor fxn. 
dim <- map(taxa_count_list,2)
dim


## Get the integrative correlation matrix from this dimension 
source(here::here("GEE_Chen2020","cor2zcor.r"))
# This also takes ~30 min to run with a large number of OTUs. 
R <- cor2zcor(n = 1, dim = dim, par = 1)[[1]]
R

# Its dimensions will be LARGE: 21828x21828... 

# Number of correlations this will have to estimate
# is it ok that the top level is not 1? 
max(R) - 1













#maybe make into a function later... find which level is the first to have only singletons 
singletons <- map_lgl(dim, function(x){sum(x > 1) > 0})
singletons
c("k","p","c","o","f","g","s")[which.min(singletons)]


# test ordering idea 
k <- c("M","M","M","M","N","N")
p <- c("Z","Z","Z","Y","X","X")
c <- c("B","D","F","A","C","E")
(df_sort_test <- data.frame(k,p,c)[sample(1:6,6),])

(df_sorted <- df_sort_test %>% 
    arrange(k,p,c))
tax2cor(df_sorted, taxa_levels = c("k","p","c"))

# test summing up over at a level. this function should take in ordered data 
(test_sum <- sum_at_taxa_level(test_sorted, "o"))









###### When creating the dimensions for the full American gut dataset, 
# there are several OTUs that are the same at the species level.... 

test_blank_fill %>% 
    select(s) %>% 
    group_by(s) %>% 
    tally() %>% 
    arrange(desc(n))

gut_full %>% 
    filter(s == "s__prausnitzii") 

