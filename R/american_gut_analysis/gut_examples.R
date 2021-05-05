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
tax2cor(test_data, ranks = c("o","f","g"))




### Examples with the gut data, running throught the pipeline of fxns in the right order 
# Load the data 
file <- here::here("Data","American_Gut","AG_100nt_taxonomy.txt")
colnames = c("OTU_id", 1:4827,"k","p","c","o","f","g","s")
gut <- read_table2(file, skip = 2, n_max = 50, col_names = colnames)

gut_taxa <- gut %>% 
    select(c("OTU_id",))


# try it out 

#replace blank levels with the OTU id to count as distinct
test_blank_fill <- blank_taxa_fill(gut)

test_list <- tax2cor(test_blank_fill, taxa_levels = c("k","p","c","o","f","g","s"))
test_list
#extract just the counts
dim <- map(test_list,2)
dim


#maybe make into a function later... find which level is the first to have only singletons 
singletons <- map_lgl(dim, function(x){sum(x > 1) > 0})
singletons
c("k","p","c","o","f","g","s")[which.min(singletons)]

#order the data 
test_sorted <- sort_by_taxa(test_blank_fill)

# test summing up over at a level. this function should take in ordered data 
test_sum <- sum_at_taxa_level(test_sorted, "o")


