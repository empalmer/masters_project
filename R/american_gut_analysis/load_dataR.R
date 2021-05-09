
library(tidyverse)


file <- here::here("Data","American_Gut","AG_100nt_taxonomy.txt")
colnames = c("OTU_id", 1:4827,"k","p","c","o","f","g","s")

gut <- read_table2(file, skip = 2, n_max = 50, col_names = colnames)
gut_taxa <- gut %>% 
    select(c("k","p","c","o","f","g","s"))



gut_large <- read_table2(file, skip = 2, col_names = colnames)








