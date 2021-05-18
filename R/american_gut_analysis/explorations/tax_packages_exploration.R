
## RAM package 
install.packages("RAM")
library(RAM)

file2 <- here::here("Data","American_Gut","AG_100nt_taxonomy2.txt")
gut_ram <- fread.OTU(file2)

gut_ram_small <- gut_ram[1:20,]

OTU.recap(list(data = gut_ram))

phylog_taxonomy(gut_ram)
data.clust(gut_ram)
factor.abundance()



## Explore existing packages 
install.packages("taxonomizr")
library(taxonomizr)
makeNewick(gut_taxa)


install.packages("taxize")
library('taxize')

install.packages("taxa")
install.packages("metacoder")
library(taxa)
library(metacoder)

ranks <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
obj <- parse_tax_data(gut_taxa, class_cols = ranks, named_by_rank = T)
