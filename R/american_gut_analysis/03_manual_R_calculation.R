library(tidyverse)
source(here::here("R","american_gut_analysis","01_taxonomy_to_correlation_list.R"))
source(here::here("GEE_Chen2020","cor2zcor.r"))
## load in ----- 
# try at species level 
#filtered_data <- read_rds(here::here("Data","American_Gut","sum_genus_filtered.rds"))
# try at genus level on selected sample sites
#filtered_data2 <- read_rds(here::here("Data","American_Gut","sum_species_level.rds"))

# try at thresholded genus level at selected sample sites 
#filtered_data <- read_rds(here::here("Data","American_Gut","genus_threshold_10.rds"))
filtered_data <- read_rds(here::here("Data","American_Gut","samples100.rds"))

## get taxa list ---- 
taxa_count_list <- tax2cor(filtered_data, taxa_levels = c("k","p","c","o","f","g"))
taxa_count_list


#extract just the counts - the dimension list to give to the cor2zcor fxn. 
dim <- map(taxa_count_list, 2)
dim


#dont need lowest level of all 1s. 
dim <- dim[1:5]
dim

# Number of levels 
K <- length(dim)
K
# Number of total OTUs 
N <- sum(dim[[K]])
N
library(Matrix)
count = dim
if (K > 1) {
    for (k in (K - 1):1) {
        J = length(dim[[k]])
        count[[k]] = rep(1, J)
        for (j in 1:J) {
            index = 0
            while (sum(dim[[k + 1]][sum(count[[k]][1:j]):(sum(count[[k]][1:j]) + index)]) <=
                   dim[[k]][j]) {
                index = index + 1
                if (sum(dim[[k + 1]][sum(count[[k]][1:j]):(sum(count[[k]][1:j]) + index -
                                                           1)]) == dim[[k]][j])
                    break
            }
            count[[k]][j] = index
        }
    }
}


index.level = rep(1, K)
for (k in K:1) {
    J = length(dim[[k]])
    squares = list()
    index.taxa = 0
    for (j in 1:J) {
        if (count[[k]][j] > 1) {
            index.taxa = index.taxa + 1
        }
        squares.new = list(matrix(
            index.level[k] + index.taxa,
            nrow = dim[[k]][j],
            ncol = dim[[k]][j]
        ))
        squares = c(squares, squares.new)
    }
    
    BD = as.matrix(bdiag(squares))
    diag(BD) = 1
    if (k < K) {
        if (sum(dim[[k]]) != N) {
            stop("Inconsistent dimensions: total number of OTUs unequal at each level")
        }
        for (s in 1:N) {
            for (t in 1:N) {
                if (structure.OTU[s, t] > 0 & structure.OTU[s, t] < BD[s, t]) {
                    BD[s, t] = structure.OTU[s, t]
                }
            }
        }
    }
    if (k > 1) {
        index.level[k - 1] = max(BD)
    }
    structure.OTU = BD
}

dim(structure.OTU)
view_otu <- cbind(filtered_data[,1:6],structure.OTU)
View(view_otu)
#write_rds(structure.OTU,here::here("Data","American_Gut","R_genus_threshold.rds") )
write_rds(structure.OTU,here::here("Data","American_Gut","R_100.rds") )
# number of correlations to estimate: 
max(structure.OTU) - 1


# Calculate the zcor matrix: 
max_structure_OTU <- max(structure.OTU) - 1
N
choose(N,2)
# for example, use n = 1 
n <- 1
zcor = matrix(nrow = n*choose(N , 2), ncol = max_structure_OTU )
              #ncol = max(structure.OTU) - 1)
#for (i in 1:(max(structure.OTU) - 1)) {
for (i in 1:max_structure_OTU ) {
    zcor[, i] = rep(as.numeric(structure.OTU[lower.tri(structure.OTU)] == i + 1),n)
}
dim(zcor)
# This is the zcor for 1 observation. 
#zcor = zcor[, colSums(zcor != 0) > 0]


# zcor for all samples: number of samples: 
# - 7 cols for the taxa and ids 
n_samples <- ncol(filtered_data) - 7 
n_samples
    
# initialize 
#check dimension 
n_samples*choose(N , 2) * max_structure_OTU
2^31

# CAN initialize. 
zcor = matrix(nrow = n_samples*choose(N , 2), ncol = max_structure_OTU )


# But this for loop won't run. 
for (i in 1:max_structure_OTU ) {
    zcor[, i] = rep(as.numeric(structure.OTU[lower.tri(structure.OTU)] == i + 1),n_samples)
}
dim(zcor)


write_rds(zcor,here::here("Data","American_Gut","zcor_100_samples.rds") )




# Will purr version work? 
# start 
lower_tri <- lower.tri(structure.OTU)
lower_tri_vals <- structure.OTU[lower_tri]



zcor_full <- map_dfc(1:53, ~rep(as.numeric(lower_tri_vals == .x + 1), n_samples))


write_rds(zcor_full,here::here("Data","American_Gut","zcor_full.rds") )

zcor_full <- read_rds(here::here("Data","American_Gut","zcor_full.rds") )

dim(zcor_full)[1]*dim(zcor_full)[2]
str(zcor_full)
zcor_mat <- as.matrix(zcor_full)
