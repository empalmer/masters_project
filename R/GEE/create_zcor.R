


source(here::here("GEE_Chen2020","cor2zcor.r"))
dim <- list(9,c(4,1,4))

# Try with only 2 clusters 
R <- cor2zcor(2,dim,c(2,2),corstr="exchangeable",otustr="exchangeable")

r <- R[[1]]
r
R[[2]]

# I dont think this fixed2Zcor will work... 
family2 <- rep(1:2, each = 36)
zcor_test <- fixed2Zcor(r, id = family2, waves = family2)
zcor_test
nrow(zcor_test)    




# From cor2zcor code 
zcor = matrix(nrow = n * choose(N * M, 2),
              ncol = max(cor.matrix) - 1)


# replace with a vector that contains the cluster length? 
for (i in 1:15)) {
    zcor[, i] = rep(as.numeric(r[lower.tri(r)] == i + 1), n)
}
zcor = zcor[, colSums(zcor != 0) > 0]


rep(7:8, 2:3)




dim=list(2)
cor2zcor(1,dim,c(2),corstr="exchangeable",otustr="exchangeable")[[1]]
cor2zcor(3,dim,c(2),corstr="exchangeable",otustr="exchangeable")[[2]]

# When there are 2 clusters, and each cluster has 4 observations
genZcor(clusz = c(4,4), waves = rep(1:4,2), corstrv = 4 )

# When there are 2 clusters and the 1st has 4 observations and the 2nd has 3
genZcor(clusz = c(4,3), waves = c(1,2,3,4,1,2,3), corstrv = 4 )





########## Create a separate zcor for each family depending on data 

data_twins_only %>% 
    group_by(FAMILY) %>% tally() %>% table()

data_twins 


data_twins_only %>% 
    group_by(FAMILY,HOST_SUBJECT_ID) %>% tally()


table_data <- data_twins_only %>% 
    select(FAMILY,HOST_SUBJECT_ID)

