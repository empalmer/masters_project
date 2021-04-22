source(here::here("GEE_Chen2020","cor2zcor.r"))



########## Create a separate zcor for each family depending on data 
data_twins_only %>% 
    group_by(FAMILY) %>% tally() %>% table()
data_twins 
data_twins_only %>% 
    group_by(FAMILY,HOST_SUBJECT_ID) %>% tally()
table_data <- data_twins_only %>% 
    select(FAMILY,HOST_SUBJECT_ID)



########## ##############################
########## ##############################
########## Reduce the zcor matrix #######
#########################################
# clusz vector with the 

# parameters for twin data: 
n <- 54
dim <- list(9,c(4,1,4))
par <- c(2,2)
max_size <- 36


library(Matrix)
adjust_zcor <- function(clusz, max_size, dim, par, n){
    #create the unadjusted zcor, if there were full 
    zcor_unadjusted <- cor2zcor(n, dim, par, corstr = "exchangeable", otustr = "exchangeable")[[2]]
    
    # number of correlation paramters they are estimating 
    # 15 for twin example 
    n.cor <- dim(zcor_unadjusted)[2]
    n.cor
    # initialize 
    zcor.reduced = data.frame()
    
    
    
    
}

for (i in 1:n) {
    cor.matrix = cor2zcor(n,
                          dim,
                          2,
                          corstr = "exchangeable",
                          otustr = "exchangeable")[[1]] - 1
    for (m in 1:M) {
        if (is.na(Z[(i - 1) * M + m])) {
            cor.matrix[m, ] = -2
            cor.matrix[, m] = -2
        }
    }
    temp = as.numeric(cor.matrix[lower.tri(cor.matrix)])
    
    if (n.cor > 1) {
        if (ncol(zcor.reduced) > 0) {
            zcor.reduced = as.matrix(zcor.reduced)
        }
        zcor.reduced = rbind(zcor.reduced, zcor[(i - 1) * length(temp) + 1:length(temp), ][temp !=
                                                                                               -2, ])
    }
    else if (n.cor == 1) {
        zcor.reduced = c(as.numeric(zcor.reduced), zcor[(i - 1) * length(temp) +
                                                            1:length(temp)][temp != -2])
    }
}




