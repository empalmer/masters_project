alphas <- normal_antibiotic_mod$geese[[2]]

alpha_df <- data.frame(alpha = c(0,alphas), id = 0:55)

R <- r_zcor_antibiotic_100[[1]]



r_vec <- as.vector(R)
alphas_0 <- c(1,alphas)
alpha_vec <- alphas_0[r_vec]
#alpha_vec[alpha_vec > 1] <- 1
alpha_mat <- matrix(alpha_vec, nrow = 164)

#install.packages("corrplot")
corrplot::corrplot(alpha_mat, method = "color", is.corr = FALSE,tl.pos='n')



library(tidyverse)
####### 
taxa_levels <- antibiotic_otu_table_100 %>% select(k,p,c,o,f,g)

R_w_taxa <- cbind(taxa_levels,R)

alpha_w_taxa <- cbind(taxa_levels,alpha_mat)

alpha_df %>% filter(alpha > 1)
which(alpha_mat > 1, arr.ind = T)

R_w_taxa[87,]

e1 <- antibiotic_otu_table_100 %>% 
    filter(f == "f__Ruminococcaceae")

e1_val <- e1 %>% 
    ungroup() %>% 
    select(-k,-p,-c,-o,-f,-g,`#OTU ID`)
cor(e1_val)


corrplot::corrplot(cor(e1_val))
