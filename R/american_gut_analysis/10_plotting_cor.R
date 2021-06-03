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
