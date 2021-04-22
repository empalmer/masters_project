### Try another package? 

install.packages("CRTgeeDR")
library(CRTgeeDR)


geeDREstimation(presence ~ obesity, 
                data = data0_sorted, 
                family = binomial("logit"),
                id = family, 
                corstr = "userdefined", 
                corr.mat = R[[1]], 
                nameTRT = cor_id)