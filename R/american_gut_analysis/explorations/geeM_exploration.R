

# Exploring other GEE packages: 

install.packages("geeM")
library(geeM)



gut_otu_repsonse_presence
colnames(gut_otu_repsonse_presence)
#start 10:35
geem(presence ~ AGE, data = gut_otu_repsonse_presence, family = binomial, id = sample_id,
     corstr = "userdefined", corr.mat = structure.OTU)


# Test this out on the paper's example? 
