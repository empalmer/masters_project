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
(paper_data_sorted <- sort_by_taxa(test_data, taxa_levels = c("o","f","g")))
tax2cor(paper_data_sorted, taxa_levels = c("o","f","g"))
# This works, but it is in a different order than the paper, so the integrative correlation matrix will be different. 
