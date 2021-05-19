#filtered_data <- read_rds(here::here("Data","American_Gut","genus_threshold_10.rds"))
#zcor_full <- read_rds(here::here("Data","American_Gut","zcor_full.rds"))

library(tidyverse)
filtered_data <- read_rds(here::here("Data","American_Gut","prevalence100.rds"))
zcor <- read_rds(here::here("Data","American_Gut","zcor_100_samples.rds"))


filtered_data <- filtered_data %>% 
    mutate(sample_id2 = factor(sample_id), 
           AGE = as.numeric(AGE), 
           presence0 = ifelse(presence,1,0))




library(geepack)
# over an hour and nothing... 

gendat <- function() {
    id <- gl(5, 4, 20)
    visit <- rep(1:4, 5)
    y <- rnorm(id)
    dat <- data.frame(y, id, visit)[c(-2,-9),]
}


#generating the design matrix for the unstructured correlation
zcor <- genZcor(clusz = table(dat$id), waves = dat$visit, corstrv=4)
# defining the Toeplitz structure 
zcor.toep<-matrix(NA, nrow(zcor),3)
zcor.toep[,1]<-apply(zcor[,c(1,4,6)],1,sum)
zcor.toep[,2]<-apply(zcor[,c(2,5)],1,sum)
zcor.toep[,3]<-zcor[,3]

zfit1 <- geese(y ~ 1,id = id, data = dat,
               corstr = "userdefined", zcor = zcor.toep)

set.seed(88)
dat<-gendat()
dat

#test with other 
fake_zcor <- genZcor(clusz = table(filtered_data$sample_id), waves = filtered_data$OTU_name, corstrv = 3)
fake_zcor
dim(fake_zcor)


geepack_fit <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, id = sample_id2, 
                      corstr = "unstructured")
summary(geepack_fit)

geepack_fit <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, id = sample_id, 
                      corstr = "userdefined", zcor = fake_zcor)


geepack_fit <- geeglm(presence ~ AGE, family = binomial, data = filtered_data, 
                      corstr = "userdefined", zcor = zcor, id = sample_id2)




library(geeM)

geem_fit <- geem(presence ~ AGE, data = filtered_data, family = binomial, id = sample_id,
     corstr = "independence")


geem(presence ~ AGE, data = filtered_data, family = binomial, id = sample_id,
     corstr = "userdefined", corr.mat = structure.OTU)


