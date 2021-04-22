#### Testing genZcor 
# Try with one observation 
clusz <- 5


gendat <- function() {
    id <- gl(5, 4, 20)
    visit <- rep(1:4, 5)
    y <- rnorm(id)
    dat <- data.frame(y, id, visit)[c(-2,-9),]
}

set.seed(88)
dat<-gendat()
dat

zcor <- genZcor(clusz = table(dat$id), waves = dat$visit, corstrv=4)
zcor
# defining the Toeplitz structure 
zcor.toep<-matrix(NA, nrow(zcor),3)
zcor.toep[,1]<-apply(zcor[,c(1,4,6)],1,sum)
zcor.toep[,2]<-apply(zcor[,c(2,5)],1,sum)
zcor.toep[,3]<-zcor[,3]

zfit1 <- geese(y ~ 1,id = id, data = dat,
               corstr = "userdefined", zcor = zcor.toep)


zfit2 <- geeglm(y ~ 1,id = id, data = dat,
                corstr = "userdefined", zcor = zcor.toep)


### Try with our data: 
remove_missing <- na.omit(data0_sorted)
clusz <- table(remove_missing$family)
sum(clusz * (clusz - 1) / 2)
myzcor <- genZcor(clusz = table(remove_missing$family),
                  waves = as.numeric(as.factor(remove_missing$cor_id)) ,
                  corstrv = 4)
myzcor
nrow(myzcor)

#test 
#2 families 
family2 <- data0_sorted %>% 
    filter(family %in% c(1,11)) %>% 
    mutate(waves = rep(sample(1:16,36,replace = T),2))
R <-  cor2zcor(1,list(9,c(4,1,4)),c(2,2),corstr="exchangeable",otustr="exchangeable")
R <- R[[1]]
dim(R)
R[lower.tri(R)]


as.numeric(cor.matrix[lower.tri(cor.matrix)] == i + 1)


length(unique(family2$cor_id))
zcor_2families <- genZcor(clusz = c(36,36), 
                          waves = as.numeric(as.factor(family2$cor_id)),
                          corstr = 4)
zcor_2families %>% head()
dim(zcor_2families)
nrow(zcor_2families)
clusz36 <- c(36,36)
sum(clusz36 * (clusz36 - 1) / 2)


geeglm(
    presence ~ obesity,
    data = family2,
    family = "binomial",
    id = family,
    corstr = "userdefined",
    zcor = zcor_2families,
    waves = 
)
