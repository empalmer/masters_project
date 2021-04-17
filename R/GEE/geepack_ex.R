
## Examples from tutorials and papers on geepack

library(geepack)
# Data included in geeglm package 
data(respiratory)



# Model from paper 
#fit with the exchangeable working correlation is obtained by
m.ex <- geeglm(outcome ~ baseline + center + sex + treat + age + I(age^2),
                    data = respiratory,
                    id = interaction(center, id),
                    family = binomial, corstr = "exchangeable")


# The design matrix for the Toeplitz working correlation was constructed by first obtaining the design matrix for the unstructured correlation using the genZcor function and then adding the appropriate columns. The use of the waves argument is not necessary because there are no missing observations.
zcor <- genZcor(clusz = c(xtabs(~ id + center, data = respiratory)),
                        waves = respiratory$visit, corstrv = 4)
zcor.toep <- matrix(NA, nrow(zcor), 3) 
zcor.toep[,1] <- apply(zcor[,c(1,4,6)], 1, sum)
zcor.toep[,2] <- apply(zcor[,c(2,5)], 1, sum)
zcor.toep[,3] <- zcor[,3]


m.toep <- geeglm(outcome ~ baseline + center + sex + treat + age + I(age^2),
                      data = respiratory, id = interaction(center, id),
                      family = binomial, corstr = "userdefined", zcor = zcor.toep)










# Other example, that includes missing data 
# from example https://cran.r-project.org/web/packages/geepack/vignettes/geepack-manual.pdf
library(geepack)
timeorder <- rep(1:5, 6)
tvar <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:6, each=5)
uuu <- rep(rnorm(6), each=5)
yvar <- 1 + 2*tvar + uuu + rnorm(length(tvar))
simdat <- data.frame(idvar, timeorder, tvar, yvar)
head(simdat,12)


#Notice that clusters of data appear together in simdat and that observations are ordered (according to timeorder) within clusters.
#We can fit a model with an AR(1) error structure as

mod1 <- geeglm(yvar~tvar, id=idvar, data=simdat, corstr="ar1")
mod1


# If observatios were not ordered according to cluster and time within cluster we
# would get the wrong result:
set.seed(123)
library(doBy)
simdatPerm <- simdat[sample(nrow(simdat)),]
simdatPerm <- orderBy(~idvar, simdatPerm)
simdatPerm <- simdatPerm[order(simdatPerm$idvar),]
head(simdatPerm)

# Notice that in simdatPerm data is ordered according to subject but the time ordering within subject is random.
# Fitting the model as before gives
mod2 <- geeglm(yvar~tvar, id=idvar, data=simdatPerm, corstr="ar1") > mod2



### R code from vignette source 'geepack-manual.Rnw'

###################################################
### code chunk number 1: geepack-manual.Rnw:16-19
###################################################
require( geepack )
prettyVersion <- packageDescription("geepack")$Version
prettyDate <- format(Sys.Date())


###################################################
### code chunk number 2: geepack-manual.Rnw:73-75
###################################################
library(geepack)
citation("geepack")


###################################################
### code chunk number 3: geepack-manual.Rnw:92-100
###################################################
library(geepack)
timeorder <- rep(1:5, 6)
tvar      <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:6, each=5)
uuu   <- rep(rnorm(6), each=5)
yvar  <- 1 + 2*tvar + uuu + rnorm(length(tvar))
simdat <- data.frame(idvar, timeorder, tvar, yvar)
head(simdat,12)


###################################################
### code chunk number 4: geepack-manual.Rnw:109-111
###################################################
mod1 <- geeglm(yvar~tvar, id=idvar, data=simdat, corstr="ar1")
mod1


###################################################
### code chunk number 5: geepack-manual.Rnw:127-133
###################################################
set.seed(123)
## library(doBy)
simdatPerm <- simdat[sample(nrow(simdat)),]
## simdatPerm <- orderBy(~idvar, simdatPerm)
simdatPerm <- simdatPerm[order(simdatPerm$idvar),]
head(simdatPerm)


###################################################
### code chunk number 6: geepack-manual.Rnw:143-145
###################################################
mod2 <- geeglm(yvar~tvar, id=idvar, data=simdatPerm, corstr="ar1")
mod2


###################################################
### code chunk number 7: geepack-manual.Rnw:152-155
###################################################
## simdatPerm2 <- orderBy(~timeorder, data=simdat)
simdatPerm2 <- simdat[order(simdat$timeorder),]
geeglm(yvar~tvar, id=idvar, data=simdatPerm2, corstr="ar1")


###################################################
### code chunk number 8: geepack-manual.Rnw:162-166
###################################################
wav <- simdatPerm$timeorder
wav
mod3 <- geeglm(yvar~tvar, id=idvar, data=simdatPerm, corstr="ar1", waves=wav)
mod3


###################################################
### code chunk number 9: geepack-manual.Rnw:175-181
###################################################
cor.fixed <- matrix(c(1    , 0.5  , 0.25,  0.125, 0.125,
                      0.5  , 1    , 0.25,  0.125, 0.125,
                      0.25 , 0.25 , 1   ,  0.5  , 0.125,
                      0.125, 0.125, 0.5  , 1    , 0.125,
                      0.125, 0.125, 0.125, 0.125, 1     ), 5, 5)
cor.fixed


###################################################
### code chunk number 10: geepack-manual.Rnw:189-191
###################################################
zcor <- fixed2Zcor(cor.fixed, id=simdatPerm$idvar, waves=simdatPerm$timeorder)
zcor


###################################################
### code chunk number 11: geepack-manual.Rnw:200-202
###################################################
mod4 <- geeglm(yvar~tvar, id=idvar, data=simdatPerm, corstr="fixed", zcor=zcor)
mod4