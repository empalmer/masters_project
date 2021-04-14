# My packages 
library(tidyverse)

# Their packages
library(geepack)
library(mvtnorm)
library(mipfp)
library(lme4)


# n represents the sample size 
# which is.... 
n <- 100
# Dim represents the list of taxonomic structures
# For example it is only 2 OTUs that have 1 correlation 
dim <- 2
# Get the structure of the correlation matrix to use 
# based off of 2 OTU correlations and 2 repeated measures
# Take the first of the output, which is in the form Gamma like in the paper
# Subtract 1 to make diagonals 0 instead of 1
cor.index <- cor2zcor(n,
                     dim,
                     par = 2,
                     corstr = "exchangeable",
                     otustr = "exchangeable")[[1]] -1

#M is ... 
# In this case 4
M <- dim(cor.index)[1]

# K is the largest entry i the cor.index 
# in this case 3 
# Represents the clusters? at least in the paper 
# or the number of correlations to estimate 
K <- max(cor.index)

# J is ... 
J <- 1000


# rho.B is the "real" values of the 3 correlations 
rho.B <- c(0.3, 0.3, 0)
temp = as.numeric(cor.index)
for (k in 1:K) {
    temp[temp == k] = rho.B[k]
}
cor.number.B = matrix(temp, nrow = M, ncol = M)
diag(cor.number.B) = 1
# The above creates the "real" correlation matrix 


# rho.N is another? real value for 3 correlations 
# *** Why are there two? 
rho.N = c(0.3, 0.3, 0)
temp = as.numeric(cor.index)
for (k in 1:K) {
    temp[temp == k] = rho.N[k]
}
cor.number.N = matrix(temp, nrow = M, ncol = M)
diag(cor.number.N) = 1
cor.number.N


# *** what is mod? 
# Mod is the number of models that will be compared 
# So we have dimensions to initialize models 
mod <- 7
# Two values for the "real" beta values 
# Still unsure of what B and N stand for 
# Oh for the two parts of the two part model??
beta.B.all <- c(0, 0.1, -0.1)
beta.N.all <- c(0, 0.05, -0.05)


## Create placeholder matrices for results 
## nrow is 9 because... ?? 
mean.beta.B = matrix(nrow = 9, ncol = mod)
mean.beta.N = matrix(nrow = 9, ncol = mod)
power = matrix(nrow = 9, ncol = mod)
quantiles.B = matrix(nrow = 9, ncol = 2 * mod)
quantiles.N = matrix(nrow = 9, ncol = 2 * mod)



# The outermost loop is s:1-3
s <- 1:3
# what is s and why is it 3?
# maybe because thats the number of correlations to estimate?? idk
# Or, i think its all possible combinations of betas? 
# Then for loop for t:1-3
t <- 1:3
#instead make combination s and t 
st <- expand.grid(s,t)
st


##### for every pair st, do the following: #####

# store s and t values for pretend going throguth the loop 
s <- 1
t <- 1

# initialize 3 matrices to store things 
# What is J?? 
beta.est.B = matrix(nrow = J, ncol = mod)
beta.est.N = matrix(nrow = J, ncol = mod)
p.value = matrix(nrow = J, ncol = mod)

## See function to get the joint probability 
jp <- get_jp_B(s,M, beta.B.all)

# See function to get the RE values?? 
# Still unsure what this part does 
# take the tth value for the N part of the model 
lRE <- get_jp_N(t,n,beta.N.all)
RE1 <- lRE[[1]]
RE2 <- lRE[[2]]

## More for loops 
# now for J in 1: J
j <- 1:J

## After done looping through the following , do these summaries
mean.beta.B[3 * (s - 1) + t, ] = colMeans(beta.est.B)
mean.beta.N[3 * (s - 1) + t, ] = colMeans(beta.est.N)
power[3 * (s - 1) + t, ] = colSums(p.value < 0.05) / J
for (i in 1:mod) {
    quantiles.B[3 * (s - 1) + t, 2 * (i - 1) + 1:2] = quantile(beta.est.B[, i], c(0.025, 0.975), na.rm =
                                                                   TRUE)
    quantiles.N[3 * (s - 1) + t, 2 * (i - 1) + 1:2] = quantile(beta.est.N[, i], c(0.025, 0.975), na.rm =
                                                                   TRUE)
}

## start looping, then do above 
## Start j loop  

# Is this X the predictors?? 
# is 100 long, and is -1 and 1
X = rbinom(n, 1, 0.5) * 2 - 1
mu = beta.N * X

# More initializing.... 
Y = matrix(nrow = M, ncol = n)
Z = matrix(nrow = M, ncol = n)
W = matrix(nrow = M, ncol = n)
W.exp = matrix(nrow = M, ncol = n)


## Simulate Y,Z,W
for (i in 1:n) {
    Y[, i] = RMultBinary(n = 1, mult.bin.dist = jp[[X[i] / 2 + 1.5]])$binary.sequences
    Z[, i] = rmvnorm(1, -3 + rep(mu[i], M), cor.number.N)
    
    Z[Y[, i] == 0, i] = NA
    W[Y[, i] == 0, i] = -5
    W[Y[, i] == 1, i] = Z[Y[, i] == 1, i]
    W.exp[Y[, i] == 0, i] = 0
    W.exp[Y[, i] == 1, i] = exp(W[Y[, i] == 1, i])
    W.exp[W.exp[, i] > 1, i] = 0.99
}


Y = as.numeric(Y)
Z = as.numeric(Z)
W = as.numeric(W)
W.exp = as.numeric(W.exp)

## Transform X, repeating it each M times
X = rep(X, each = M)
id = rep(c(1:n), each = M)


# calculate the zcor, which is the 2nd argument in cor2zcor
# this format is needed for the gee function 
zcor = cor2zcor(n, dim, 2, corstr = "exchangeable", otustr = "exchangeable")[[2]]
zcor = as.matrix(zcor)
n.cor = dim(zcor)[2]

# calculate this "reduced" z-cor
# Unsure what this is... 
zcor.reduced = data.frame()
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



# Now FINALLY the Model fitting 
# fit 1 is the 0 part of the papers method 
# fit 
fit1 = geeglm(
    Y ~ X,
    family = "binomial",
    id = id,
    corstr = "userdefined",
    zcor = zcor
)
fit2 = geeglm(
    Z ~ X,
    family = "gaussian",
    id = id,
    corstr = "userdefined",
    zcor = zcor.reduced
)
p.value[j, 1] = summary(fit1)$coefficients[2, 4]
beta.est.B[j, 1] = summary(fit1)$coefficients[2, 1]
beta.est.N[j, 1] = NA
p.value[j, 2] = summary(fit2)$coefficients[2, 4]
beta.est.B[j, 2] = NA
beta.est.N[j, 2] = summary(fit2)$coefficients[2, 1]
cauchy.test = 0.5 * tan((0.5 - p.value[j, 1]) * pi) + 0.5 * tan((0.5 - p.value[j, 2]) *
                                                                    pi)
p.value[j, 3] = 1 - pcauchy(cauchy.test)
beta.est.B[j, 3] = summary(fit1)$coefficients[2, 1]
beta.est.N[j, 3] = summary(fit2)$coefficients[2, 1]




fit4.B = glm(Y ~ X, family = "binomial")
fit4.N = lm(Z ~ X)
cauchy.test = 0.5 * tan((0.5 - summary(fit4.B)$coefficients[2, 4]) * pi) +
    0.5 * tan((0.5 - summary(fit4.N)$coefficients[2, 4]) * pi)
p.value[j, 4] = 1 - pcauchy(cauchy.test)
beta.est.B[j, 4] = summary(fit4.B)$coefficients[2, 1]
beta.est.N[j, 4] = summary(fit4.N)$coefficients[2, 1]


fit5 = geeglm(
    W ~ X,
    family = "gaussian",
    id = id,
    corstr = "userdefined",
    zcor = zcor
)
p.value[j, 5] = summary(fit5)$coefficients[2, 4]
beta.est.B[j, 5] = summary(fit5)$coefficients[2, 1]
beta.est.N[j, 5] = summary(fit5)$coefficients[2, 1]


fit6 = lm(W ~ X)
p.value[j, 6] = summary(fit6)$coefficients[2, 4]
beta.est.B[j, 6] = summary(fit6)$coefficients[2, 1]
beta.est.N[j, 6] = summary(fit6)$coefficients[2, 1]

fit7 = lmer(W ~ X + (1 | RE1) + (1 | RE2))
F.summary = as.numeric(anova(fit7)[1, 4])
p.value[j, 7] = 1 - pchisq(F.summary, 1)
beta.est.B[j, 7] = summary(fit7)$coefficients[2, 1]
beta.est.N[j, 7] = summary(fit7)$coefficients[2, 1]