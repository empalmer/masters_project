set.seed(2)
library(geepack)
library(mvtnorm)
library(mipfp)
library(lme4)




n = 1000
n = 100
dim = 2
cor.index = cor2zcor(n, dim, 2, corstr = "exchangeable", otustr = "exchangeable")[[1]] -
    1
M = dim(cor.index)[1]
K = max(cor.index)
J = 100

rho.B = c(0.3, 0.3, 0)
temp = as.numeric(cor.index)
for (k in 1:K) {
    temp[temp == k] = rho.B[k]
}
cor.number.B = matrix(temp, nrow = M, ncol = M)
diag(cor.number.B) = 1

rho.N = c(0.3, 0.3, 0)
temp = as.numeric(cor.index)
for (k in 1:K) {
    temp[temp == k] = rho.N[k]
}
cor.number.N = matrix(temp, nrow = M, ncol = M)
diag(cor.number.N) = 1

mod = 7
beta.B.all = c(0, 0.1, -0.1)
beta.N.all = c(0, 0.05, -0.05)
mean.beta.B = matrix(nrow = 9, ncol = mod)
mean.beta.N = matrix(nrow = 9, ncol = mod)
power = matrix(nrow = 9, ncol = mod)
quantiles.B = matrix(nrow = 9, ncol = 2 * mod)
quantiles.N = matrix(nrow = 9, ncol = 2 * mod)
for (s in 1:3) {
    for (t in 1:3) {
        beta.est.B = matrix(nrow = J, ncol = mod)
        beta.est.N = matrix(nrow = J, ncol = mod)
        p.value = matrix(nrow = J, ncol = mod)
        
        
        beta.B = beta.B.all[s]
        p1 = exp(-beta.B) / (1 + exp(-beta.B))
        p2 = exp(beta.B) / (1 + exp(beta.B))
        jp = list(
            ObtainMultBinaryDist(corr = cor.number.B, marg.probs = rep(p1, M)),
            ObtainMultBinaryDist(corr = cor.number.B, marg.probs = rep(p2, M))
        )
        beta.N = beta.N.all[t]
        RE1 = rep(c(1:n) * 2 - 2, each = 4) + rep(rep(c(1:2), each = 2), n)
        RE2 = rep(c(1:n) * 2 - 2, each = 4) + rep(rep(c(1:2), 2), n)
        for (j in 1:J) {
            X = rbinom(n, 1, 0.5) * 2 - 1
            mu = beta.N * X
            Y = matrix(nrow = M, ncol = n)
            Z = matrix(nrow = M, ncol = n)
            W = matrix(nrow = M, ncol = n)
            W.exp = matrix(nrow = M, ncol = n)
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
            X = rep(X, each = M)
            id = rep(c(1:n), each = M)
            
            zcor = cor2zcor(n, dim, 2, corstr = "exchangeable", otustr = "exchangeable")[[2]]
            zcor = as.matrix(zcor)
            n.cor = dim(zcor)[2]
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
        }
        mean.beta.B[3 * (s - 1) + t, ] = colMeans(beta.est.B)
        mean.beta.N[3 * (s - 1) + t, ] = colMeans(beta.est.N)
        power[3 * (s - 1) + t, ] = colSums(p.value < 0.05) / J
        for (i in 1:mod) {
            quantiles.B[3 * (s - 1) + t, 2 * (i - 1) + 1:2] = quantile(beta.est.B[, i], c(0.025, 0.975), na.rm =
                                                                           TRUE)
            quantiles.N[3 * (s - 1) + t, 2 * (i - 1) + 1:2] = quantile(beta.est.N[, i], c(0.025, 0.975), na.rm =
                                                                           TRUE)
        }
        
    }
}
