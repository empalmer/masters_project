
## A bunch of functions to make 


get_jp_B <- function(s,M,beta.B.all){
    # choose the true value of beta to simulate 
    # take the sth value
    beta.B = beta.B.all[s]
    # probabilities based off of values 
    # make probabilities for - and + beta.B? Why? 
    p1 = exp(-beta.B) / (1 + exp(-beta.B))
    p2 = exp(beta.B) / (1 + exp(beta.B))
    # From package mipfp 
    # Simulate joint distribution of K multivariate binary variables
    # repeat M times for ... 
    jp = list(
        ObtainMultBinaryDist(corr = cor.number.B, marg.probs = rep(p1, M)),
        ObtainMultBinaryDist(corr = cor.number.B, marg.probs = rep(p2, M))
    )
    # This gives an error, not sure why...
    
    return(jp)
}


assign_N_vals <- function(t, n, beta.N.all){
    beta.N = beta.N.all[t]
    # what is RE? 
    RE1 = rep(c(1:n) * 2 - 2, each = 4) + rep(rep(c(1:2), each = 2), n)
    RE2 = rep(c(1:n) * 2 - 2, each = 4) + rep(rep(c(1:2), 2), n)
    
    return(list(RE1,RE2))
}