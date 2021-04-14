#This function produces correlation structure matrix (first object in return list) and corresponding userdefined working correlation zcor requested by geepack (second object in return list). 
#arguments:
#n: sample size 
#dim: a list of taxonomic structure for OTU correlations. For example, Figure 1 in our paper is list(6,c(3,3),c(2,1,3),c(1,1,1,1,1,1)). The last vector of 1s can be omitted.
#par: other correlations. For example, c(2,2) for two twins and two correlated timepoints.
#corstr: correlation structure for other correlations, including exchangeable, independence and unstructured. Default is exchangeable.
#otustr: correlation structure for OTU correlations, including exchangeable and independence. Default is exchangeable.

library(Matrix)
cor2zcor = function(n,
                    dim,
                    par,
                    corstr = "exchangeable",
                    otustr = "exchangeable") {
    if (corstr == "independence" & otustr == "independence") {
        stop("zcor has zero dimension")
    }
    K = length(dim)
    N = sum(dim[[K]])
    count = dim
    if (K > 1) {
        for (k in (K - 1):1) {
            J = length(dim[[k]])
            count[[k]] = rep(1, J)
            for (j in 1:J) {
                index = 0
                while (sum(dim[[k + 1]][sum(count[[k]][1:j]):(sum(count[[k]][1:j]) + index)]) <=
                       dim[[k]][j]) {
                    index = index + 1
                    if (sum(dim[[k + 1]][sum(count[[k]][1:j]):(sum(count[[k]][1:j]) + index -
                                                               1)]) == dim[[k]][j])
                        break
                }
                count[[k]][j] = index
            }
        }
    }
    
    
    index.level = rep(1, K)
    for (k in K:1) {
        J = length(dim[[k]])
        squares = list()
        index.taxa = 0
        for (j in 1:J) {
            if (count[[k]][j] > 1) {
                index.taxa = index.taxa + 1
            }
            squares.new = list(matrix(
                index.level[k] + index.taxa,
                nrow = dim[[k]][j],
                ncol = dim[[k]][j]
            ))
            squares = c(squares, squares.new)
        }
        
        BD = as.matrix(bdiag(squares))
        diag(BD) = 1
        if (k < K) {
            if (sum(dim[[k]]) != N) {
                stop("Inconsistent dimensions: total number of OTUs unequal at each level")
            }
            for (s in 1:N) {
                for (t in 1:N) {
                    if (structure.OTU[s, t] > 0 & structure.OTU[s, t] < BD[s, t]) {
                        BD[s, t] = structure.OTU[s, t]
                    }
                }
            }
        }
        if (k > 1) {
            index.level[k - 1] = max(BD)
        }
        structure.OTU = BD
    }
    
    
    
    L = length(par)
    M = prod(par)
    block = as.matrix(1)
    for (l in 1:L) {
        structure = diag(par[l])
        if (corstr == "unstructured") {
            structure[lower.tri(structure)] = c(2:(1 + choose(par[l], 2)))
            structure = structure + t(structure) - diag(par[l])
        }
        else if (corstr == "exchangeable") {
            structure[lower.tri(structure)] = rep(2, choose(par[l], 2))
            structure = structure + t(structure) - diag(par[l])
        }
        else if (corstr == "independence") {
            structure = diag(par[l])
        }
        else{
            stop("Correlation structure undefined")
        }
        full = matrix(nrow = prod(par[1:l]), ncol = prod(par[1:l]))
        for (s in 1:par[l]) {
            for (t in 1:par[l]) {
                b = dim(block)[1]
                full[b * (s - 1) + 1:b, b * (t - 1) + 1:b] = max(block) * (structure[s, t] -
                                                                               1) + block
                if (corstr == "independence" & s != t) {
                    full[b * (s - 1) + 1:b, b * (t - 1) + 1:b] = matrix(0, nrow = b, ncol =
                                                                            b)
                }
            }
        }
        block = full
    }
    
    cor.matrix = matrix(nrow = N * M, ncol = N * M)
    for (s in 1:N) {
        for (t in 1:N) {
            cor.matrix[M * (s - 1) + 1:M, M * (t - 1) + 1:M] = max(full) * (structure.OTU[s, t] -
                                                                                1) + full
            if (s != t) {
                if (otustr == "independence") {
                    cor.matrix[M * (s - 1) + 1:M, M * (t - 1) + 1:M] = matrix(0, nrow = M, ncol =
                                                                                  M)
                }
                else if (corstr == "independence") {
                    cor.matrix[M * (s - 1) + 1:M, M * (t - 1) + 1:M] = structure.OTU[s, t] *
                        full
                }
            }
        }
    }
    
    # What does this code mean? 
    # N is dim of gamma matrix 
    # M is dim of omega matrix
    # This gives the output for the correlation matrix for geeglm
    zcor = matrix(nrow = n * choose(N * M, 2),
                  ncol = max(cor.matrix) - 1)
    for (i in 1:(max(cor.matrix) - 1)) {
        zcor[, i] = rep(as.numeric(cor.matrix[lower.tri(cor.matrix)] == i + 1), n)
    }
    zcor = zcor[, colSums(zcor != 0) > 0]
    return(list(cor.matrix, zcor))
}


# example
dim=list(6,c(3,3),c(2,1,3))
cor2zcor(1,dim,c(2,2),corstr="exchangeable",otustr="exchangeable")[[1]]
print(cor2zcor(1,dim,c(2,2),corstr="exchangeable",otustr="exchangeable")[[1]])
