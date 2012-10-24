#
# calculate.R.R
#
# copyright (c) 2010-2012 - GBIC/Johannes Lab: Pariya Behrouzi, Danny Arends
# last modified Oct, 2012
# first written Oct, 2012
# 
# calculate.R - function to calculate something
#

calculate.R = function(y, sigma, lower.upper, gibbs.iter=1000, verbose = FALSE){
  if(missing(y)) stop("argument \"y\" is missing, with no default")
  S <- 0
  p <- ncol(y)
  n <- nrow(y)
  if(missing(sigma)){
    warning("sigma is taken matrix I")
    sigma <- diag(rep(1,p))
  }
  if(missing(lower.upper)) lower.upper <- calculate.lower.upper(y)
  
  for(i in 1:n){
    if(verbose) cat("calc.R -> ",i,"/",n,"\n")
    s <- proc.time()
    z <- rtmvnorm(n=gibbs.iter, sigma=sigma, lower=lower.upper$lower[i,], upper=lower.upper$upper[i,], algorithm="gibbsR")
    if(verbose) cat("calculated Z in",(proc.time()-s)[3]," seconds\n")
    summed <- matrix(0, p, p)
    for(iteration in 1: gibbs.iter){ 
      if(verbose) cat("Gibbs:",iteration,"/",gibbs.iter,"\n")
      summed <- summed + (z[iteration,] %*% t(z[iteration,]))
    }
    if(verbose) cat("calculating new S \n")
    S = S + (summed / gibbs.iter)
  }
  return(cov2cor(S))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
}
