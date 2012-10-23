#
# calculate.R.R
#
# copyright (c) 2010-2012 - GBIC/Johannes Lab: Pariya Behrouzi, Danny Arends
# last modified Oct, 2012
# first written Oct, 2012
# 
# calculate.R
#

calculate.R = function(y,sigma,lower.upper,gibbs.iter=1000){
  S=0
  p <- ncol(y)
  n <- nrow(y)
  lower.upper<-calculate.lower.upper(y)
  cl <- startCluster(rep("localhost",20))
  items <- list(1:n)
  for(i in 1:n){
  
    cat("calc.R -> ",i,"/",n,"\n")
    s <- proc.time()
    z <- rtmvnorm(n=gibbs.iter, sigma=sigma, lower=lower.upper$lower[i,], upper=lower.upper$upper[i,], algorithm="gibbsR")
    cat("calculated Z in",(proc.time()-s)[3]," seconds\n")
    summed <- matrix(0, p, p)
    for(iteration in 1: gibbs.iter){ 
      cat("Gibbs:",iteration,"/",gibbs.iter,"\n")
      summed <- summed + (z[iteration,] %*% t(z[iteration,]))
    }
    cat("calculated S.hlp\n")
    S.hlp = summed / gibbs.iter
    
    #bb <- apply(z,1,function(zi){zi %*% t(zi)})
    #S.hlp =   #= apply(apply(z,1,function(zi){zi %*% t(zi)}),1,mean)
    S = S + matrix(S.hlp,ncol=p)
  }
  return(cov2cor(S))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
}
