#
# gnet.R
#
# copyright (c) 2010-2012 - GBIC/Johannes Lab: Pariya Behrouzi, Danny Arends
# last modified Oct, 2012
# first written Oct, 2012
# 
# gnet
#

gnet = function(y, n.rho=10, em.iter=10, gibbs.iter=1000, n.cores=1){
  cat("Starting with Gnet\n")
  p = ncol(y)
  rho = seq(0.01,.9,length=n.rho)
  cutoffs = calculate.cutoffs(y)
  lower.upper = calculate.lower.upper(y, cutoffs)
  cat("Cut-offs calculated\n")
  results = list(Theta=list(), loglik=NULL, rho=rho)
  if(n.cores==1){
    for(i in 1:n.rho){ #PARRALEL  20 x NODES
      R.gl <- calculate.Theta(i,y,rho,lower.upper, em.iter, gibbs.iter)
      results$Theta[[i]] <- Matrix(R.gl$wi)
      results$loglik     = c(results$loglik,R.gl$loglik)
    }
  }else{
    eff.cores <- min(n.cores, n.rho)
    if(eff.cores != n.cores) warning("Reduced the number of cores to",eff.cores,"otherwise we'd crash some cores without work")
    cl <- makeCluster(rep("localhost",eff.cores))
    clusterEvalQ(cl, library(base))
    clusterEvalQ(cl, library(tmvtnorm))
    clusterEvalQ(cl, library(glasso))
    clusterEvalQ(cl, library(simone))
    clusterEvalQ(cl, library(Matrix))
    clusterEvalQ(cl, source("R/gnet.R"))
    clusterEvalQ(cl, source("R/calculate.cutoffs.R"))
    clusterEvalQ(cl, source("R/calculate.lower.upper.R"))
    clusterEvalQ(cl, source("R/calculate.R.R"))
    allR.gl <- parLapply(cl, 1:n.rho, get("calculate.Theta"), y=y, rho=rho, lower.upper=lower.upper, em.iter=em.iter, gibbs.iter=gibbs.iter)
    for(i in 1:n.rho){
      results$Theta[[i]] <- Matrix(allR.gl[[i]]$wi)
      results$loglik     = c(results$loglik,allR.gl[[i]]$loglik)
    }
    stopCluster(cl)
  }
  return(results)
}

calculate.Theta <- function(i, y, rho, lower.upper, em.iter = 10, gibbs.iter=1000){
  cat("Calculation of Rho", i,"\n")
  p = ncol(y)
  sigma.ini = diag(rep(1,p))
  for (j in 1:em.iter){
    s <- proc.time()
    cat("Start of iteration",j,"\n")
    R    <- calculate.R(y, sigma.ini, lower.upper, gibbs.iter)
    R.gl <- glasso(R, rho[i], penalize.diagonal=FALSE)
    sigma.ini <- R.gl$w
    cat("Iteration",j,"completed after",(proc.time()-s)[3],"seconds\n")
  }
  cat("Rho",i,"completed\n")
  return(R.gl)
}
