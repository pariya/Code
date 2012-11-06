library(base)
library(tmvtnorm)
library(glasso)
library(simone)


# main function

gnet = function(y, n.rho=10, em.iter=100, c.em.iter=1, gibbs.iter=1000, n.cores=1){
	cat("starting with net\n")
	p = ncol(y)
	rho = seq(0.01,.9,length=n.rho)
	cutoffs = calculate.cutoffs(y)
	lower.upper = calculate.lower.upper(y,cutoffs)
	cat("lower.upper calculated\n")

	results = list(Theta=list(), loglik=NULL, rho=rho)
	if(n.cores==1){
		for (i in 1:n.rho){
				R.gl <- calculate.EM.new(i, y, rho, lower.upper, em.iter, c.em.iter, gibbs.iter)
				results$Theta[[i]] = R.gl$wi
				results$loglik = c(results$loglik,R.gl$loglik)
	}
   }else{
	eff.cores <- min(n.cores,n.rho)
	if(eff.cores != n.cores) warning("Reduced the number of cores to",eff.cores,"otherwise we'd crash some cores without work")
	cl <- makeCluster(rep("localhost", eff.cores))
	clusterEvalQ(cl, library(base)) # in each cluster we should introduce libraries and functions. We apply it by clusterEvalQ argument
	clusterEvalQ(cl, library(tmvtnorm))
	clusterEvalQ(cl, library(glasso))
	clusterEvalQ(cl, library(simone))
	clusterEvalQ(cl, library(Matrix))
	clusterEvalQ(cl, source("codes/gnetMultiCores.R"))
	clusterEvalQ(cl, source("codes/calculate.R.R"))
	clusterEvalQ(cl, source("codes/calculate.lower.upper.R"))
	clusterEvalQ(cl, source("codes/calculate.cutoffs.R"))
	allR.gl <- parLapply(cl, 1:n.rho, get("calculate.EM.new"), y=y, rho=rho, lower.upper=lower.upper, em.iter=em.iter, c.em.iter=c.em.iter, gibbs.iter=gibbs.iter) #get ?	
	for(i in 1:n.rho){
		results$Theta[[i]] <- allR.gl[[i]]$wi
		results$loglik		= c(results$loglik, allR.gl[[i]]$loglik)
	}
	stopCluster(cl)
  }
  save(results, file=paste("ResultsAfterIteration",c.em.iter,".Rdata",sep=""))
	return(results)
}


## calculation EM for each rho (i=1,...,n.rho)
calculate.EM = function(i, y, rho, lower.upper, em.iter=100, c.em.iter=1, gibbs.iter=1000 ){
	cat("start calculation Rho", i, "\n")
	p = ncol(y)
	sigma.ini = diag(rep(1,p))
	if(missing(lower.upper)) lower.upper <- calculate.lower.upper(y)
	#LOAD results from (c.em.iter-1)
	#DO 1 iteration
	#Save again
	for( j in 1:em.iter) {
		cat("start of em.iter", j, "/", em.iter, "\n")
		s <- proc.time()
		#Dont forget to save this here before we do a real analysis
    #save(list(j, y, sigma.ini, lower.upper, gibbs.iter),file=paste("_debug_rho_",i,"i_",j,".Rdata",sep=""))
		R 	<- calculate.R(y,sigma.ini,lower.upper,gibbs.iter)
		#cat("R completed\n")
		R.gl <- glasso(R, rho[i], penalize.diagonal=FALSE)
		#print(R.gl$w)
    if(det(R.gl$w) <= 0){
       cat("Glasso is messing up we need a positive defined matrix, fixing using nearPD\n")
       R.gl$w <- nearPD(R.gl$w,keepDiag=TRUE)$mat 
		}
    sigma.ini <- R.gl$w 
		cat("Iteration", j, "completed after", (proc.time()-s)[3], "seconds\n")
	}
	cat("Rho", i, "completed\n")
	return(R.gl)	
	
}

## calculation EM for each rho (i=1,...,n.rho)
calculate.EM.new = function(i, y, rho, lower.upper, em.iter=100, c.em.iter=1, gibbs.iter=1000 ){
	p = ncol(y)
	sigma.ini = diag(rep(1,p))
	if(missing(lower.upper)) lower.upper <- calculate.lower.upper(y)
	#LOAD results from (c.em.iter-1)
	#DO 1 iteration
	#Save again
	if((c.em.iter-1) != 0){
    cat("Loading in old iteration results !!!\n")
    load(paste("iter_rho",i,"_res",(c.em.iter-1),".Rdata",sep=""))
  }
	cat("start of em.iter", c.em.iter, "/", em.iter, "\n")
	s <- proc.time()
	#Dont forget to save this here before we do a real analysis
  #save(list(j, y, sigma.ini, lower.upper, gibbs.iter),file=paste("_debug_rho_",i,"i_",j,".Rdata",sep=""))
	R 	<- calculate.R(y,sigma.ini,lower.upper,gibbs.iter)
	R.gl <- glasso(R, rho[i], penalize.diagonal=FALSE)
	
  if(det(R.gl$w) <= 0){
     cat("Glasso is messing up we need a positive defined matrix, fixing using nearPD\n")
     R.gl$w <- nearPD(R.gl$w,keepDiag=TRUE)$mat 
	}
  sigma.ini <- R.gl$w 
  cat("Iteration", c.em.iter, "completed after", (proc.time()-s)[3], "seconds\n")
	save(sigma.ini, file=paste("iter_rho",i,"_res",(c.em.iter),".Rdata",sep=""))
	#if((c.em.iter-1) != 0) file.delete(paste("iter_rho",i,"_res",(c.em.iter-1),".Rdata",sep=""))
	return(R.gl)
}
