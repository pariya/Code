Code
====

By Pariya
<<<<<<HEAD
NEW
========
OLD
>>>>>>>>

================

Pariya Behrouzi 
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
