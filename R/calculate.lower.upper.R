#
# calculate.lower.upper.R
#
# copyright (c) 2010-2012 - GBIC/Johannes Lab: Pariya Behrouzi, Danny Arends
# last modified Oct, 2012
# first written Oct, 2012
# 
# calculate.lower.upper
#

calculate.lower.upper = function(y,cutoffs){
  cutoffs<-calculate.cutoffs(y)
  levels<-unique(sort(unlist(y)))
  n.levels<-length(levels)
  n<-nrow(y)
  p<-ncol(y)
  lower = matrix(nrow=n,ncol=p)
  upper = matrix(nrow=n,ncol=p)
  for (i in 1:n){
    sel <- match(y[i,],levels)
    lower[i,]<-apply(cbind(sel,cutoffs),1,function(x){x[x[1]+1]})
    upper[i,]<-apply(cbind(sel,cutoffs),1,function(x){x[x[1]+2]})
  }
  return(list(lower=lower,upper=upper))
}
