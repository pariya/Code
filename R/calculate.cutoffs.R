#
# calculate.cutoffs.R
#
# copyright (c) 2010-2012 - GBIC/Johannes Lab: Pariya Behrouzi, Danny Arends
# last modified Oct, 2012
# first written Oct, 2012
# 
# calculate.cutoffs
#

calculate.cutoffs = function(y){
  p<-ncol(y)
  n<-nrow(y)
  levels<-unique(sort(unlist(y)))
  n.levels<-length(levels)
  q<-matrix(nrow=p,ncol=n.levels)
  for(i in 1:p){
    X=factor(y[,i],level=levels)
    No<-tabulate(X)
    q[i,]<-qnorm(cumsum(No)/n)
  }
  q<-cbind(-Inf,q)
  return(q)
}
