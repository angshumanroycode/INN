innc<-function(Xlist=NULL,Dlist=NULL){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,dist)
  Dlist=lapply(Dlist,as.matrix)
  ns=as.numeric(sapply(Dlist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  result=INNC(Dlist,n,length(Dlist))
  names(result)=c("TCsum","TCmax")
  return(result)
}

inns<-function(Xlist=NULL,Dlist=NULL){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,dist)
  Dlist=lapply(Dlist,as.matrix)
  ns=as.numeric(sapply(Dlist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  result=INNS(Dlist,n,length(Dlist))
  names(result)=c("TSsum","TSmax")
  return(result)
}

inn2<-function(Xlist=NULL,Dlist=NULL){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,dist)
  Dlist=lapply(Dlist,as.matrix)
  ns=as.numeric(sapply(Dlist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  d=length(Dlist)
  Dlist.sq=lapply(Dlist,function(X) X^2)
  Dlist.sq.total=list(Reduce("+",Dlist.sq))[rep(1,d)]
  Dlist.sq.rest=mapply("-",Dlist.sq.total,Dlist.sq,SIMPLIFY=FALSE)
  result=INN2(Dlist.sq,Dlist.sq.rest,n,d)
  names(result)=c("T2sum","T2max")
  return(result)
}

inn2c<-function(Dlist.sq,n,d){
  Dlist.sq.total=list(Reduce("+",Dlist.sq))[rep(1,d)]
  Dlist.sq.rest=mapply("-",Dlist.sq.total,Dlist.sq,SIMPLIFY=FALSE)
  return(INN2(Dlist.sq,Dlist.sq.rest,n,d))
}

innc.test<-function(Xlist=NULL,Dlist=NULL,B=100,alpha=0.05){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,dist)
  Dlist=lapply(Dlist,as.matrix)
  ns=as.numeric(sapply(Dlist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  result0=INNCtest(Dlist,n,length(Dlist),B)
  permval=result0[[1]]
  val=result0[[2]]
  pval=result0[[3]]
  cutoff=apply(permval,2,function(x) as.numeric(quantile(x,1-alpha)))
  result=list(val[1],val[2],cutoff[1],cutoff[2],pval[1],pval[2])
  names(result)=c("TCsum.stat","TCmax.stat","TCsum.cutoff","TCmax.cutoff",
                  "TCsum.pvalue","TCmax.pvalue")
  return(result)
}

inns.test<-function(Xlist=NULL,Dlist=NULL,B=100,alpha=0.05){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,dist)
  Dlist=lapply(Dlist,as.matrix)
  ns=as.numeric(sapply(Dlist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  result0=INNStest(Dlist,n,length(Dlist),B)
  permval=result0[[1]]
  val=result0[[2]]
  pval=result0[[3]]
  cutoff=apply(permval,2,function(x) as.numeric(quantile(x,1-alpha)))
  result=list(val[1],val[2],cutoff[1],cutoff[2],pval[1],pval[2])
  names(result)=c("TSsum.stat","TSmax.stat","TSsum.cutoff","TSmax.cutoff",
                  "TSsum.pvalue","TSmax.pvalue")
  return(result)
}

inn2.test<-function(Xlist=NULL,Dlist=NULL,B=100,alpha=0.05){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,dist)
  Dlist.sq=lapply(Dlist,function(X) (as.matrix(X))^2)
  ns=as.numeric(sapply(Dlist.sq,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  result0=INN2test(Dlist.sq,n,length(Dlist.sq),B)
  permval=result0[[1]]
  val=result0[[2]]
  pval=result0[[3]]
  cutoff=apply(permval,2,function(x) as.numeric(quantile(x,1-alpha)))
  result=list(val[1],val[2],cutoff[1],cutoff[2],pval[1],pval[2])
  names(result)=c("T2sum.stat","T2max.stat","T2sum.cutoff","T2max.cutoff",
                  "T2sum.pvalue","T2max.pvalue")
  return(result)
}