### helper functions ####
#1. generate simulations for one population
sim<-function(x,nrun){y<-replicate(nrun,{tt<-sample(x,n,replace=F);c(sum(tt==0)/n,sum(tt==1)/n,sum(tt==2)/n)})
c(rowMeans(y),apply(y,1,sd),apply(y,1,function(x)quantile(x,p=0.95,na.rm=T)),
  apply(y,1,function(x)quantile(x,p=0.05,na.rm=T)),apply(y,1,function(x)quantile(x,p=0.975,na.rm=T)),
  apply(y,1,function(x)quantile(x,p=0.025,na.rm=T)))}

#2. generate median allele ratios for a given number of samples for one depth value
dp.cov<-function(cov.i,nsamp){
  if(cov.i==0){return(rep(NA,length(nsamp)))}
  if(cov.i==1){return(rep(1,length(nsamp)))} else {
    unlist(lapply(nsamp,function(z,cov.i){
      reads<-replicate(z,rbinom(cov.i,1,p=0.5))
      tt<-apply(reads,2,table)
      if(is.list(tt)){
        md<-median(do.call(cbind,tt)[1,]/cov.i)
      } else if(is.matrix(tt)){
        md<- median(proportions(apply(reads,2,table),2)[1,],na.rm = T)
      } else {
        md<-NA
      }
      return(md)
    },cov.i))
  }
}

#3 make a given vector of colors transparent to a desired opacity
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}

#4 calculate expected proportions of heterozygotes from real data
ex.prop<-function(rs){n<-unlist(c(rs[4]))
p<-(rs[1]+rs[2]/2)/n
ob<-unlist(c(rs[1:3]))
eX<-unlist(c((p^2) * n,2*p*(1-p) * n,((1-p)^2) * n))
delta <- rs[2]-(2*p*(1-p) * n)
pval<-suppressWarnings(fisher.test(cbind(ob,eX)))$p.value
return(c(eX/n,pval,delta))}

