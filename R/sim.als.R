## Helper functions ##
#1. generate simulations for one population
sim<-function(x,nrun,n){y<-replicate(nrun,{tt<-sample(x,n,replace=F);c(sum(tt==0)/n,sum(tt==1)/n,sum(tt==2)/n)})
c(rowMeans(y),apply(y,1,sd),apply(y,1,function(x)quantile(x,p=0.95,na.rm=T)),
  apply(y,1,function(x)quantile(x,p=0.05,na.rm=T)),apply(y,1,function(x)quantile(x,p=0.975,na.rm=T)),
  apply(y,1,function(x)quantile(x,p=0.025,na.rm=T)))}
#2. generate median allele ratios for a given number of samples for one depth value
dp.cov<-function(depth,sam,sims){
  dout<-lapply(sam,function(x,depth,sims){y<-replicate(n=sims,{
    Allele1<-rbinom(n = x,size = depth,prob=0.5) # Binomial sampling of number of reads supporting Allele1
    Dev<- ((depth/2) - Allele1) / depth # deviation from expectation of 0.5
    Devsum<-abs(mean((Dev)))
  });return(mean(abs(y)))},depth=depth,sims=sims)
  dout<-simplify2array(dout)
  return(dout)
}

#2. generate Z-score confidence for allele ratios for a given number of samples for one depth value
dp.covZ<-function(nsamp,bias,depth,sims,p){
  dout<-lapply(bias,function(y,nsamp,depth,sims,p){rp<-replicate(n=sims,{
    Allele1<-rbinom(n = nsamp,size = depth,prob=y)
    Zscore<- ((depth*p) - Allele1) / sqrt(depth*p*(1-p))
    Zsum<-sum(Zscore)
    pval<-pnorm(q = Zsum, mean = 0, sd = sqrt(nsamp),lower.tail = F)
  })
  tmp<-length(which(rp <= 0.05))/sims
  return(tmp)
  },depth=depth,sims=sims,nsamp=nsamp,p=p) #;return(length(which(rp <= 0.05))/sims)
  dout<-simplify2array(dout)
  dout<-cbind(nsamp,bias,dout)
  dout<-dout[which(dout[,3]>=0.95),]
  return(dout)
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

#5 plot depth vs samples
plot.svd <- function(MR,cols=c("#1C86EE", "#00BFFF", "#DAA520", "#FF0000")){
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
  MR<-MR[-1,-1]
  colfunc <- colorRampPalette(cols,bias=ifelse(any(dim(MR))>30,.2,1))
  qt<-quantile(MR,na.rm = T,p=seq(.1,1,ifelse(length(MR)<10000,1/length(MR),0.0001)))
  qt<-c(-Inf,qt,Inf)
  cl<-colfunc(length(qt))
  tcl<-cl[cut(MR,breaks=qt)]
  #tcl[is.na(tcl)]<-cl[1]
  mat<-matrix(tcl,nrow=nrow(MR),ncol=ncol(MR))
  mat<-apply(t(mat),2,rev)
  ld <- as.raster(mat,nrow=nrow(mat))
  colfunc <- colorRampPalette(cols)
  cl<-colfunc(length(qt))
  legend_image <- as.raster(matrix(rev(cl), ncol=1))
  layout(matrix(1:2,ncol=2), widths = c(3,1),heights = c(1,1))
  par(mar=c(4,4,3,0))
  plot(0,type="n",xlim = range(as.numeric(colnames(MR))),ylim=range(as.numeric(rownames(MR))),xlab="Number of samples",ylab="Median depth coverage",frame=F)
  rasterImage(ld,2,2,range(as.numeric(colnames(MR)))[2],range(as.numeric(rownames(MR)))[2])
  par(mar=c(3,0.5,7,1))
  plot(c(0,2),c(0,1.3),type = 'n', axes = F,xlab = '', ylab = '')#, main = 'legend title'
  text(x=0.6,y=1.1,labels = 'Deviation from 0.5',cex=0.7,font=2)
  text(x=0.5, y = c(0.1,1), labels = c("Low","High"),cex=0.7)
  rasterImage(legend_image, 0, 0.1, 0.2,1)
}



#' Simulate Allele Frequencies
#'
#' This function simulates allele frequencies of a desired population size
#'  under HWE
#' @param n desired populations size (set this value same as your actual
#'  population size for an accurate simulation)
#' @param nrun number of simulations to run on each allele frequency.
#' The higher this number, the closer the simulations will be to the
#' theoretical values (at the cost of computer power); 10000 is an optimal value.
#' @param res desired resolution of the theoretical allele frequency
#' @param plot logical. whether to plot the simulation
#'
#' @return A list of two matrices:
#' 1. allele_freqs: theoretical allele frequency
#' 2. simulated_freqs: simulated frequencies at different confidence intervals
#'
#' @author Piyal Karunarathne, Pascal Milesi
#'
#'
#' @examples
#' \dontrun{alleles <- sim.als(n=200,nrun=1000,res=0.001,plot=TRUE)}
#'
#' @importFrom stats fisher.test median quantile rbinom sd smooth.spline
#' @importFrom graphics legend lines
#' @export
sim.als<-function(n=500,nrun=10000,res=0.001,plot=TRUE){
  #n/(1-n(e^2))
  nsim <- round(abs(n/(1-n*(0.05^2))),0)
  p2<-seq(0,1,res)
  p <- sqrt(p2)
  q <- 1-p
  Hex <- 2*p*q
  q2 <- q^2
  control <- p2+Hex+q2
  oo <- cbind(p,q,p2,Hex,q2,0,control)
  o0<-oo[,3:5]*nsim
  pops <- apply(o0,1,function(x){unlist(mapply(rep,c(0,1,2),x))})
  outs<-lapply_pb(pops,sim,nrun=nrun,n=n)
  out<-data.frame(do.call(rbind,outs))
  colnames(out)<-c("p2","het","q2","p2_sd","het_sd","q2_sd","p2_95","het_95",
                   "q2_95","p2_05","het_05","q2_05","p2_975","het_975","q2_975",
                   "p2_025","het_025","q2_025")
  if(plot){
    plot(out$het~out$q2,type="n",ylim=c(0,1),xlab="Proportion of homozygote alternative",ylab="Proportion of heterozygote")
    lines(oo[,4]~oo[,5],col="red")
    lines(smooth.spline(y=c(out$het_95),x=c(c(c(out$q2_95+ c(out$het_95/2))/c(out$het_95+out$q2_95+out$p2_95))^2)),col="green",lty=2)
    lines(smooth.spline(y=c(out$het_05),x=c(c(c(out$q2_05+ c(out$het_05/2))/c(out$het_05+out$q2_05+out$p2_05))^2)),col="green",lty=2)
    lines(smooth.spline(y=c(out$het_975),x=c(c(c(out$q2_975+ c(out$het_975/2))/c(out$het_975+out$q2_975+out$p2_975))^2)),col="blue",lty=2)
    lines(smooth.spline(y=c(out$het_025),x=c(c(c(out$q2_025+ c(out$het_025/2))/c(out$het_025+out$q2_025+out$p2_025))^2)),col="blue",lty=2)
    legend("topright",lty=c(1,2,2),col=c("red","green","blue"),legend=c("Expected under HWE",
                                                                        "Expected at 95/5 % CI","Expected at 97.5/2.5 CI"),cex=0.7,bty="n")
  }
  return(list(allele_freqs=oo,simulated_feqs=out))
}


#' Simulate median allele ratios for varying no. of samples and depth coverage
#'
#' This function will simulate the expected median allele ratios under HWE
#' for given ranges of no. of samples and depth coverage values.
#' This is useful if you need to find the cutoff values of allele ratios for
#' different no. of samples and depth of coverage values in your data set.
#'
#' @param cov.len max value of depth of coverage to be simulated
#' @param sam.len maximum no. of samples to be simulated
#' @param nsims numerical. no. of simulations to be done for each combination of samples and depth
#' depth and no. samples ranges
#' @param plot logical. Whether to plot the output (a plot of no. samples
#' vs median depth of coverage colored by median allele ratios)
#' @param col character. Two colors to add to the gradient
#'
#' @return A matrix of median allele ratios where rows are the number of
#' samples and columns are depth of coverage values
#'
#' @author Pascal Milesi, Piyal Karunarathne
#'
#' @examples
#' \dontrun{depthVsSample(cov.len=50,sam.len=100)}
#'
#' @importFrom stats fisher.test median quantile rbinom sd smooth.spline
#' @importFrom grDevices as.raster colorRampPalette
#' @importFrom graphics layout par rasterImage text
#' @export
depthVsSample<-function(cov.len=100,sam.len=100,nsims=1000,plot=TRUE,col=c("#1C86EE","#00BFFF","#DAA520","#FF0000")){
  sims<-nsims
  rdepth<-seq(2,cov.len,by=1)
  nsamp<-seq(2,sam.len,by=1)
  alsim<-matrix(NA,nrow=cov.len,ncol=sam.len)

  for(i in seq_along(rdepth)){
    pb <- txtProgressBar(min = 0, max = length(rdepth), style = 3, width = 50, char = "=")
    setTxtProgressBar(pb, i)
    tm<-dp.cov(depth=rdepth[i],sam=nsamp,sims=sims)
    alsim[i+1,]<-c(NA,tm)
  }
  close(pb)
  colnames(alsim)<-1:sam.len
  rownames(alsim)<-1:cov.len

  if(plot){
    plot.svd(MR=alsim,cols = col)
  }
  return(alsim)
}


#' Simulate and plot detection power of bias in allele ratios
#'
#' This function simulates 95% confidence level Z-score based detection power
#' of allele biases for a given number of samples and a range of depths
#'
#' @param Dlist numerical. vector of depths values to be tested
#' @param sam numerical. number of samples
#' @param intensity numerical. frequency of bias
#' @param nsims numerical. number of simulations to be done for each sample
#' @param p numerical. expected allele ratio (0.5 for data with known
#' sequencing biases)
#' @param plot logical. plot the output
#'
#' @return Returns a list of detection probability values for the given range of
#' samples and depth
#'
#' @author Pascal Milesi, Piyal Karunarathne
#'
#' @importFrom grDevices hcl.colors
#'
#' @export
power.bias<-function(Dlist=c(2,4,8,16),sam=100,intensity=0.005,nsims=1000,p=0.5,plot=TRUE){
  bias<-seq(0,0.5,intensity)
  sims=nsims
  d<-list()
  for(i in seq_along(Dlist)){
    pb <- txtProgressBar(min = 0, max = length(Dlist), style = 3, width = 50, char = "=")
    setTxtProgressBar(pb, i)
    depth=Dlist[i]
    saml<-lapply(1:sam,FUN=dp.covZ,bias=bias,depth=depth,sims=sims,p=p)
    saml<-do.call(rbind,saml)
    d[[i]]<-saml[,-3]
  }
  close(pb)
  names(d)<-Dlist
  if(plot){
    cl<-hcl.colors(length(Dlist),"dark 3")
    plot(c(0,0.5)~c(0,sam),typ="n",ylab = "Simutlated Allele Ratio",xlab = "Number of heterozygotes", main=paste0("Detection power of bias in allele ratio \n for the expected ratio of ",p))
    for(j in seq_along(Dlist)){
      sub<-unlist(by(d[[j]][,2],d[[j]][,1],max,simplify = F))
      lines(smooth.spline(sub~as.numeric(names(sub)),df=20),lwd = 2,col=cl[j])
    }
    legend("bottomright",legend = Dlist,lwd=2,col=cl[1:length(Dlist)],title = "Depth",bty="n",cex=.8)
  }
  return(d)
}




