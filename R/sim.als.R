## Helper functions ##
#1. generate simulations for one population
sim<-function(x,nrun,n){y<-replicate(nrun,{tt<-sample(x,n,replace=F);c(sum(tt==0)/n,sum(tt==1)/n,sum(tt==2)/n)})
c(rowMeans(y),apply(y,1,sd),apply(y,1,function(x)quantile(x,p=0.95,na.rm=T)),
  apply(y,1,function(x)quantile(x,p=0.05,na.rm=T)),apply(y,1,function(x)quantile(x,p=0.975,na.rm=T)),
  apply(y,1,function(x)quantile(x,p=0.025,na.rm=T)))}
#2. generate median allele ratios for a given number of samples for one depth value
dp.cov<-function(cov.i,nsamp){
  if(cov.i==0){return(rep(NA,length(nsamp)))}
  if(cov.i==1){return(rep(1,length(nsamp)))} else {
    unlist(lapply(nsamp,function(z,cov.i){
      reads<-replicate(z,rbinom(cov.i,1,prob=0.5))
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

#5 plot depth vs samples
plot.svd <- function(MR,cols=c("red","cyan")){
  colfunc <- colorRampPalette(cols)
  cols<-makeTransparent(colfunc(10),alpha = 0.7)

  MR2 <- MR
  MR2[which(MR2>0.5)]<-1-MR2[which(MR2>0.5)]
  qtt <-quantile(MR2,p=0.05,na.rm = T)
  qtt1 <- quantile(MR2,p=0.01,na.rm = T)
  mn <- mean(MR2,na.rm = T)
  md <- median(MR2,na.rm = T)
  MR3 <- MR2
  MR3[which(MR2>=md)]<-cols[10]
  MR3[which(MR2<md)]<-cols[10]
  MR3[which(MR2<mn) ]<-cols[9]
  MR3[which(MR2<qtt) ]<-cols[2]
  MR3[which(MR2<=qtt1)]<-cols[1]

  mat<-apply(t(MR3),2,rev)#rotate matrix -90
  ld <- as.raster(mat,nrow=nrow(mat))
  legend_image <- as.raster(matrix(rev(cols), ncol=1))
  layout(matrix(1:2,ncol=2), widths = c(3,1),heights = c(1,1))
  par(mar=c(4,4,3,0))
  plot(0,type="n",xlim = range(range(as.numeric(rownames(MR)))),ylim=range(range(as.numeric(colnames(MR)))),xlab="Number of samples",ylab="Median depth coverage",frame=F)
  rasterImage(ld,0,0,1000,200)
  par(mar=c(3,0.5,7,1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')#, main = 'legend title'
  text(x=0.9, y = c(0.1,1), labels = c("low confidence","high confidence"),cex=0.6)
  rasterImage(legend_image, 0, 0.1, 0.2,1)
}



#' Simulate Allele Frequencies
#'
#' This function simulates allele frequencies of a desired population size under HWE
#' @param n desired populations size (set this value same as your actual population size for an accurate simulation)
#' @param nrun number of simulations to run on each allele frequency. The higher this number, the closer the simulations will be to the theoretical values (at the cost of computer power); 10000 is an optimal value.
#' @param res desired resolution of the theoretical allele frequency
#' @param plot whether to plot the simulation
#'
#' @return A list of two matrices: 1. allele_freqs: theoretical allele frequency, 2. simulated_freqs: simulated frequencies at different confidence intervals
#'
#' @author Piyal Karunarathne, Pascal Milesi
#'
#' @references <add reference>
#'
#'
#' @examples
#' alleles <- sim.als(n=200,nrun=1000,res=0.001)
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
  outs<-lapply(pops,sim,nrun=nrun,n=n)
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
#' This function will simulate the expected median allele ratios under HWE for given ranges of no. of samples and depth coverage values.
#' This is useful if you need to find the cutoff values of allele ratios for different no. of samples and depth of coverage values in your dataset.
#'
#' @param cov.len max value of depth of coverage to be simulated
#' @param sam.len maximum no. of samples to be simulated
#' @param incr a vector of two integers indicating ncrement size for both depth and no. samples ranges
#' @param plot logical. Whether to plot the output (a plot of no. samles vs median depth coverage colored with median allele ratios)
#' @param plot.cols string. Two colors to add to the gradient
#'
#' @return A matrix of median allele rations where rows are the number of samples and columns are depth of coverage values
#'
#' @author Pascal Milesi, Piyal Karunarathne
#'
#' @references <add reference>
#'
#' @examples
#' depthVsSample(cov.len=50,sam.len=100)
#'
#' @importFrom stats fisher.test median quantile rbinom sd smooth.spline
#' @importFrom grDevices as.raster colorRampPalette
#' @importFrom graphics layout par rasterImage text
#' @export
depthVsSample<-function(cov.len=400,sam.len=1000,incr=c(1,1),plot=TRUE,plot.cols=c("red","cyan")){
  cov<- seq(1,cov.len,incr[1])
  nsamp<- seq(1,sam.len,incr[2])
  MR<-NULL
  for(i in cov){
    tmp<-dp.cov(cov.i=i,nsamp)
    MR <- cbind(MR,tmp)
  }
  colnames(MR)<-cov
  rownames(MR)<-nsamp
  if(plot){
    plot.svd(MR,cols=plot.cols)
  }
  return(MR)
  #saveRDS(MR,paste0(dirout,"/sim_SampVsDepth_",i,".rds"),compress = "gzip")
}



#' Identify significantly different heterozygotes from SNPs data
#'
#' This function will recognize the SNPs that are significantly different from the expected under HWE and plot potential duplicated snps
#'
#' @param dup.info duplication info table generated from filtered vcfs using the function dup.snp.info
#' @param plot logical. Whether to plot the identified duplicated snps with the expected values
#'
#' @return A matrix of expected heterozygote proportions from the observed data with p-value indicating significantly deviating snps, thus duplicates.
#' If enabled, the function also plots the recognized duplicates in red and expected (singletons) in black on a prop. of homozygote alt Vs. prop. of heterozygote plot
#' p2 observed homozygote reference
#' het expected
#' @examples
#' data(dup.info)
#' duplicates<-sig.hets(dup.info,plot=TRUE)
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats fisher.test median quantile rbinom sd smooth.spline
#' @export
sig.hets<-function(dup.info,plot=TRUE){
  d<-dup.info[,c("NHomFreq","NumHet","NHomRare","truNsample")]
  colnames(d)<-c("h1","het","h2","truNsample")
  df<-data.frame(t(apply(d,1,ex.prop)))
  colnames(df)<-c("p2","het","q2","pval","delta")
  if(plot){
    cols<-makeTransparent(c("red","black"),alpha=0.5)
    d$Color <- cols[2]
    d$Color [which(df$pval < 0.05 & df$delta > 0 )]<- cols[1]#& df$delta > 0
    plot(dup.info$PropHet~dup.info$PropHomRare, pch=19, cex=0.2,col=d$Color,
         xlab="Proportion of Alternate Homozygotes",ylab="Proportion of Heterozygotes")
    #lines(smooth.spline(y=out2$`1`,x=out2$`2`,spar = spar),col="green") # expected from HWE with simulation
    lines((smm<-smooth.spline(df$het~df$q2)),col="blue")
  }
  return(df)
}






