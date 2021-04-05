#' Simulate Allele Frequencies
#'
#' This function simulates allele frequencies of a desired population size under HWE
#' @param n desired populations size (set this value same as your actual population size for an accurate simulation)
#' @param nrun number of simulations to run on each allele frequency. The higher this number, the closer the simulations will be to the theoretical values (at the cost of computer power); 10000 is an optimal value.
#' @param res desired resolution of the theoretical allele frequency
#' @param plot whether to plot the simulation
#'
#' @return A list of two: allele_freqs: theoretical allele frequency, simulated_freqs: simulated frequencies at different confidence intervals
#'
#' @examples
#' alleles <- sim.als(n=500,nruns=10000,res=0.001)
#'
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
  outs<-lapply(pops,sim,nrun=nrun)
  out<-data.frame(do.call(rbind,outs))
  colnames(out)<-c("p2","het","q2","p2_sd","het_sd","q2_sd","p2_95","het_95",
                   "q2_95","p2_05","het_05","q2_05","p2_975","het_975","q2_975",
                   "p2_025","het_025","q2_025")
  if(plot){
    plot(out2$het~out2$q2,type="n",ylim=c(0,1),xlab="Proportion of homozygote alternative",ylab="Proportion of heterozygote")
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
#' @param cov.len Max value of depth of coverage to be simulated
#' @param sam.len Maximum no. of samples to be simulated
#' @param incr A vector of two integers indicating ncrement size for both depth and no. samples ranges
#'
#' @return A matrix of median allele rations where rows are the number of samples and columns are depth of coverage values
#'
#' @examples
#' depthVsSample(cov.len=20,sam.len=100)
#'
#' @export
depthVsSample<-function(cov.len=400,sam.len=1000,incr=c(1,1)){
  cov<- seq(1,cov.len,incr[1])
  nsamp<- seq(1,sam.len,incr[2])
  MR<-NULL
  for(i in cov){
    tmp<-dp.cov(cov.i=i,nsamp)
    MR <- cbind(MR,tmp)
  }
  colnames(MR)<-cov
  rownames(MR)<-nsamp
  return(MR)
  #saveRDS(MR,paste0(dirout,"/sim_SampVsDepth_",i,".rds"),compress = "gzip")
}

#' Identify significantly different heterozygotes from SNPs data
#'
#' This function will recognize the SNPs that are significantly different from the expected under HWE and plot potential duplicated snips
#'
#' @param dup.info duplication info table generated from filtered vcfs using the function dup.snp.info
#' @param plot logical. Wheather to plot the identified duplicated snps with the expected values
#'
#' @examples
#' duplicates<-sig.hets(dup.info,plot=T)
#'
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






