################################################################################
################# depth vs samples simulation ##########################
nsimul<- 1000 # Number of simulations
var<-seq(2,200,1) # vector of paramater space (e.g. for all depth value from 2 to 50 with a step of 2)
#out2<-NULL
out3<-matrix(NA,length(var)^2,3)
k<-1
for(j in var){ # j will take in turn all values in var
  for( i in var){ # i will take in turn all values in var
    out<-replicate(nsimul,{ # replicates what is into brackets nsimul time
      nsamp<-i # number of samples
      depth<-j # Depth of coverage (Number of reads covering a given position)
      Allele1<-rbinom(n = nsamp,size = depth,p=0.5) # Binomial sampling of number of reads supporting Allele1
      Dev<- ((depth/2) - Allele1) / depth # deviation from expectation of 0.5
      # Zscore<- ((depth/2) - Allele1) / sqrt(depth*0.25) #
      Devsum<-abs(mean((Dev))) # average absolute deviation for nsamp individuals
      #pval<-pnorm(q = Zsum, mean = 0, sd = sqrt(nsamp))
      return((Devsum))}
    )
    out3[k,1]<-j
    out3[k,2]<-i
    out3[k,3]<-mean(abs(out)) #Average across the number of simulation
    #out2[k]<-mean(abs(out))
    k<-k+1 # k is called a pointer, in that case it allow filling the matrix in the right row
  }
}

#plotting
colfunc <- colorRampPalette(c("dodgerblue2","lightgoldenrod","orange","orangered","orangered","red","red"))
color<-colfunc(nrow(out3))
tcl<-cut(out3[,3],nrow(out3))
#plot(out3[,1]~out3[,2],col=color[tcl],pch=19, ylab = "Depth of coverage",xlab = "Number heterozygotes",las =1,xlim=c(0,100),ylim=c(0,100))
layout(matrix(1:2,ncol=2), widths = c(3,1),heights = c(1,1))
par(mar=c(4,4,3,0))
plot(out3[,1]~out3[,2],col=color[tcl],pch=19, ylab = "Depth of coverage",xlab = "Number heterozygotes",las =1,xlim=c(0,100),ylim=c(0,100),frame=F)
# legend
par(mar=c(3,0.5,7,1))
colfunc <- colorRampPalette(c("dodgerblue2","lightgoldenrod","orange","red"))
color<-colfunc(nrow(out3))
legend_image <- as.raster(matrix(rev(color), ncol=1))
plot(c(0,2),c(0,1.5),type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.6,y=1.1,labels = 'Deviation from 0.5',cex=0.7,font=2)
text(x=0.5, y = c(0.1,1), labels = c("Low","High"),cex=0.7)
rasterImage(legend_image, 0, 0.1, 0.2,1)


#### confidence level #########
fun<-function(x){ # create a function with all the stuffs inside
  nsimul<- 1000 # Number of simulations
  var1<-seq(1,100,2) # n sample
  var2<-seq(0,0.5,0.005) # intensity of bias (from 0 to 0.5 step = 0.01)
  out3<-matrix(NA,length(var1)*length(var2),3) # output matrix, full of NA so far
  k<-1 # pointer
  for(j in var1){ # j will take in turn all values in var1
    for( i in var2){ # i will take in turn all values in var2
      out<-replicate(nsimul,{ # replicates what is into brackets nsimul time
        nsamp<-j # number of samples
        depth<-x # Depth of coverage (Number of reads covering a given position);

        # When i ≠ 0.5 you simulate a bias from the 50% - 50% expectation.
        Allele1<-rbinom(n = nsamp,size = depth,p=i) # Binomial sampling of number of reads supporting Allele1.

        p <- 0.5 # I write it like this for you to have the correct formalization
        Zscore<- ((depth*p) - Allele1) / sqrt(depth*p*(1-p)) # "True Z-score: (mean(x) - xi) / sd (x); it follows a normal distribution with mean = 0 and standard deviation = 1

        Zsum<-sum(Zscore) # Sum of the deviation for nsamp individuals. It follows a normal distribution of mean = sum of the mean of the individual distributiona (nsamp *0) and variance = sum of the individuals variance (nsamp * 1)

        pval<-pnorm(q = Zsum, mean = 0, sd = sqrt(nsamp),lower.tail = F) # test for the significance of the deviation (this is only possible because we can predict the distribution of Zsum, its mean and its standard deviation : sqrt(variance))

        return((pval))} # out is now a vector of nsimul pvalues.
      )
      out3[k,1]<-j
      out3[k,2]<-i
      out3[k,3]<-length(which(out <= 0.05))/nsimul # how many pvalues are < 0.05 over the 2000 simulations
      k<-k+1 # k is called a pointer, in that case it allow filling the matrix in the right row
    }
  }
  outdepth<-out3[which(round(out3[,3],2)>=0.95),] # subset of values of Allele ratio and nsamp for which the fraction of significant pvalues is > 0.95 and < 0.96; will be used for plotting.
  return(outdepth)
}
test<-sapply(c(4,8,16,32,64), function(x) fun(x)) # apply function "fun" with x (thus depth) varying from 2 to 100 by step 2.
plot(c(0,0.5)~c(0,100),col="white",ylab = "Simutlated Allele Ratio",xlab = "Number of heterozygotes", main="Detection power of bias in allele ratio")
for(l in 1:length(test)){
  sub<-unlist(by(test[[l]][,2],test[[l]][,1],max,simplify = F))
  #points(sub~as.numeric(names(sub)),col=l,pch=20,cex=0.4)
  #value<-loess(sub~as.numeric(names(sub)))
  #lines(predict(value)~as.numeric(names(sub)))
  lines(smooth.spline(sub~as.numeric(names(sub)),df=20),col=l,lty=l,lwd = 2)
}
legend("bottomright",legend = c(4,8,16,32,64),lwd=2,lty=c(1:length(test)),col=c(1:length(test)),title = "Depth",bty="n",cex=.8)


