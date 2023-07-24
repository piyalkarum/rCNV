################################################################################
######### All the data and analysis included in the rCNV publication ###########
#-------------------------------------------------------------------------------

library(rCNV)

### Excess of heterozygotes plot
# simulate allele frequencies
tt<-sim.als(n = 500, nrun = 10000, res = 0.001, plot = FALSE)
cl<-hcl.colors(5,palette = "plasma")
png("exHet.sim.png",h=400,w=450,bg="transparent")
plot(tt$simulated_feqs$het~tt$allele_freqs[,1],xlab="q (allele frequency)",
     ylab="Proportion of heterozygotes",typ="l", ylim=c(0,0.7),col=4)
points(tt$simulated_feqs$het_975~tt$allele_freqs[,1],typ="l",col="tomato")
legend("topleft",legend = c("Expected under HWE","95% CI"),lty=1,
       col=c(4,"tomato"),bty = "n", cex=1)
dev.off()


################### Data analysis #############################################
#1. Chinook salmon (data from Larson et al. 2014)
htc1<-readRDS("salmon_AD.corr.outRem.rds") #Allele depth table generated with hetTgen
htn0<-readRDS("salmon_AD.nor.corr.outRem.rds") #Normalized allele depth with cpm.normal
# generating metrics table
# the output plot is not included in the paper
AI<-allele.info(htc1,x.norm = htn0,plot.allele.cov=T)

# get deviants from all the SNPs
dv<-dupGet(AI, test = c("z.05","chi.05"))

# filter putative CNVs from deviants with unsupervised clustering
cv<-cnv(AI,test = c("z.05","chi.05"),filter = "kmeans")

# Alternative allele frequency
pp<-(AI$NHet+((AI$NHomAlt)*2))/(2*AI$Nsamp)

# custom colors for the plots (also used in all the following plots)
cl<-rCNV:::makeTransparent(c(1,hcl.colors(5,palette = "dark 3")[-2]),alpha = 0.6)

# plot of Z-score and deviants
plot(-scale(AI$z.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,cex=0.4,
     col=cl[ifelse(AI$z.05<0.05/nrow(AI),1,2)],bty="n")
# plot of Chi-square and deviants
plot(-scale(AI$chi.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,cex=0.4,
     col=cl[ifelse(AI$chi.05<0.05/nrow(AI),3,2)],bty="n")
# plot of excess of heterozygotes and deviants
plot(pp~AI$propHet,pch=20,cex=0.4,col=cl[ifelse(AI$eH.pval < 0.05/nrow(AI) &
                                                  AI$eH.delta > 0 ,4,2)],bty="n")

# add all statistics on one plot with median ratio vs proportion of heterozygotes
clr<-rep(cl[2],nrow(AI))
clr[which(AI$z.05<0.05/nrow(AI))]<-cl[1]
clr[which(AI$chi.05<0.05/nrow(AI))]<-cl[3]
clr[which(AI$eH.pval < 0.05/nrow(AI) & AI$eH.delta > 0)]<-cl[4]
plot(AI$medRatio~AI$propHet,pch=20,cex=0.4,bty="n",col=clr)


### density distribution plot of read-depth coverage
tt<-rCNV:::apply_pb(htn0[,-c(1:4)],2,function(x){
  tmp<-data.frame(stringr::str_split_fixed(x,",",n=2))
  tm<-apply(tmp,1,function(y){sum(as.numeric(y))})
  return(tm)
})
bv<-data.frame(htn0[,1:4],tt)
sing<-bv[match(cv$ID[which(cv$dup.stat=="singlet")],bv$ID),]
dub<-bv[match(cv$ID[which(cv$dup.stat=="duplicated")],bv$ID),]
ss<-apply(sing[,-c(1:4)],1,median)
dd<-apply(dub[,-c(1:4)],1,median)
rr<-max(density(ss,adjust = 2)$y)
plot(0,typ="n",xlim=c(0,120),ylim=c(0,rr),bty="n")
polygon(density(ss,adjust = 2),col=cl[2],lwd=1.5)
polygon(density(dd,adjust = 2),col=cl[1],lwd=1.5)


################################################################################
#2. American lobster (data from Dorant et al. 2020)
htc<-readRDS("lobster_AD.corr.outRem.rds") #Allele depth table generated with hetTgen
htn0<-readRDS("lobster_AD.nor.corr.outRem.rds") #Normalized allele depth with cpm.normal

# generating allele metrics from allele depth
AI<-allele.info(htc,x.norm = htn0,plot.allele.cov=T)

# Flagging deviants with p=0.5 for expected allele frequency
dv<-dupGet(AI,test = c("z.05","chi.05"))
dv.all<-dupGet(AI,test=c("z.all","chi.all")) #to compare with the z. & chi.all for ddRDAseq+RAPTURE

# filtering putative CNVs with p=0.5 for expected allele frequency
# with unsupervised clustering
cv<-cnv(AI,test = c("z.05","chi.05"), filter = "kmeans")
cv.all<-cnv(AI,test=c("z.all","chi.all"),filter = "intersection") #to compare with the z. & chi.all for ddRDAseq+RAPTURE

# Alternative allele frequence across all samples for all SNPs
pp<-(AI$NHet+((AI$NHomAlt)*2))/(2*AI$Nsamp)

# plotting of Z-score and deviants
plot(-scale(AI$z.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,cex=0.4,
     col=cl[ifelse(AI$z.05<0.05/nrow(AI),4,2)],bty="n")
# plotting of chi-square and deviants
plot(-scale(AI$chi.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,cex=0.4,
     col=cl[ifelse(AI$chi.05<0.05/nrow(AI),4,2)],bty="n")
# plotting deviants from excess of heterozygotes
plot(pp~AI$propHet,pch=20,cex=0.4,
     col=cl[ifelse(AI$eH.pval < 0.05/nrow(AI) & AI$eH.delta > 0 ,4,2)],
     bty="n")
# combined plot using the in-build dup.plot function
dup.plot(cv,col=cl[c(4,2)],ylim=c(0,0.8))

### plotting density distribution plot of read-depth coverage
tt<-rCNV:::apply_pb(htn0[,-c(1:4)],2,function(x){ #total depth of loci
  tmp<-data.frame(stringr::str_split_fixed(x,",",n=2))
  tm<-apply(tmp,1,function(y){sum(as.numeric(y))})
  return(tm)
})
bv<-data.frame(htn0[,1:4],tt)
sing<-bv[match(cv$ID[which(cv$dup.stat=="singlet")],bv$ID),]
dub<-bv[match(cv$ID[which(cv$dup.stat=="duplicated")],bv$ID),]
ss<-apply(sing[,-c(1:4)],1,median)
dd<-apply(dub[,-c(1:4)],1,median)
rr<-max(density(ss,adjust = 2)$y)
yd<-max(dd)+10
plot(0,typ="n",xlim=c(0,yd),ylim=c(0,rr),bty="n")
polygon(density(ss,adjust = 2),col=cl[2],lwd=1.5)
polygon(density(dd,adjust = 2),col=cl[1],lwd=1.5)

### plot Dorant detection
#comparison
orig.lobster<-readRDS("/Users/piyalkarunarathne/Desktop/UPPSALA/R/rCNV/tst/lobster/Lob.Dorant.original.rds")

dorant.table<-data.frame(cbind(AI[match(orig.lobster$ID,AI$ID),],orig.lobster$dup.stat))
dorant.table<-na.omit(dorant.table)

duplicates<-dorant.table[dorant.table$orig.lobster.dup.stat=="duplicated",]
singletons<-dorant.table[dorant.table$orig.lobster.dup.stat=="singlet",]

plot(dorant.table$medRatio~dorant.table$propHet,type="n",xlab="Proportion of heterozygotes",ylab="Median ratio")
points(singletons$medRatio~singletons$propHet,pch=19,cex=.5)
points(duplicates$medRatio~duplicates$propHet,pch=19,cex=.5,col=2)

plot(dorant.table$medRatio~dorant.table$propHet,type="n",xlab="Proportion of heterozygotes",ylab="Median ratio")
points(duplicates$medRatio~duplicates$propHet,pch=19,cex=.5,col=2)
points(singletons$medRatio~singletons$propHet,pch=19,cex=.5)


################################################################################

#3. Norway Spruce (subset of data from Chen et al. 2019)

# allele metrics table generated with allele.info function
AI<-readRDS("PA_AI_outlRem.rds")
info<-readRDS("sprucestatsscaffold.rds") # sorted list of loci on P.abies genome

# Flagging deviants with p determined from probe-bias for expected allele frequency
# this is especially important for exome-capture data like Norway spruce dataset
# where probe-biases are possible
dv<-dupGet(AI,c("z.all","chi.all"),plot=T)

# filtering putative CNVs with p determined from probe-bias for expected allele
# frequency
cv<-cnv(AI,c("z.all","chi.het"),filter = "kmeans",plot=T)

# proportion of alternative allele frequency for plotting
pp<-(AI$NHet+((AI$NHomAlt)*2))/(2*AI$Nsamp)

# plotting of Z-score and deviants
plot(-scale(AI$z.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,
     cex=0.4,col=cl[ifelse(AI$z.05<0.05/nrow(AI),4,2)],bty="n",ylim=c(-10,10))
# plotting of chi-square and deviants
plot(-scale(AI$chi.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,
     cex=0.4,col=cl[ifelse(AI$chi.05<0.05/nrow(AI),4,2)],bty="n",ylim=c(-30,0))
# plotting of deviants with excess of heterozygotes
plot(pp~AI$propHet,pch=20,cex=0.4,
     col=cl[ifelse(AI$eH.pval < 0.05/nrow(AI) & AI$eH.delta > 0 ,4,2)],bty="n")
# combined plot
dup.plot(cv,col=cl[c(4,2)],ylim=c(0,0.8))

### density distribution plot of read-depth coverage
bv<-readRDS("PA_ADto_CorrNor_outlRem.rds") # total allele depth table

sing<-bv[match(paste0(cv$CHROM,".",cv$POS)[which(cv$dup.stat=="singlet")],
               paste0(bv$CHROM,".",bv$POS)),]
dub<-bv[match(paste0(cv$CHROM,".",cv$POS)[which(cv$dup.stat=="duplicated")],
              paste0(bv$CHROM,".",bv$POS)),]
ss<-apply(sing[,-c(1:4)],1,median)
dd<-apply(dub[,-c(1:4)],1,median)
rr<-max(density(ss,adjust = 2)$y)
plot(0,typ="n",xlim=c(0,120),ylim=c(0,rr),bty="n")
polygon(density(ss,adjust = 2),col=cl[2],lwd=1.5)
polygon(density(dd,adjust = 2),col=cl[1],lwd=1.5)


###### plotting the enrichment putative duplicates along the genome ###########
AI2<-AI[match(info$CHROM.POS,paste0(AI$CHROM,".",AI$POS)),]

dv<-dupGet(AI,test = c("z.05","chi.het"),plot = F)
cv<-cnv(AI,test = c("z.05","chi.het"),filter = "kmeans",plot = F)
# selecting only the positions for which the relative position on the genome
# is available
cv2<-cv[match(info$CHROM.POS,paste0(cv$CHROM,".",cv$POS)),]
cv.i<-cnv(AI,test = c("z.all","chi.all"),filter = "intersection",plot = F)
cv.i2<-cv.i[match(info$CHROM.POS,paste0(cv.i$CHROM,".",cv.i$POS)),]

bv<-bv[match(paste0(cv2$CHROM,".",cv2$POS),paste0(bv$CHROM,".",bv$POS)),-c(1:4)]
nb2<-rowMeans(bv) # mean allele depth per SNP
nb2[cv2$dup.stat=="singlet"]<-0

# generate a sliding window for the plot to visualize the enrichment of CNVs
wnd=500
mps<-NULL
for(j in 1:length(nb2)){
  pb <- txtProgressBar(min = 0, max = length(nb2), style = 3, width = 50, char = "=")
  setTxtProgressBar(pb, j)
  frg<-j
  if(j+wnd<=length(nb2)){
    mps[j]<-mean(nb2[j:(j+wnd)])
  } else {mps[j]<-mean(nb2[j:(j-wnd)])}
}
mps<-data.frame(as.numeric(mps),as.numeric(rownames(cv2)))
mps<-mps[order(mps[,2]),]
thr<-quantile(mps[,1],p=.75)
thr2<-quantile(mps[,1],p=.95)
mps[mps[,1]<thr,1]<-0
mps[mps[,1]>thr2,1]<-mps[mps[,1]>thr2,1]+thr2
dp2<-dp+quantile(mps[,1],p=.95)

plot(dp2~rownames(cv2),pch=19,cex=0.4,col=cl[factor(cv2$dup.stat)],bty="n",ylim = c(0,max(dp)))
lines(mps[,2],mps[,1],col="dodgerblue")
################################################################################

# 4. Anopheles gambiae Mosquito (data from MaleriaGen)
htc<-readRDS("Ag_AD.corr.40plus40.rds") # allele depth table
htn0<-readRDS("Ag_AD.nor.corr.40plus40.rds") # normalized allele depth
# generating allelel metrics
AI<-allele.info(htc,x.norm = htn0,plot.allele.cov = F)

# Flagging deviants with p=0.5 for expected allele frequency
dv<-dupGet(AI,c("z.05","chi.05"))

# filtering putative CNVs with p=0.5 for expected allele frequency
cv<-cnv(AI,c("z.05","chi.05"),filter = "kmeans")

# proportion of alternative allele frequency for plotting
pp<-(AI$NHet+((AI$NHomAlt)*2))/(2*AI$Nsamp)

# plotting of Z-score and deviants
plot(-scale(AI$z.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,
     cex=0.4,col=cl[ifelse(AI$z.05<0.05/nrow(AI),4,2)],bty="n")
# plotting of chi-square and deviants
plot(-scale(AI$chi.05.sum/AI$NHet,center = T)~AI$propHet,pch=20,
     cex=0.4,col=cl[ifelse(AI$chi.05<0.05/nrow(AI),4,2)],bty="n")
# plotting deviants with excess of heterozygotes
plot(pp~AI$propHet,pch=20,cex=0.4,
     col=cl[ifelse(AI$eH.pval < 0.05/nrow(AI) & AI$eH.delta > 0 ,4,2)],bty="n")
# combined plot
dup.plot(cv,col=cl[c(4,2)],ylim=c(0,0.8))

## density distribution plot of read-depth coverage
tt<-rCNV:::apply_pb(htn0[,-c(1:4)],2,function(x){ # total depth of loci
  tmp<-data.frame(stringr::str_split_fixed(x,",",n=2))
  tm<-apply(tmp,1,function(y){sum(as.numeric(y))})
  return(tm)
})
bv<-data.frame(htn0[,1:4],tt)
bv[is.na(bv)]<-0
sing<-bv[match(cv$ID[which(cv$dup.stat=="singlet")],bv$ID),]
dub<-bv[match(cv$ID[which(cv$dup.stat=="duplicated")],bv$ID),]
ss<-apply(sing[,-c(1:4)],1,max)
dd<-apply(dub[,-c(1:4)],1,max)
rr<-max(density(ss,adjust = 2,na.rm = T)$y)
plot(0,typ="n",xlim=c(20,100),ylim=c(0,rr),bty="n")
polygon(density(ss,adjust = 2,na.rm = T),col=cl[2],lwd=1.5)
polygon(density(dd,adjust = 2,na.rm = T),col=cl[1],lwd=1.5)

###### plotting the enrichment putative duplicates along the genome ###########
rmn<-apply(tt,1,mean)
nb<-bv[match(cv$POS,bv$POS),]
nbr<-apply(nb[,-c(1:4)],1,function(x)mean(x))

nb2<-nbr
nb2[cv$dup.stat=="singlet"]<-0

# generate a sliding window for the plot to visualize the enrichment of CNVs
wnd=100 # window size
mps<-NULL
for(j in 1:length(nb2)){
  pb <- txtProgressBar(min = 0, max = length(nb2), style = 3, width = 50, char = "=")
  setTxtProgressBar(pb, j)
  frg<-j
  if(j+wnd<=length(nb2)){
    mps[j]<-mean(nb2[j:(j+wnd)])
  } else {mps[j]<-mean(nb2[j:(j-wnd)])}
}

thr<-quantile(mps,p=.75)
thr2<-quantile(mps,p=.95)
mps[mps<thr]<-0
plot(nbr~nb$POS,pch=19,cex=0.4,col=cl[factor(cv$dup.stat)],bty="n",ylim = c(0,max(nbr)))
points(mps~nb$POS,col="dodgerblue",type = "l")
################################################################################

################ depth vs sample plot (Fig. S1) #############################
#----------------------------------------------------------------------------
out3<-read.table("/Users/piyalkarunarathne/Desktop/UPPSALA/R/rCNV/tst/simul.2.txt",h=F)
#plotting
colfunc <- colorRampPalette(c("dodgerblue2","dodgerblue2","lightgoldenrod","orange","red"))
color<-colfunc(nrow(out3))
tcl<-cut(log(out3[,3]),nrow(out3))

#plot(out3[,1]~out3[,2],col=color[tcl],pch=19, ylab = "Depth of coverage",xlab = "Number heterozygotes",las =1,xlim=c(0,100),ylim=c(0,100))
layout(matrix(1:2,ncol=2), widths = c(3,1),heights = c(1,1))
par(mar=c(4,4,3,0))
plot(out3[,1]~out3[,2],col=color[tcl],pch=19, ylab = "Depth of coverage",xlab = "Number heterozygotes",las =1,xlim=c(0,100),ylim=c(0,100),frame=F)
# legend
par(mar=c(3,0.5,7,1))

color<-colfunc(nrow(out3))
legend_image <- as.raster(matrix(rev(color), ncol=1))
plot(c(0,2),c(0,1.5),type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.6,y=1.1,labels = 'Deviation from 0.5',cex=0.7,font=2)
text(x=0.5, y = c(0.1,1), labels = c("Low","High"),cex=0.7)
rasterImage(legend_image, 0, 0.1, 0.2,1)



