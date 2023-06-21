############ Jay's plot of CNV lines ########


cnvlines<-readRDS("/Users/piyalkaru/Desktop/UPP/manuscripts/rCNV/submits/MER/revision1/CNV_line.rds")

cnvll<-split.data.frame(cnvlines,f=cnvlines$CN)

par(mfrow=c(1,2))
par(mar=c(5,6,2,1))
plot(cnvlines$Deviation~cnvlines$pDup,type="n",xlab="Proportion of Duplicates (H)",ylab="Deviation from expected ratio (0.5)")

#cl<-rev(hcl.colors(length(cnvll),palette = "Blues 3"))
cl<-colorRampPalette(colors=c(1,2))(length(cnvll))
cltp<-rCNV:::makeTransparent(cl,alpha = .5)

for(i in length(cnvll):1){
  tw<-cnvll[[i]]
  tw_ng<-tw[tw$Deviation<=0,]
  tw_ps<-tw[tw$Deviation>=0,]

  ss_ng<-smooth.spline(x=tw_ng[,2],y=tw_ng[,1],spar = 2,all.knots = T)
  ss_ps<-smooth.spline(x=tw_ps[,2],y=tw_ps[,1],spar = 2,all.knots = T)

  polygon(c(ss_ng$x,rev(ss_ps$x)),c(ss_ng$y,rev(ss_ps$y)),col=cltp[i],border=cl[i])
}
legend(x=0,y=.4,pch=19,bty="n",col=cltp,border=cl,legend=c(unique(cnvlines$CN)),title="No. of CNVs",horiz = T,y.intersp = .7,x.intersp = .5,cex=.9)


## fig iib
ttlines<-readRDS("/Users/piyalkaru/Desktop/UPP/manuscripts/rCNV/submits/MER/revision1/fig2b.rds")
ttlines$NC<-factor(ttlines$NC,labels=cl)
ttlines$NC<-as.character(ttlines$NC)
ttsmooth<-smooth.spline(x=ttlines$pDup,y=ttlines$Deviation,spar=1)

par(mar=c(6,6,2,1))
plot(ttlines$Deviation~ttlines$pDup,col=ttlines$NC,pch=20,cex=.2,xlab="Proportion of duplicates (H)",ylab="Variance of number of alternative alleles \nn heterozygotes (p = 0.5)")

