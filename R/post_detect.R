#helpers
# 1. moving window average
wind<-function(xx,dd){
  d<-dd[dd[,1]==xx,c(2,8)]
  xx<-unlist(xx)
  tmp<-cbind(min(d[,1]):max(d[,1]),d[match(min(d[,1]):max(d[,1]),d[,1]),2])
  tmp[is.na(tmp)]<-0
  tmp[,2]
}

#' Validate detected duplicates
#'
#' This function will validate the detected duplicated-SNPs using a moving
#' window approach (see details)
#'
#' @param d.detect a data frame of detected SNPs of duplicates and singlets
#'  (output of \code{dupGet})
#' @param window.size numerical. a single value of the desired moving window
#'  size (default \code{100} bp)
#'
#' @details Chromosome positions correctly ordered according to a reference
#' sequence is necessary for this function to work properly. Therefore, this
#' function is still in development for non-mapped reference sequences.
#'
#' @return A data frame of scaffold names and their average presence in the
#' scaffold.
#'
#' @author Piyal Karunarathne
#'
#' @export
dup.validate<-function(d.detect,window.size=100){
  nm<-unique(d.detect$Scaffold)
  gg<-lapply(nm,wind,dd=d.detect)
  names(gg)<-nm
  means<-lapply_pb(nm,function(x,mw,gg){yy<-unlist(gg[names(gg)==x])
  ll<-mw+(mw/2)
  if(length(yy)>ll){
    if(length(yy)>20000){
      y.list<-split(yy, ceiling(seq_along(yy)/10000))
      yl<-lapply(y.list,function(v,mw){
        vv<-unlist(v)
        if(length(vv)>ll){
          dpp<-NULL
          for(i in 1:(length(vv)-mw)){
            tmp<-unlist(vv)[i:(i+mw)]
            ss<-sum(tmp=="singleton")
            dd<-sum(tmp=="duplicated")
            dpp[i]<-dd/(dd+ss)
          }
          dpp[dpp==0]<-NA
          return(cbind(mean(dpp,na.rm = T),length(yy),sum(vv=="duplicated"),sum(vv=="singleton")))
        }
      },mw=mw)
      dp<-do.call(rbind,yl)
      dp<-cbind(x,dp)
    } else {
      dp<-NULL
      for(i in 1:(length(yy)-mw)){
        tmp<-unlist(yy)[i:(i+mw)]
        ss<-sum(tmp=="singleton")
        dd<-sum(tmp=="duplicated")
        dp[i]<-dd/(dd+ss)
      }
      dp[dp==0]<-NA
      dp<-mean(dp,na.rm = T)
      dp<-cbind(x,dp,length(yy),sum(yy=="duplicated"),sum(yy=="singleton"))
    }
  } else {
    dp<-sum(yy=="duplicated")/(sum(yy=="duplicated")+sum(yy=="singleton"))
    dp<-cbind(x,dp,length(yy),sum(yy=="duplicated"),sum(yy=="singleton"))
  }
  return(dp)
  },mw=window.size,gg=gg)
  dup.ratio<-do.call(rbind,means)
  dup.ratio[dup.ratio=="NaN"]<-0
  colnames(dup.ratio)<-c("Scaffold","dupl.ratio","scaf.length","no.duplicates","no.singletons")
  return(data.frame(dup.ratio))
}


#' Calculate population-wise Vst
#'
#' This function calculates Vst (variant fixation index) for populations given
#'  a list of duplicated loci
#'
#' @param AD data frame of total allele depth values of (duplicated, if
#' \code{id.list} is not provided) SNPs
#' @param pops character. A vector of population names for each individual.
#'  Must be the same length as the number of samples in AD
#' @param id.list character. A vector of duplicated SNP IDs. Must match the IDs
#'  in the AD data frame
#' @param qGraph logical. Plot the network plot based on Vst values
#' (see details)
#' @param \dots additional arguments passed to \code{qgraph}
#'
#' @importFrom qgraph qgraph
#' @importFrom grDevices boxplot.stats
#' @importFrom graphics barplot
#' @importFrom stats var
#' @importFrom utils combn
#'
#' @return Returns a matrix of pairwise Vst values for populations
#'
#' @details Vst is calculated with the following equation
#' \deqn{V_{T} = \frac{ V_{S} }{V_{T}}} where VT is the variance of normalized
#'  read depths among all individuals from the two populations and VS is the
#'  average of the variance within each population, weighed for population size
#'   (see reference for more details)
#' See \code{qgraph} help for details on qgraph output
#'
#' @references
#' Redon, Richard, et al. Global variation in copy number in the human genome.
#' nature 444.7118 (2006)
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{data(alleleINF)
#' data(ADtable)
#' DD<-dupGet(alleleINF)
#' ds<-DD[DD$dup.stat=="duplicated",]
#' ad<-ADtable[match(paste0(ds$CHROM,".",ds$POS),paste0(ADtable$CHROM,".",ADtable$POS)),]
#' vst(ad,pops=substr(colnames(ad)[-c(1:4)],1,11))}
#'
#' @export
vst<-function(AD,pops,id.list=NULL,qGraph=TRUE,...){
  if(!is.null(id.list)){
    AD<-AD[match(id.list,AD$ID),]
  }
  nm<-colnames(data.frame(AD))[-c(1:4)]
  pop<-na.omit(unique(pops))
  AD<-AD[,-c(1:4)]
  if(is.character(AD[1,1])){
    tm<-apply(AD,2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))))}))})
  }
  AD<-tm
  AD[AD==0]<-NA
  tmp<-data.frame(ind=nm,pop=pops,t(AD))
  # Vst - for CNVs
  # Vt-Vs/Vt
  # VT is the variance of normalized read depths among all individuals from the two populations and VS is the average of the variance within each population, weighed for population size
  Vst<-combn_pb(pop,2,function(x){
    jj<-tmp[tmp$pop==x[1],-c(1:2)]
    kk<-tmp[tmp$pop==x[2],-c(1:2)]
    ft<-ncol(jj)
    ll<-lapply(1:ft, function(y){
      vt<-var(c(jj[,y],kk[,y]),na.rm=T)
      vs<-(var(jj[,y],na.rm=T)*length(na.omit(jj[,y]))+
             var(kk[,y],na.rm=T)*length(na.omit(kk[,y])))/(length(na.omit(jj[,y]))+length(na.omit(kk[,y])))
      return((vt-vs)/vt)
    })
    vst<-mean(unlist(ll),na.rm = T)
    return(matrix(vst,dimnames=list(x[1],x[2])))
  },simplify = F)
  mt<-matrix(NA,nrow=length(pop),ncol=length(pop))
  dimnames(mt)<-list(pop,pop)
  for(i in 1:length(Vst)){
    mt[colnames(Vst[[i]]),rownames(Vst[[i]])]<-Vst[[i]]
  }
  if(qGraph){
    qgraph::qgraph(1/mt,layout="spring", ...=...)
  }
  return(mt)
}




