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
#' This function will validate the detected duplicated-SNPs using a moving window approach (see details)
#'
#' @param d.detect a data frame of detected SNPs of duplicates/singletons (output of dup.detect)
#' @param window.size numerical. a single value of the desired moving window size (default=100 bp)
#'
#' @return a data frame of scaffold names and their average presence in the scaffold.
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' print("to be added")
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


