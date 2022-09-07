## old function
allele.freq1<-function(gtt,verbose=TRUE){
  gs<-gtt[,-c(1:4)]
  if(verbose){
    tmp<-apply_pb(gs,1,function(x){
      x<-as.character(x)
      tl<-strsplit(x,"/")
      tt<-unlist(lapply(tl,function(y){length(which(y==1))/2}))
      return(tt)
    })
  } else {
    tmp<-apply_pb(gs,1,function(x){
      x<-as.character(x)
      tl<-strsplit(x,"/")
      tt<-unlist(lapply(tl,function(y){length(which(y==1))/2}))
      return(tt)
    })
  }
  tmp<-t(tmp)
  tmp[gs=="./."]<-NaN
  colnames(tmp)<-colnames(gs)
  rownames(tmp)<-paste0(gtt[,1],".",gtt[,2])
  tmp<-data.frame(gtt[,1:4],tmp)
  return(tmp)
}


#' Generate allele frequency table for individuals or populations
#'
#' Get alternative allele frequency across all individuals per SNP from the
#'  genotype or allele depth tables
#'
#' @param gtt a list or data frame of genotype and/or allele depth table produced from \code{hetTgen} (or similar)
#' @param f.typ character. type of allele frequency to be calculated (individual \code{"ind"} or population \code{"pop"})
#' @param verbose logical. whether to show the progress of the analysis
#'
#' @details If the allele frequencies to be calculated for populations from both getnotype table and the allele depth table, they must be provided in a list with element names \code{AD} for allele depth table and \code{GT} for the genotype table. See the examples.
#'
#' @return Returns a data frame or a list (if both genotype and allele depth used)
#' of allele frequencies
#'
#' @author Piyal Karunarathne
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' het.table<-hetTgen(vcf,"GT")
#' ad.table<-hetTgen(vcf,"AD")
#'
#' # for individual based AF
#' frQ<-allele.freq(het.table,f.typ="ind")
#'
#' #for population-wise and both allele depth and genotype tables
#' \dontrun{frQ<-allele.freq(list(AD=ad.table,GT=het.table),f.typ="pop")}
#'
#' @export
allele.freq<-function(gtt,f.typ=c("pop","ind"),verbose=TRUE){
  ty<-match.arg(f.typ,several.ok = TRUE)
  if(length(ty)>1){stop("Please select one output type for f.typ=")}
  ### individuals ###
  if(ty=="ind"){
    if(!inherits(gtt,"list")){gs<-gtt[,-c(1:4)]} else {stop("Provide a data frame")}
    if(grepl("/||",gs[1,1],fixed = FALSE)){
      if(verbose){
        message("genotype table provided\ncalculating allele frequency from genotype")
        tmp<-apply_pb(gs,1,function(x){
          apply(stringr::str_split_fixed(x,pattern=c("/||"),n=4),1,function(x)sum(x==1)/2)
        })
      } else {
        tmp<-apply(gs,1,function(x){
          apply(stringr::str_split_fixed(x,pattern=c("/||"),n=4),1,function(x)sum(x==1)/2)
        })
      }
      tmp<-t(tmp)
      tmp[gs=="./."]<-NaN
      colnames(tmp)<-colnames(gs)
      rownames(tmp)<-paste0(gtt[,1],".",gtt[,2])
      tmp<-data.frame(gtt[,1:4],tmp)
      return(tmp)
    } else if (grepl(",",gs[1,1],fixed = TRUE)){

      if(verbose){
        message("allele depth table provided\ncalculating allele frequency from allele depth")
        tmp<-apply_pb(gs,1,function(snp1){
          y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
          y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
          proportions(as.matrix(y),margin = 1)[,2]
        })
      } else {
        tmp<-apply(gs,1,function(snp1){
          y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
          y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
          proportions(as.matrix(y),margin = 1)[,2]
        })
      }
      tmp<-t(tmp)
      tmp[gs=="./."]<-NaN
      colnames(tmp)<-colnames(gs)
      rownames(tmp)<-paste0(gtt[,1],".",gtt[,2])
      tmp<-data.frame(gtt[,1:4],tmp)
      return(tmp)
    }
    ### populations ####
  } else if(ty=="pop"){
    ## for when gtt is provided as a list
    if(inherits(gtt,"list")){
      ### for when both GT and AD provided
      if(length(gtt)>1){
        ad<-gtt$AD[,-c(1:4)]
        gt<-gtt$GT[,-c(1:4)]
        ## for GT ##
        if(verbose){
          message("AF from genotype")
          tmp.g<-apply_pb(gt,1,function(x){
            tt<-stringr::str_split_fixed(x,pattern=c("/||"),n=4)
            sum(tt==1)/(sum(tt==1)+sum(tt==0))
          })
        } else {
          tmp.g<-apply(gt,1,function(x){
            tt<-stringr::str_split_fixed(x,pattern=c("/||"),n=4)
            sum(tt==1)/(sum(tt==1)+sum(tt==0))
          })
        }
        ## for AD ##
        if(verbose){
          message("AF from allele depth")
          tmp.a<-apply_pb(ad,1,function(snp1){
            y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
            y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            mean(proportions(as.matrix(y),margin = 1)[,2],na.rm = TRUE)
          })
        } else {
          tmp.a<-apply(ad,1,function(snp1){
            y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
            y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            mean(proportions(as.matrix(y),margin = 1)[,2],na.rm = TRUE)
          })
        }
        return(list(AD=tmp.a,GT=tmp.g))

        ### for when only one type of data provided
      } else if(length(gtt)==1){
        gs<-gtt[1][,-c(1:4)]
        if(grepl("/||",gs[1,1],fixed = FALSE)){
          ## for GT
          if(verbose){
            message("genotype table provided\ncalculating allele frequency from genotype")
            tmp<-apply_pb(gs,1,function(x){
              tt<-apply(stringr::str_split_fixed(x,pattern=c("/||"),n=4),1,function(x)sum(x==1)/2)
              return(mean(tt,na.rm = TRUE))
            })
          } else {
            tmp<-apply(gs,1,function(x){
              tt<-apply(stringr::str_split_fixed(x,pattern=c("/||"),n=4),1,function(x)sum(x==1)/2)
              return(mean(tt,na.rm = TRUE))
            })
          }
          return(tmp)
        } else if(grepl(",",gs[1,1],fixed = TRUE)){
          ## for AD
          if(verbose){
            message("allele depth table provided\ncalculating allele frequency from allele depth")
            tmp<-apply_pb(gs,1,function(snp1){
              y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
              y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
              mean(proportions(as.matrix(y),margin = 1)[,2],na.rm = TRUE)
            })
          } else {
            tmp<-apply(gs,1,function(snp1){
              y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
              y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
              mean(proportions(as.matrix(y),margin = 1)[,2],na.rm = TRUE)
            })
          }
          return(tmp)
        }
      }
    } else {
      ## for when gtt is provided as a data frame
      if(any(colnames(gtt)=="CHROM")){gs<-gtt[,-c(1:4)]} else {gs<-gtt}
      if(grepl("/||",gs[1,1],fixed = FALSE)){
        if(verbose){
          message("genotype table provided\ncalculating allele frequency from genotype")
          tmp<-apply_pb(gs,1,function(x){
            tt<-apply(stringr::str_split_fixed(x,pattern=c("/||"),n=4),1,function(x)sum(x==1)/2)
            return(mean(tt,na.rm = TRUE))
          })
        } else {
          tmp<-apply(gs,1,function(x){
            tt<-apply(stringr::str_split_fixed(x,pattern=c("/||"),n=4),1,function(x)sum(x==1)/2)
            return(mean(tt,na.rm = TRUE))
          })
        }
        return(tmp)
      } else if(grepl(",",gs[1,1],fixed = TRUE)){
        if(verbose){
          message("allele depth table provided\ncalculating allele frequency from allele depth")
          tmp<-apply_pb(gs,1,function(snp1){
            y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
            y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            mean(proportions(as.matrix(y),margin = 1)[,2],na.rm = TRUE)
          })
        } else {
          tmp<-apply(gs,1,function(snp1){
            y<-data.frame(stringr::str_split_fixed(snp1,pattern=",",n=2))
            y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            mean(proportions(as.matrix(y),margin = 1)[,2],na.rm = TRUE)
          })
        }
        return(as.data.frame(tmp))
      }
    }
  }
}


allele.freq0<-function(gtt,f.typ=c("pop","ind"),verbose=TRUE){
  ty<-match.arg(f.typ,several.ok = T)
  if(length(ty)>1){stop("Please select one output type for f.typ=")}
  if(ty=="ind"){
    if(!inherits(gtt,"list")){gs<-gtt[,-c(1:4)]}
    if(grepl("/",gs[1,1],fixed = T)){
      if(verbose){
        message("genotype table provided\ncalculating allele frequency from genotype")
        tmp<-apply_pb(gs,1,function(x){
          x<-as.character(x)
          tl<-strsplit(x,"/") # use also pipes or anything else
          tt<-unlist(lapply(tl,function(y){length(which(y==1))/2}))
          return(tt)
        })
      } else {
        tmp<-apply(gs,1,function(x){
          x<-as.character(x)
          tl<-strsplit(x,"/")
          tt<-unlist(lapply(tl,function(y){length(which(y==1))/2}))
          return(tt)
        })
      }
      tmp<-t(tmp)
      tmp[gs=="./."]<-NaN
      colnames(tmp)<-colnames(gs)
      rownames(tmp)<-paste0(gtt[,1],".",gtt[,2])
      tmp<-data.frame(gtt[,1:4],tmp)
      return(tmp)
    } else if (grepl(",",gs[1,1],fixed = T)){

      if(verbose){
        message("allele depth table provided\ncalculating allele frequency from allele depth")
        tmp<-apply_pb(gs,1,function(snp1){
          y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
          tt<-proportions(as.matrix(y),margin = 1)[,2]
          return(tt)
        })
      } else {
        tmp<-apply(gs,1,function(snp1){
          y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
          tt<-proportions(as.matrix(y),margin = 1)[,2]
          return(tt)
        })
      }
      tmp<-t(tmp)
      tmp[gs=="./."]<-NaN
      colnames(tmp)<-colnames(gs)
      rownames(tmp)<-paste0(gtt[,1],".",gtt[,2])
      tmp<-data.frame(gtt[,1:4],tmp)
      return(tmp)
    }

  } else if(ty=="pop"){
    if(inherits(gtt,"list")){
      if(length(gtt)>1){ad<-gtt$AD[,-c(1:4)]
      gt<-gtt$GT[,-c(1:4)]
      #### for GT ###

      if(verbose){
        message("AF from genotype")
        tmp.g<-apply_pb(gt,1,function(x){
          x<-as.character(x)
          tl<-strsplit(x,"/")
          tt<-mean(unlist(lapply(tl,function(y){length(which(y==1))/2})),na.rm = T)
          return(tt)
        })
      } else {
        tmp.g<-apply(gt,1,function(x){
          x<-as.character(x)
          tl<-strsplit(x,"/")
          tt<-mean(unlist(lapply(tl,function(y){length(which(y==1))/2})),na.rm = T)
          return(tt)
        })
      }
      #########
      ### for AD ###

      if(verbose){
        message("AF from allele depth")
        tmp.a<-apply_pb(ad,1,function(snp1){
          y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
          tt<-colMeans(proportions(as.matrix(y),margin = 1),na.rm = T)[2]
          return(tt)
        })
      } else {
        tmp.a<-apply(ad,1,function(snp1){
          y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
          tt<-colMeans(proportions(as.matrix(y),margin = 1),na.rm = T)[2]
          return(tt)
        })
      }
      return(list(AD=tmp.a,GT=tmp.g))

      } else if(length(gtt)==1){ gs<-gtt[1][,-c(1:4)]
      if(grepl("/",gs[1,1],fixed = T)){

        if(verbose){
          message("genotype table provided\ncalculating allele frequency from genotype")
          tmp<-apply_pb(gs,1,function(x){
            x<-as.character(x)
            tl<-strsplit(x,"/")
            tt<-mean(unlist(lapply(tl,function(y){length(which(y==1))/2})),na.rm = T)
            return(tt)
          })
        } else {
          tmp<-apply(gs,1,function(x){
            x<-as.character(x)
            tl<-strsplit(x,"/")
            tt<-mean(unlist(lapply(tl,function(y){length(which(y==1))/2})),na.rm = T)
            return(tt)
          })
        }
        return(tmp)
      } else if(grepl(",",gs[1,1],fixed = T)){

        if(verbose){
          message("allele depth table provided\ncalculating allele frequency from allele depth")
          tmp<-apply_pb(gs,1,function(snp1){
            y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            tt<-colMeans(proportions(as.matrix(y),margin = 1),na.rm = T)[2]
            return(tt)
          })
        } else {
          tmp<-apply(gs,1,function(snp1){
            y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            tt<-colMeans(proportions(as.matrix(y),margin = 1),na.rm = T)[2]
            return(tt)
          })
        }
        return(tmp)
      }
      }
    } else {
      if(any(colnames(gtt)=="CHROM")){gs<-gtt[,-c(1:4)]} else {gs<-gtt}
      if(grepl("/",gs[1,1],fixed = T)){

        if(verbose){
          message("genotype table provided\ncalculating allele frequency from genotype")
          tmp<-apply_pb(gs,1,function(x){
            x<-as.character(x)
            tl<-strsplit(x,"/")
            tt<-mean(unlist(lapply(tl,function(y){length(which(y==1))/2})),na.rm = T)
            return(tt)
          })
        } else {
          tmp<-apply(gs,1,function(x){
            x<-as.character(x)
            tl<-strsplit(x,"/")
            tt<-mean(unlist(lapply(tl,function(y){length(which(y==1))/2})),na.rm = T)
            return(tt)
          })
        }
        return(tmp)
      } else if(grepl(",",gs[1,1],fixed = T)){

        if(verbose){
          message("allele depth table provided\ncalculating allele frequency from allele depth")
          tmp<-apply_pb(gs,1,function(snp1){
            y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            tt<-colMeans(proportions(as.matrix(y),margin = 1),na.rm = T)[2]
            return(tt)
          })
        } else {
          tmp<-apply(gs,1,function(snp1){
            y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
            tt<-colMeans(proportions(as.matrix(y),margin = 1),na.rm = T)[2]
            return(tt)
          })
        }
        return(as.data.frame(tmp))
      }
    }
  }
}
