#helpers
# 1. moving window average
wind<-function(xx,dd){
  d<-dd[dd[,1]==xx,c("POS","dup.stat")]
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
#' @param d.detect a data frame of detected duplicates or deviants from the outputs of \code{dupGet} or \code{cnv}
#'  (output of \code{dupGet})
#' @param window.size numerical. a single value of the desired moving window
#'  size (default \code{100} bp)
#' @param scaf.size numerical. scaffold size to be checked. i.e. the split size of the fragment to be checked with the specified window size.
#' default=10000
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
dup.validate<-function(d.detect,window.size=100, scaf.size=10000){
  nm<-unique(d.detect[,1])
  gg<-lapply(nm,wind,dd=d.detect)
  names(gg)<-nm
  means<-lapply_pb(nm,function(x,mw,gg){yy<-unlist(gg[names(gg)==x])
  ll<-mw+(mw/2)
  if(length(yy)>ll){
    if(length(yy)>2*scaf.size){
      y.list<-split(yy, ceiling(seq_along(yy)/scaf.size))
      yl<-lapply(y.list,function(v,mw){
        vv<-unlist(v)
        if(length(vv)>ll){
          dpp<-NULL
          for(i in 1:(length(vv)-mw)){
            tmp<-unlist(vv)[i:(i+mw)]
            ss<-sum(tmp=="non-cnv" | tmp=="non-deviant")
            dd<-sum(tmp=="cnv" | tmp=="deviant")
            dpp[i]<-dd/(dd+ss)
          }
          dpp[dpp==0]<-NA
          return(cbind(mean(dpp,na.rm = TRUE),length(yy),sum(vv=="cnv" | vv=="deviant"),sum(vv=="non-cnv" | vv=="non-deviant")))
        }
      },mw=mw)
      dp<-do.call(rbind,yl)
      dp<-cbind(paste0(x,".",1:length(y.list)),dp)
    } else {
      dp<-NULL
      for(i in 1:(length(yy)-mw)){
        tmp<-unname(unlist(yy)[i:(i+mw)])
        ss<-sum(tmp=="non-cnv" | tmp=="non-deviant")
        dd<-sum(tmp=="cnv" | tmp=="deviant")
        dp[i]<-dd/(dd+ss)
      }
      dp[dp==0]<-NA
      dp<-mean(dp,na.rm = TRUE)
      dp<-cbind(x,dp,length(yy),sum(yy=="cnv" | yy=="deviant"),sum(yy=="non-cnv" | yy=="non-deviant"))
    }
  } else {
    dp<-sum(yy=="cnv" | yy=="deviant")/(sum(yy=="cnv" | yy=="deviant")+sum(yy=="non-cnv" | yy=="non-deviant"))
    dp<-cbind(x,dp,length(yy),sum(yy=="cnv" | yy=="deviant"),sum(yy=="non-cnv" | yy=="non-deviant"))
  }
  return(dp)
  },mw=window.size,gg=gg)
  dup.ratio<-do.call(rbind,means)
  dup.ratio[dup.ratio=="NaN"]<-0
  colnames(dup.ratio)<-c("CHROM","dp.ratio","CHROM.length","no.duplicates","no.singlets")
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
#' @param verbose logical. show progress
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
#' ds<-DD[DD$dup.stat=="deviant",]
#' ad<-ADtable[match(paste0(ds$CHROM,".",ds$POS),paste0(ADtable$CHROM,".",ADtable$POS)),]
#' vst(ad,pops=substr(colnames(ad)[-c(1:4)],1,11))}
#'
#' @export
vst<-function(AD,pops,id.list=NULL,qGraph=TRUE,verbose=TRUE,...){
  if(!is.null(id.list)){
    AD<-AD[match(id.list,AD$ID),]
  }
  nm<-colnames(data.frame(AD))[-c(1:4)]
  pop<-na.omit(unique(pops))
  AD<-AD[,-c(1:4)]
  if(is.character(AD[1,1])){
    tm<-apply(AD,2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))))}))})
    AD<-tm
  }

  AD[AD==0]<-NA
  tmp<-data.frame(ind=nm,pop=pops,t(AD))
  # Vst - for CNVs
  # Vt-Vs/Vt
  # VT is the variance of normalized read depths among all individuals from the two populations and VS is the average of the variance within each population, weighed for population size
 if(verbose){
   Vst<-combn_pb(pop,2,function(x){
     jj<-tmp[tmp$pop==x[1],-c(1:2)]
     kk<-tmp[tmp$pop==x[2],-c(1:2)]
     ft<-ncol(jj)
     ll<-lapply(1:ft, function(y){
       vt<-var(c(jj[,y],kk[,y]),na.rm=TRUE)
       vs<-(var(jj[,y],na.rm=TRUE)*length(na.omit(jj[,y]))+
              var(kk[,y],na.rm=TRUE)*length(na.omit(kk[,y])))/(length(na.omit(jj[,y]))+length(na.omit(kk[,y])))
       return((vt-vs)/vt)
     })
     vst<-mean(unlist(ll),na.rm = TRUE)
     return(matrix(vst,dimnames=list(x[1],x[2])))
   },simplify = FALSE)
 } else {
   Vst<-combn(pop,2,function(x){
     jj<-tmp[tmp$pop==x[1],-c(1:2)]
     kk<-tmp[tmp$pop==x[2],-c(1:2)]
     ft<-ncol(jj)
     ll<-lapply(1:ft, function(y){
       vt<-var(c(jj[,y],kk[,y]),na.rm=TRUE)
       vs<-(var(jj[,y],na.rm=TRUE)*length(na.omit(jj[,y]))+
              var(kk[,y],na.rm=TRUE)*length(na.omit(kk[,y])))/(length(na.omit(jj[,y]))+length(na.omit(kk[,y])))
       return((vt-vs)/vt)
     })
     vst<-mean(unlist(ll),na.rm = TRUE)
     return(matrix(vst,dimnames=list(x[1],x[2])))
   },simplify = FALSE)
 }
  mt<-matrix(NA,nrow=length(pop),ncol=length(pop))
  dimnames(mt)<-list(pop,pop)
  for(i in seq_along(Vst)){
    mt[colnames(Vst[[i]]),rownames(Vst[[i]])]<-Vst[[i]]
  }
  if(qGraph){
    suppressWarnings(qgraph::qgraph(1/mt,layout="spring", ...=...))
  }
  return(mt)
}


#' Run permutation on Vst
#'
#' This function runs a permutation test on Vst calculation
#'
#' @param AD data frame of total allele depth values of SNPs
#' @param pops character. A vector of population names for each individual.
#' Must be the same length as the number of samples in AD
#' @param nperm numeric. Number of permutations to perform
#' @param qGraph logical. Plot the network plot based on observed Vst values
#' (see \code{vst()} help page for more details)
#' @param histogram logical. plots the distribution histogram of permuted vst values vs. observed values
#'
#' @importFrom qgraph qgraph
#' @importFrom graphics hist
#'
#' @return Returns a list with observed vst values, an array of permuted vst values and the p-values for the permutation test
#'
#'
#' @author Piyal Karunarathne, Jorge Cortés-Miranda (email:)
#'
#' @examples
#' \dontrun{data(alleleINF)
#' data(ADtable)
#' DD<-dupGet(alleleINF)
#' ds<-DD[DD$dup.stat=="deviant",]
#' ad<-ADtable[match(paste0(ds$CHROM,".",ds$POS),paste0(ADtable$CHROM,".",ADtable$POS)),]
#' vstPermutation(ad,pops=substr(colnames(ad)[-c(1:4)],1,11))}
#'
#' @export
vstPermutation<-function(AD,pops,nperm=100,histogram=TRUE,qGraph=TRUE){
  npop<-length(unique(pops))
  message("Permuting Vst")
  perm_vst<-sapply_pb(1:nperm,function(x){
    ind_cnv <- AD[-1:-4]
    ind_test <- sample(ind_cnv, size = length(ind_cnv), replace = F)
    dd_test <- cbind(AD[1:4], ind_test)
    Vs_sim <-  vst(dd_test, pops= pops, qGraph = FALSE,verbose = FALSE)
  },simplify="array")

  message("Calculating ovserved Vst")
  obs_vst<-vst(AD, pops, qGraph = qGraph,verbose = TRUE)

  out_mat<-matrix(0, nrow = npop, ncol = npop)
  for(i in 1:nperm){
    matrix1=perm_vst[,,i]
    comparison_matrix <- matrix(0, nrow = nrow(matrix1), ncol = ncol(matrix1))
    comparison_matrix[!is.na(matrix1) & !is.na(obs_vst)] <- as.numeric(matrix1[!is.na(matrix1)] > obs_vst[!is.na(obs_vst)])
    out_mat<-out_mat+comparison_matrix
  }
  p_mat<-out_mat/nperm
  p_mat[upper.tri(p_mat)] <- NA

  if(histogram){
    h_dat<-perm_vst[2,1,1:nperm]
    hist(h_dat, main = "Vst Permutation Test Histogram",
         xlab = "VST", ylab = "Frequency",
         col = "dodgerblue", border = F)
    abline(v = obs_vst[2,1], col = 2, lwd = 2)
    legend("topright", legend = c("Observed Vst", "Permutation Vst"),
           col = c(2, "dodgerblue"), lwd = 2)
  }
  return(list(observed_vst=obs_vst,permuted_vst=perm_vst,pvalue=p_mat))
}

