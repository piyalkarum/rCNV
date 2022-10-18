## deprecated (older versions) functions
get.pvals0<-function(x,df,p.cal){
  snp1<-df[x,-c(1:4)]
  y<-data.frame(do.call(rbind,strsplit(as.character(unlist(snp1)),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
  rr1<-y[,2]/rowSums(y,na.rm = TRUE)
  snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==TRUE),]
  homalt<-sum(rr1==1,na.rm=TRUE)
  homref<-sum(rr1==0,na.rm=TRUE)
  Nsamp <- nrow(snp1het)+homalt+homref
  if(nrow(snp1het)>3){
    propHet<-nrow(na.omit(snp1het))/length(na.omit(rr1))
    medRatio<-median(proportions(as.matrix(snp1het),margin = 1)[,2],na.rm = TRUE)
    p.sum<-p.cal[x,2]
    p.05<-0.5
    p.all<-p.cal[x,1]
    n<-unname(rowSums(snp1het,na.rm = TRUE))
    chi.het<-sum((n*p.sum-snp1het[,2])^2/(n*p.sum),na.rm = TRUE)
    chi.het.sum<-chi.het
    chi.05<-sum((n*p.05-snp1het[,2])^2/(n*p.05),na.rm = TRUE)
    chi.0.5.sum<-chi.05
    chi.all<-sum((n*p.all-snp1het[,2])^2/(n*p.all))
    chi.all.sum<-chi.all
    z <- (n*p.sum-snp1het[,2])/sqrt(n*p.sum*(1-p.sum))
    z.het.sum<-sum(z,na.rm = TRUE)
    z<-pnorm(sum(z,na.rm = TRUE),0,sqrt(nrow(snp1het)))
    z.05 <- (n*p.05-snp1het[,2])/sqrt(n*p.05*(1-p.05))
    z.05.sum<-sum(z.05,na.rm = TRUE)
    z.05<-pnorm(sum(z.05,na.rm = TRUE),0,sqrt(nrow(snp1het)))
    z.all<- (n*p.all-snp1het[,2])/sqrt(n*p.all*(1-p.all))
    z.all.sum<-sum(z.all,na.rm = TRUE)
    z.all<-pnorm(sum(z.all,na.rm = TRUE),0,sqrt(nrow(snp1het)))
    ll<-data.frame(NHet=nrow(snp1het),propHet,medRatio,NHomRef=homref,NHomAlt=homalt,propHomAlt=homalt/Nsamp,Nsamp,
                   pAll=p.all,pHet=p.sum,fis=1-(nrow(snp1het)/(2*(homref+(nrow(snp1het)/2))*(homalt+(nrow(snp1het)/2)))),
                   z.het=ifelse(z>0.5, (1-z)*2, z*2),
                   z.05=ifelse(z.05>0.5, (1-z.05)*2, z.05*2),
                   z.all=ifelse(z.all>0.5, (1-z.all)*2, z.all*2),
                   chi.het=pchisq(chi.het,nrow(snp1het)-1,lower.tail=F),
                   chi.05=pchisq(chi.05,nrow(snp1het)-1,lower.tail = F),
                   chi.all=pchisq(chi.all,nrow(snp1het)-1,lower.tail = F),
                   z.het.sum,z.05.sum,z.all.sum,chi.het.sum,chi.0.5.sum,chi.all.sum)
  } else {
    ll<-NA
  }
  return(ll)
}




gt.format0 <- function(gt,info,format=c("benv","bpass"),snp.subset=NULL,verbose=FALSE) {
  if(is.character(gt)){
    gt <-as.data.frame(fread(gt))
    gts <-gt[,-c(1,2)]
  } else {
    gts<-gt[,-c(1:4)]
  }
  if(is.character(info)){
    pop.col<-NULL
    for(i in seq_along(info)){
      pop.col[grep(info[i],colnames(gts))]<-info[i]
    }
    info<-data.frame(population=pop.col)
  }
  rownames(gts)<-paste(gt$CHROM,gt$POS,sep=".")
  pp<-na.omit(unique(info$population))
  infos<-as.character(info$population)
  format<-match.arg(format,several.ok = TRUE)
  if(any(format=="bpass")){
    if(verbose){
      message("formating BayPass")
      pgt.b<-lapply_pb(pp,function(x,gts,info){
        tm <- as.data.frame(gts[,which(infos==x)])
        gtt <- lapply(1:nrow(tm),function(y,tm){gg(as.character(tm[y,]))},tm=tm)
        gf <- do.call(rbind,gtt)
        return(gf)
      },gts=gts,info=infos)
    } else {
      pgt.b<-lapply(pp,function(x,gts,info){
        tm <- as.data.frame(gts[,which(infos==x)])
        gtt <- lapply(1:nrow(tm),function(y,tm){gg(as.character(tm[y,]))},tm=tm)
        gf <- do.call(rbind,gtt)
        return(gf)
      },gts=gts,info=infos)
    }
    pgt.b<-do.call(cbind,pgt.b)
    nm <- unlist(lapply(pp,FUN=function(x)c(paste0(x,"~1"),paste0(x,"~2"))))
    colnames(pgt.b)<-nm
    pgt.b <- as.matrix(pgt.b)
    pgt.b[which(is.na(pgt.b))] <- 0
    rownames(pgt.b)<-paste0(gt$CHROM,".",gt$POS)
    if(!is.null(snp.subset)){
      rn<-sample(1:snp.subset,nrow(gts),replace = T)
      chu.b<-lapply(1:snp.subset,function(nn,pgt.b,rn){
        tmp<-pgt.b[which(rn==nn),]
        return(tmp)
      },pgt.b=pgt.b,rn=rn)
    } else { chu.b <- NULL}
  }

  if(any(format=="benv")){
    if(verbose){
      message("formating BayEnv")
      pgt.e<-lapply_pb(pp,function(x,gts,info){
        tm <- as.data.frame(gts[,which(infos==x)])
        gtt <- lapply(1:nrow(tm),function(y,tm){gg(as.character(tm[y,]))},tm=tm)
        gf <- t(do.call(cbind,gtt))
        colnames(gf)<-x
        return(gf)
      },gts=gts,info=infos)
    } else {
      pgt.e<-lapply(pp,function(x,gts,info){
        tm <- as.data.frame(gts[,which(infos==x)])
        gtt <- lapply(1:nrow(tm),function(y,tm){gg(as.character(tm[y,]))},tm=tm)
        gf <- t(do.call(cbind,gtt))
        colnames(gf)<-x
        return(gf)
      },gts=gts,info=infos)
    }

    pgt.e<-do.call(cbind,pgt.e)
    nm <- unlist(lapply(rownames(gts),FUN=function(x)c(paste0(x,"~1"),paste0(x,"~2"))))
    rownames(pgt.e)<-nm
    pgt.e <- as.matrix(pgt.e)
    pgt.e[which(is.na(pgt.e))] <- 0
    if(!is.null(snp.subset)){
      rn<-sample(1:snp.subset,nrow(gts),replace = T)
      snps<-rownames(gts)
      chu.e<-lapply(1:snp.subset,function(nn,pgt.e,snps,rn){
        sset<-snps[rn==nn]
        tmp0<-NULL
        for(k in seq_along(sset)){
          tmp<-pgt.e[grep(sset[k],rownames(pgt.e)),]
          tmp0<-rbind(tmp0,tmp)
        }
        return(tmp0)
      },pgt.e=pgt.e,snps=snps,rn=rn)
    } else { chu.e <- NULL}
  }
  return(list(baypass=pgt.b,bayenv=pgt.e,sub.bp=chu.b,sub.be=chu.e,pop=as.character(pp)))
}


het.sity_old<-function(ind){
  O<-length(which(ind=="0/0"))+length(which(ind=="1/1"))
  N<-length(which(ind=="0/0"))+length(which(ind=="1/1"))+length(which(ind=="0/1"))+length(which(ind=="1/0"))
  alle<-unname(unlist(lapply(ind,strsplit,split="/")))
  p<-length(which(alle=="0"))/(length(which(alle=="0"))+length(which(alle=="1")))
  q<-length(which(alle=="1"))/(length(which(alle=="0"))+length(which(alle=="1")))
  E<-(p^2+q^2)*N
  FF<-(O-E)/(N-E)
  return(c(O,E,N,FF))
}



cpm.normal_old <-function(het.table, method=c("TMM","TMMex","MedR","QN","pca"),logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10,verbose=TRUE,plot=TRUE){
  method<-match.arg(method)
  if(length(method)>1){method="TMM"}
  if(verbose){
    message("calculating normalization factor")
    tdep<-apply_pb(het.table[,-c(1:4)],2,function(tmp){
      tmp <- stringr::str_split_fixed(tmp, ",", 2L)
      tt <- as.integer(tmp[,1]) + as.integer(tmp[,2])
      return(tt)
    })
  } else {
    tdep<-apply(het.table[,-c(1:4)],2,function(tmp){
      tmp <- stringr::str_split_fixed(tmp, ",", 2L)
      tt <- as.integer(tmp[,1]) + as.integer(tmp[,2])
      return(tt)
    })
  }

  #find and warn about outliers
  ot<-boxplot.stats(colSums(tdep,na.rm = T))$out
  cl<-rep("dodgerblue",ncol(tdep))
  ot.ind<-which(colnames(tdep) %in% names(ot))
  cl[ot.ind]<-2
  if(length(ot)>0){
    if(plot){barplot(colSums(tdep,na.rm = T),col=cl,border=NA,xlab="sample",ylab="total depth")}
    message("OUTLIERS DETECTED\nConsider removing the samples:")
    cat(colnames(tdep)[ot.ind])}
  if(method=="TMM" | method=="TMMex"){
    nf<-norm.fact(tdep,method = method,logratioTrim=logratioTrim,sumTrim=sumTrim,Weighting=Weighting,Acutoff=Acutoff)
    if(verbose){
      message("\ncalculating normalized depth")
      out<-apply_pb(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/(nf[,1]*nf[,2]))*1e6,2)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    } else {
      out<-apply(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/(nf[,1]*nf[,2]))*1e6,2)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    }
  } else if(method=="MedR") {
    pseudo<- apply(tdep,1,function(xx){exp(mean(log(as.numeric(xx)[as.numeric(xx)>0])))})
    nf<-  apply(tdep,2,function(xx){median(as.numeric(xx)/pseudo,na.rm=T)})
    if(verbose){
      message("\ncalculating normalized depth")
      out<-apply_pb(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/nf),0)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    } else {
      out<-apply(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/nf),0)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})

    }
  } else if(method=="QN"){
    out<-t(do.call(cbind,quantile_normalisation(tdep,het.table,verbose=verbose)))
  }

  else if(method=="pca"){
    if(verbose){message("\ncalculating normalized depth")}
    new.mat <- t(tdep) ### check the direction to confirm if this step need to be done
    colmean <- colMeans(new.mat)
    colsd <- apply(new.mat, 2, sd)
    new.mat <- apply(new.mat, 2, scale,scale = TRUE) #### essential before SVD
    test.la.svd <- La.svd(new.mat)
    u <- test.la.svd$u
    d <- test.la.svd$d
    vt <- test.la.svd$vt
    ## optimal PCs to remove with d values
    ddl<-NULL
    for(i in seq_along(1:50)){
      if(i<50) ddl[i]<-d[i+1]-d[i]
    }
    ## modified Kaiser's Rule: Sample variation of eigen values smaller than 0.7 should be kept (i.e., the first eigen value < 0.7)
    rmpc<-min(which(abs(ddl)<0.7))
    #plot(d[1:50],pch=19,cex=0.5)
    #points(rmpc,d[rmpc],cex=1.5,col="red")
    d[1:rmpc] <- 0
    out <- u %*% diag(d) %*% vt
    out <- apply(out, 1, FUN = function(x){round(x*colsd + colmean,0)})#### back transform to depth matrix
    out[out<0]<-0
    out<-t(out)

    ht<-het.table[,-c(1:4)]
    tout<-NULL
    for(i in 1:ncol(ht)){
      if(verbose){
        pb <- txtProgressBar(min = 0, max = ncol(ht), style = 3, width = 50, char = "=")
        setTxtProgressBar(pb, i)
      }
      tmp <- stringr::str_split_fixed(ht[,i], ",", 2L)
      tt<-matrix(NA,nrow = nrow(tmp),ncol = 2)
      tt[,1]<-as.integer(tmp[,1])
      tt[,2]<-as.integer(tmp[,2])
      tt <- proportions(tt,margin = 1)
      tt[is.na(tt)]<-0
      tt<-tt*out[i,]
      tt<-paste0(round(tt[,1],0),",",round(tt[,2],0))
      tout<-cbind(tout,tt)
    }
    out<-t(tout)
  }
  out<-data.frame(het.table[,c(1:4)],t(out))
  colnames(out)<-colnames(het.table)
  return(list(AD=out,outliers=data.frame(column=(ot.ind+4),colnames(tdep)[ot.ind])))
}


cpm.normal0<-function(het.table, method=c("TMM","TMMex","MedR","QN"),logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10,verbose=TRUE,plot=TRUE){
  method<-match.arg(method)
  if(length(method)>1){method="TMM"}
  if(verbose){
    message("calculating normalization factor")
    tdep<-apply_pb(het.table[,-c(1:4)],2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))),na.rm=TRUE)}))})
  } else {
    tdep<-apply(het.table[,-c(1:4)],2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))),na.rm=TRUE)}))})
  }
  #find and warn about outliers
  ot<-boxplot.stats(colSums(tdep,na.rm = T))$out
  cl<-rep("dodgerblue",ncol(tdep))
  ot.ind<-which(colnames(tdep) %in% names(ot))
  cl[ot.ind]<-2
  if(length(ot)>0){
    if(plot){barplot(colSums(tdep,na.rm = T),col=cl,border=NA,xlab="sample",ylab="total depth")}
    message("OUTLIERS DETECTED\nConsider removing the samples:")
    cat(colnames(tdep)[ot.ind])}
  if(method=="TMM" | method=="TMMex"){
    nf<-norm.fact(tdep,method = method,logratioTrim=logratioTrim,sumTrim=sumTrim,Weighting=Weighting,Acutoff=Acutoff)
    if(verbose){
      message("\ncalculating normalized depth")
      out<-apply_pb(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/(nf[,1]*nf[,2]))*1e6,2)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    } else {
      out<-apply(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/(nf[,1]*nf[,2]))*1e6,2)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    }
  } else if(method=="MedR") {
    pseudo<- apply(tdep,1,function(xx){exp(mean(log(as.numeric(xx)[as.numeric(xx)>0])))})
    nf<-  apply(tdep,2,function(xx){median(as.numeric(xx)/pseudo,na.rm=T)})
    if(verbose){
      message("\ncalculating normalized depth")
      out<-apply_pb(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/nf),0)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    } else {
      out<-apply(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/nf),0)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})

    }
  } else if(method=="QN"){
    out<-t(do.call(cbind,quantile_normalisation(tdep,het.table,verbose=verbose)))
  }
  out<-data.frame(het.table[,c(1:4)],t(out))
  colnames(out)<-colnames(het.table)
  return(list(AD=out,outliers=data.frame(column=(ot.ind+4),colnames(tdep)[ot.ind])))
}

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



depthVsSample2<-function(cov.len=400,sam.len=1000,incr=c(1,1)){
  cov<- seq(1,cov.len,incr[1])
  nsamp<- seq(1,sam.len,incr[2])
  MR<-lapply_pb(cov,function(x,nsamp){
    tmp<-dp.cov2(cov.i=x,nsamp)
    return(tmp)
  },nsamp=nsamp)
  names(MR)<-cov
  return(MR)
}


depthVsSample0<-function(cov.len=400,sam.len=1000,incr=c(1,1),plot=TRUE,plot.cols=c("red","cyan")){
  cov<- seq(1,cov.len,incr[1])
  nsamp<- seq(1,sam.len,incr[2])
  MR<-lapply_pb(cov,function(x,nsamp){
    tmp<-dp.cov0(cov.i=x,nsamp)
    return(tmp)
  },nsamp=nsamp)
  MR<-do.call(cbind,MR)
  colnames(MR)<-cov
  rownames(MR)<-nsamp
  if(plot){
    plot.svd0(MR,cols=plot.cols)
  }
  return(MR)
  #saveRDS(MR,paste0(dirout,"/sim_SampVsDepth_",i,".rds"),compress = "gzip")
}



plot.svd0 <- function(MR,cols=c("red","cyan")){
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
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



dp.cov0<-function(cov.i,nsamp){
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
