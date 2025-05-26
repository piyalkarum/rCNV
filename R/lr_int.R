##### internal functions of likelihood ratio calculations ######
nb_stats <- function(tot.tab,cl=cl){
  ##### correlation between mean depth and sd, calculate parameter for nb distribution
  colstat <- parApply(tot.tab, 2, FUN = function(x){
    mean <- mean(x[x> quantile(x,0.05) & x < quantile(x,0.95)])
    sd <- sd(x[x> quantile(x,0.05) & x < quantile(x,0.95)])
    size <- fitdistr(x, "Negative Binomial")$estimate["size"]
    return(data.frame(mean,sd,size))
  },cl=cl)

  par(mfrow =c(4,5))
  colstat <- do.call(rbind,colstat)
  plot(colstat$mean,colstat$sd,xlab = "mean",ylab = "sd") ## check correlation


  par(mfrow =c(4,5))
  fit <- lm(sd~mean,data = colstat)
  colstat$mean2 <- 2*colstat$mean              ## expected mean and depth for N=4
  colstat$sd2 <- predict(fit,newdata = data.frame(mean = colstat$mean2))
  colstat$size2 <- colstat$mean2^2/(colstat$sd2^2-colstat$mean2)       ## expected mean and depth for N=4

  ############# check if nb distribution fit well, output first 20, can be disabled
  for (i in 1:19) {
    hist(tot.tab[,i], breaks = max(tot.tab[,i]), freq = FALSE, col = "gray", main = "Negative Binomial Fit",
         xlab = "Counts", xlim = c(0, 200))
    suppressWarnings(lines(0:200, dnbinom(0:200, size = colstat$size[i], mu = colstat$mean[i]),
                           col = "red", lwd = 2))
    suppressWarnings(lines(0:200, dnbinom(0:200, size = colstat$size2[i], mu = colstat$mean2[i]),
                           col = "black", lwd = 2))
  }
  return(colstat)
}

###### calculate likelihood ratio with nb distribution modelling
cal_depth_lld_indi2 <- function(tot.tab,nb_stat){
  ### N = 2
  dep_class <- unique(tot.tab)
  prob <- dnbinom(dep_class, size = nb_stat$size, mu = nb_stat$mean)
  names(prob) <- as.character(dep_class)
  prob_each2 <- prob[as.character(tot.tab)]
  ### N = 4
  prob <- dnbinom(dep_class, size = nb_stat$size2, mu = nb_stat$mean2)
  names(prob) <- as.character(dep_class)
  prob_each4 <- prob[as.character(tot.tab)]
  return(log(prob_each4) - log(prob_each2))
}


###### calculate likelihood ratio based on allelic ratio
cal_geno_lld2 <- function(ad.tab,inb,nb_stat){
  y <- do.call(rbind,lapply(ad.tab,FUN=function(geno){
    print(geno)
    tem <- as.numeric(unlist(strsplit(geno,split = ",")))
    if(any(tem==1) & abs(tem[1]-tem[2]) >10) {
      tem[tem==1] <- 0
    }
    if (tem[1] == 0 & tem[2]!=0 ) tem[3] = 2  ##### consider calculate from genotype
    if (tem[1] != 0 & tem[2]==0 ) tem[3] = 0
    if (tem[1] != 0 & tem[2]!=0 ) tem[3] = 1
    if (tem[1] == 0 & tem[2]==0 ) tem[3] = -9
    return(tem)
  }))

  #### n=2
  p.alt <- sum(y[,3][y[,3]!= -9])/2/sum(y[,3]!= -9)
  p.ref <- 1- p.alt
  pG.alt.hom <- p.alt^2*(1-inb) + p.alt*inb   #### consider apply population info
  pG.ref.hom <- p.ref^2*(1-inb) + p.ref*inb
  pG.het <- 2*p.alt*p.ref*(1-inb)
  lld2 <- apply(y[,1:2], 1, FUN = function(dep) {
    if(sum(dep) == 0 ){ lld = 1} else{
      ##### genotype likelihood p(D|G), choose(sum(dep),dep[1]) is going to be cancel out
      sampling <- choose(sum(dep),dep[1])
      pDG.ref.hom <- 0.98^dep[1]*0.02^dep[2]
      pDG.alt.hom <- 0.02^dep[1]*0.98^dep[2]
      pDG.het <- 0.5^dep[1]*0.5^dep[2]

      ##### P(D|G) * P(G)
      p.alt.hom <- pG.alt.hom * pDG.alt.hom
      p.ref.hom <- pG.ref.hom * pDG.ref.hom
      p.het <- pG.het * pDG.het
      p.all <- p.alt.hom + p.ref.hom + p.het

      ##### expected genotype frequency from posterior probability
      p.alt.hom <- p.alt.hom/(p.all)
      p.ref.hom <- p.ref.hom/(p.all)
      p.het <- p.het/(p.all)

      lld <- (pDG.ref.hom * p.ref.hom +
                pDG.het * p.het +
                pDG.alt.hom *p.alt.hom) * sampling
    }
  })

  ##### N=4
  tem <- y[,1:2]
  tem[,1] <- tem[,1]/nb_stat$mean
  tem[,2] <- tem[,2]/nb_stat$mean
  p.alt <- sum(tem[,2])/sum(tem)
  p.ref <- 1- p.alt
  ##### genotype likelihood P(G)
  pG.aaaa <- p.alt^4*(1-inb) + p.alt*inb
  pG.Aaaa <- 4*p.ref*p.alt^3*(1-inb)
  pG.AAaa <- 6*p.ref^2*p.alt^2*(1-inb)
  pG.AAAa <- 4*p.ref^3*p.alt*(1-inb)
  pG.AAAA <- p.ref^4*(1-inb) + p.ref*inb
  lld4 <- apply(y[,1:2], 1, FUN = function(dep) {
    if(sum(dep) == 0 ){ lld = 1} else{
      ##### genotype likelihood P(D|G)
      sampling <- choose(sum(dep),dep[1])
      pDG.aaaa <- 0.04^dep[1]*0.96^dep[2]
      pDG.Aaaa <- 0.27^dep[1]*0.73^dep[2]
      pDG.AAaa <- 0.5^dep[1]*0.5^dep[2]
      pDG.AAAa <- 0.73^dep[1]*0.27^dep[2]
      pDG.AAAA <- 0.96^dep[1]*0.04^dep[2]

      #### P(D|G) * P(G)
      p.aaaa <- pG.aaaa * pDG.aaaa
      p.Aaaa <- pG.Aaaa * pDG.Aaaa
      p.AAaa <- pG.AAaa * pDG.AAaa
      p.AAAa <- pG.AAAa * pDG.AAAa
      p.AAAA <- pG.AAAA * pDG.AAAA
      ##### expected frequencies
      p.all <- p.aaaa + p.Aaaa + p.AAaa + p.AAAa + p.AAAA
      p.aaaa <- p.aaaa/p.all
      p.Aaaa <- p.Aaaa/p.all
      p.AAaa <- p.AAaa/p.all
      p.AAAa <- p.AAAa/p.all
      p.AAAA <- p.AAAA/p.all


      lld <- (pDG.AAAA * p.AAAA +
                pDG.AAAa * p.AAAa +
                pDG.AAaa * p.AAaa +
                pDG.Aaaa * p.Aaaa +
                pDG.aaaa * p.aaaa) * sampling
    }
  })
  return(log(lld4) - log(lld2))
}


sample_and_quantile <- function(x, samples,nrep, quan = 0.95) {
  tem <- replicate(nrep,sum(sample(samples,x)))
  return(quantile(tem,quan))
}

##### not used for now
# sim.thre <- function(nb_stat, fis, f = 0.1, rep = 100000, p = 0.95, cl = NULL) {
#   # Create cluster if none provided
#   if(is.null(cl)) {
#     cl <- parallel::makeCluster(parallel::detectCores() - 1)
#     on.exit(parallel::stopCluster(cl))  # Ensure cluster is stopped
#   }
#   parallel::clusterExport(cl, varlist = c("cal_geno_lld2","geno_parser"))
#   sim.lhr <- parallel::parLapply(cl = cl, X = 1:rep, fun = function(i, nb_stat, f) {
#     N <- nrow(nb_stat)
#     mean_dep <- ceiling(nb_stat$mean)
#     N4_order <- sample(1:N, ceiling(N * f), replace = FALSE)
#     N4_0.75 <- sample(N4_order,ceiling(N * f) * 0.5)
#     N2_order <- setdiff(1:N, N4_order)
#
#     dep <- c()
#     if (all(nb_stat$size > 0 & nb_stat$size2 > 0)) {
#       dep[N4_order] <- rnbinom(length(N4_order), nb_stat$size2[N4_order], nb_stat$p2[N4_order])
#       dep[N2_order] <- rnbinom(length(N2_order), nb_stat$size[N2_order], nb_stat$p[N2_order])
#       dep[dep<0] <- 0
#       dep_lhr <- log(dnbinom(dep, nb_stat$size2, nb_stat$p2)) -
#         log(dnbinom(dep, nb_stat$size, nb_stat$p))
#     } else {
#       dep[N4_order] <- rnorm(length(N4_order), nb_stat$mean2[N4_order], nb_stat$sd2[N4_order])
#       dep[N2_order] <- rnorm(length(N2_order), nb_stat$mean[N2_order], nb_stat$sd[N2_order])
#       dep[dep<0] <- 0
#       dep4_lhr <- pnorm(dep + 0.5, nb_stat$mean2, nb_stat$sd2) - pnorm(dep - 0.5, nb_stat$mean2, nb_stat$sd2)
#       dep2_lhr <- pnorm(dep + 0.5, nb_stat$mean, nb_stat$sd) - pnorm(dep - 0.5, nb_stat$mean, nb_stat$sd)
#       dep4_lhr[dep4_lhr == 0] <- 1e-20
#       dep2_lhr[dep2_lhr == 0] <- 1e-20
#       dep_lhr <- log(dep4_lhr) - log(dep2_lhr)
#     }
#
#
#     ref <- rbinom(N, ceiling(dep), 0.5)
#     ref[N4_0.75] <- rbinom(length(N4_0.75), ceiling(dep[N4_0.75]), 0.75)
#     alt <- ceiling(dep) - ref
#     ad <- paste(ref,alt,sep = ",")
#     rat_lhr <- cal_geno_lld2(ad.tab = ad,inb = 0,nb_stat = data.frame(mean = rep(1,length(ad))))
#     lhr.sim <- 2*(dep_lhr + rat_lhr)
#     lhr.sim[lhr.sim < 0] <- 0
#     return(sum(lhr.sim))
#   }, nb_stat = nb_stat, f = f)
#
#   return(mean(unlist(sim.lhr), p,na.rm = T))
# }
#

