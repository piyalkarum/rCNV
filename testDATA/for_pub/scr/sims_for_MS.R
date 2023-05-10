############################ SIMULATIONS ###################################
#___________________________________________________________________________

# this script contains all the simulations and their plots included in the 
# manuscript

generateHets_CN<-function(pAlt,N_pop,CN){
  geno <- c()
  population<-c()
  ### calculate the allelic ratio for each individuals carrying k alternative allele
  for (k in 1:(2*CN+1)) {
    geno[k] <- as.integer(N_pop*pAlt^(CN*2-k+1)*(1-pAlt)^(k-1)*choose(CN*2,k-1))
    population <- c(population,rep((k-1)/(CN*2),geno[k]))
  }
  depth <- 1000
  ar <- rep(-9,N_pop)
  for (i in which(population > 0 & population <1)){
    ### calculate the deviation (measured as the variance of number of alternative alleles)
    ar[i] <-(population[i]*2*CN-CN)^2
    ### calculate the deviation (measured as the difference between expected and observed allelic ratio)
    #ar[i] <-(population[i]-0.5)
  }
  return(ar)
}

### setup colors, steps and population size
color <- colorRampPalette(colors=c(1,2))(5)
color <- rCNV:::makeTransparent(color,alpha = .5)
x <- seq(0.0001,0.9999,0.0001)
N_pop=200000
### calculate mean deviation of SNPs with the duplicated allele carrying 2 copies
DevTab1<-data.frame(matrix(NA,ncol=3))
pAlt <- 0.5
for (t in 1:length(x)) {
  pDup <- x[t] # proportion of duplicated alleles
  ### number of individuals with different number of copies
  n1 <- round(N_pop*pDup^2,0)
  n2 <- round(N_pop*(1-pDup)^2,0)
  n3 <- N_pop-n1-n2
  ### calculate deviation for individuals carrying the same number of copies
  Dev_Homo_Dup <- generateHets_CN(pAlt,N_pop=n1,2)
  Dev_Homo_Sing <- generateHets_CN(pAlt,N_pop=n2,1)
  Dev_Het <- generateHets_CN(pAlt,N_pop=n3,1.5)
  ### keep only heterozygotes
  z1 <- ((Dev_Homo_Dup[Dev_Homo_Dup != -9]))
  z2 <- ((Dev_Homo_Sing[Dev_Homo_Sing != -9]))
  z3 <- ((Dev_Het[Dev_Het != -9]))
  temp <- c(z1,z2,z3)
  temp <- temp[!is.na(temp)]
  DevTab1[t,1] <- length(temp)/N_pop # proportion of apparent heterozygotes
  DevTab1[t,2] <- pDup
  DevTab1[t,3] <- sum(temp)/length(temp) # mean deviation
}

### calculate mean deviation of SNPs with the duplicated allele carrying 3 copies
DevTab2<-data.frame(matrix(NA,ncol=3))
pAlt <- 0.5
for (t in 1:length(x)) {
  pDup <- x[t] # proportion of duplicated alleles
  ### number of individuals with different number of copies
  n1 <- round(N_pop*pDup^2,0)
  n2 <- round(N_pop*(1-pDup)^2,0)
  n3 <- N_pop-n1-n2
  ### calculate deviation for individuals carrying the same number of copies
  Dev_Homo_Dup <- generateHets_CN(pAlt,N_pop=n1,3)
  Dev_Homo_Sing <- generateHets_CN(pAlt,N_pop=n2,1)
  Dev_Het <- generateHets_CN(pAlt,N_pop=n3,2)
  ### keep only heterozygotes
  z1 <- ((Dev_Homo_Dup[Dev_Homo_Dup != -9]))
  z2 <- ((Dev_Homo_Sing[Dev_Homo_Sing != -9]))
  z3 <- ((Dev_Het[Dev_Het != -9]))
  temp <- c(z1,z2,z3)
  temp <- temp[!is.na(temp)]
  DevTab2[t,1] <- length(temp)/N_pop # proportion of apparent heterozygotes
  DevTab2[t,2] <- pDup
  DevTab2[t,3] <- sum(temp)/length(temp) # mean deviation
}

### calculate mean deviation of SNPs with the duplicated allele carrying 4 copies
DevTab3<-data.frame(matrix(NA,ncol=3))
pAlt <- 0.5
for (t in 1:length(x)) {
  pDup <- x[t] # proportion of duplicated alleles
  ### number of individuals with different number of copies
  n1 <- round(N_pop*pDup^2,0)
  n2 <- round(N_pop*(1-pDup)^2,0)
  n3 <- N_pop-n1-n2
  ### calculate deviation for individuals carrying the same number of copies
  Dev_Homo_Dup <- generateHets_CN(pAlt,N_pop=n1,4)
  Dev_Homo_Sing <- generateHets_CN(pAlt,N_pop=n2,1)
  Dev_Het <- generateHets_CN(pAlt,N_pop=n3,2.5)
  ### keep only heterozygotes
  z1 <- ((Dev_Homo_Dup[Dev_Homo_Dup != -9]))
  z2 <- ((Dev_Homo_Sing[Dev_Homo_Sing != -9]))
  z3 <- ((Dev_Het[Dev_Het != -9]))
  temp <- c(z1,z2,z3)
  temp <- temp[!is.na(temp)]
  DevTab3[t,1] <- length(temp)/N_pop # proportion of apparent heterozygotes
  DevTab3[t,2] <- pDup
  DevTab3[t,3] <- sum(temp)/length(temp) # mean deviation
}

### calculate mean deviation of SNPs with the duplicated allele carrying 2 copies
DevTab4<-data.frame(matrix(NA,ncol=3))
pAlt <- 0.5
for (t in 1:length(x)) {
  pDup <- x[t] # proportion of duplicated alleles
  ### number of individuals with different number of copies
  n1 <- round(N_pop*pDup^2,0)
  n2 <- round(N_pop*(1-pDup)^2,0)
  n3 <- N_pop-n1-n2
  ### calculate deviation for individuals carrying the same number of copies
  Dev_Homo_Dup <- generateHets_CN(pAlt,N_pop=n1,5)
  Dev_Homo_Sing <- generateHets_CN(pAlt,N_pop=n2,1)
  Dev_Het <- generateHets_CN(pAlt,N_pop=n3,3)
  ### keep only heterozygotes
  z1 <- ((Dev_Homo_Dup[Dev_Homo_Dup != -9]))
  z2 <- ((Dev_Homo_Sing[Dev_Homo_Sing != -9]))
  z3 <- ((Dev_Het[Dev_Het != -9]))
  temp <- c(z1,z2,z3)
  temp <- temp[!is.na(temp)]
  DevTab4[t,1] <- length(temp)/N_pop # proportion of apparent heterozygotes
  DevTab4[t,2] <- pDup
  DevTab4[t,3] <- sum(temp)/length(temp) # mean deviation
}

all_CNV <- rbind(DevTab1,DevTab2,DevTab3,DevTab4)
all_CNV <- as.data.frame(all_CNV)
all_CNV$CN <- rep(2:5,each=9999)
colnames(all_CNV) <- c("pHet","pDup","Deviation","CN")



ggplot() + geom_point(aes(pDup,Deviation,stroke = 0,col=as.factor(CN)),data = all_CNV[all_CNV$CN==5,]) +
  geom_point(aes(pDup,Deviation,col=as.factor(CN)),stroke = 0,data = all_CNV[all_CNV$CN==4,]) +
  geom_point(aes(pDup,Deviation,col=as.factor(CN)),stroke = 0,data = all_CNV[all_CNV$CN==3,]) +
  geom_point(aes(pDup,Deviation,col=as.factor(CN)),stroke = 0,data = all_CNV[all_CNV$CN==2,]) +
  labs(x = "Proportion of duplicates (H)", y = "Expected variance in number of alternative alleles \nin heterozygotes (p = 0.5)",
       color = "Number of \ncopies (C)") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_color_manual(values = color)

