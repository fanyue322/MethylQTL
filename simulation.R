library(Matrix)
library(MASS)

cal_K<-function(n,m,geno)
{
  genotype=t(geno$hap1)+t(geno$hap2)
  N=nrow(genotype)
  geno=matrix(0,ncol=n,nrow=N)
  for(i in 1:N)
  {
    idx=which(genotype[i,]!=-9)
    geno_mean=mean(genotype[i,idx])
    tmp=genotype[i,]
    if(length(idx)!=n)
    {
      tmp[-idx]=geno_mean
    }
    tmp=tmp-geno_mean
    tmp=tmp/sd(tmp)*sqrt(n/(n-1))
    geno[i,]=tmp
  }
  K=matrix(0,nrow=n,ncol=n)
  N=nrow(geno)
  cout=0
  for(i in 1:N)
  {
    if(length(which(is.na(geno[i,])))==0)
    {
      K=K+ geno[i,]%*%t(geno[i,])
      cout=cout+1
    }
  }
  
  K=K/N
  return(K)
}


SimuGeno<-function(n,m,MAF)
{
  load(paste0('genotype_',n,'.RData'))   #############genotype data from WTCCC
  idx=which(maf>0.5)
  for(i in idx)
  {
    maf[i]=1-maf[i]
  }
  idx=intersect(which(maf<(MAF+0.05)),which(maf>(MAF-0.05)))
  idxx=sample(idx,m,replace=TRUE)
  genotype=genotype[idxx,]
  geno <- list()
  geno[[1]]<-matrix(0,nrow=n,ncol=m)
  geno[[2]]<-matrix(0,nrow=n,ncol=m)

  for(i in 1:m)
  {
    tmp=genotype[i,]
    idx1=which(tmp==1)
    idx2=which(tmp==2)
    tmp1=rep(0,n)
    tmp2=rep(0,n)
    tmp1[idx2]=1
    tmp2[idx2]=1
    rtmp=1
    for(j in idx1)
    {
      rtmp=runif(1,0,1)
      if(rtmp>0.5)
      {
        tmp1[j]=1
      }else{
        tmp2[j]=1
      }
    }
    geno[[1]][,i]<-tmp1
    geno[[2]][,i]<-tmp2
  }
  names(geno) <- c('hap1', 'hap2')
  return(geno)
}


logit <- function(x) {
  log(x/(1-x))
}
logistic <- function(x) {
  exp(x)/(exp(x)+1)
}

MergeLowreads <- function(Merge, thres=20) {
  while(Merge$r[1] < thres) {
    targ = which(Merge$x == Merge$x[1])[2]
    Merge[targ, 1:2] <- Merge[targ, 1:2] + Merge[1, 1:2]
    Merge = Merge[with(Merge, order(r))[-1], ]
  }
  return(Merge)
}

SimuRead <- function(geno, m1, m2, NBsize, MeanTR, K=NULL,
                     h=0, rho=0, pve=0.1, sigma=1, mu=0.5,
                     FixTR=FALSE, EvenSP=FALSE) {
  m <- m1 + m2
  n <- nrow(geno[[1]])
  g1 <- scale(geno[[1]], scale=FALSE)
  g2 <- scale(geno[[2]], scale=FALSE)
  
  total <- sapply(1:m, function(i) {
    rnbinom(n=n, size=NBsize, mu=MeanTR)
  })
  if (FixTR != FALSE) {
    total <- matrix(FixTR, nrow=n, ncol=m)
  }
  data = list()
  data[[1]] <- total
  data[[2]] <- 0
  
  # r1
  data[[3]] <- apply(total, 2, function(x) {
    qj <- rbeta(n=n, shape1=10, shape2=10)
    r1 <- sapply(x, function(r) rbinom(n=1, size=r, prob=qj))
    return(r1)
  })
  #r2
  data[[4]] <- total - data[[3]]
  #  tmp=which(data[[4]]==0)
  
  
  
  if (EvenSP == TRUE) {
    data[[3]] <- apply(total, 2, function(x) {
      r1 <- sapply(x, function(r) round(r*0.5))
      return(r1)
    })
    data[[4]] <- total - data[[3]]
  }
  
  ### Simulate methylated read counts
  asm = c(rep(1, m1), rep(0, m2))
  
  Sigma = matrix(rho*(1-pve-h)*sigma, nrow=2, ncol=2)
  diag(Sigma) <- rep((1-pve-h)*sigma, 2)
  err = mvrnorm(n=n*m, mu=c(0, 0), Sigma=Sigma)
  error = list(hap1 = matrix(err[, 1], nrow=n, ncol=m),
               hap2 = matrix(err[, 2], nrow=n, ncol=m))
  
  if (!is.null(K)) {
    # Genetic correlation
    K <- h*sigma*(K)
    g <- mvrnorm(n=m, mu=rep(0, n), Sigma=K)
    
    error[[1]] <- error[[1]] + t(g)
    error[[2]] <- error[[2]] + t(g)
  }
  
  beta <- sapply(1:m, function(i) {
    b <- sqrt(pve*sigma/(max(var(g1[, i]),var(g2[, i]))))
    return(c(mu, b*asm[i]))
  })
  
  # y1
  data[[5]] <- sapply(1:m, function(i) {
    pi1 <- logistic( logit(mu) + g1[, i] * beta[2, i] + error[[1]][, i])
    return(sapply(1:n, function(j) {
      return(rbinom(1, data[[3]][j, i], pi1[j]))
    }) )
  })
  # y2
  data[[6]] <- sapply(1:m, function(i) {
    pi2 <- logistic( logit(mu) + g2[, i] * beta[2, i] + error[[2]][, i])
    return(sapply(1:n, function(j) {
      return(rbinom(1, data[[4]][j, i], pi2[j]))
    }) )
  })
  
  data[[2]] <- data[[5]] + data[[6]]
  
  names(data) <- c('r', 'y', 'r1', 'r2', 'y1', 'y2')
  return(data)
}

  #####################
n=50
h=3/17 ############ 0.3 for individual level
rho=0
pve=0.1
sigma=0.7
mu=0.5
TR=20
###################
  m1=1000
  m2=9000
  m=m1+m2
  geno <- SimuGeno(n, m=m, MAF=0.3) #########change MAF here
  K<-cal_K(n,m,geno)
  data <- SimuRead(geno, m1=m1, m2=m2, NBsize=3, MeanTR=TR,
                   K = K, h=h,
                   rho=rho, pve=pve, sigma=sigma, mu=mu,
                   FixTR=FALSE)
} 
