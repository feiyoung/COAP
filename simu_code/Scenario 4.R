

rm(list=ls())
library(pracma)
source("simu_util_zqz.R")
library(COAP)
library(doParallel)
library(foreach)
n <- 100  ; p <- 200; d <- 50; q <- 5; rank0 <- 6  
N <- 200
tr_sta<-function(A_hat,A)
{
  return(sum(diag(t(A)%*%A_hat%*%ginv(t(A_hat)%*%A_hat)%*%t(A_hat)%*%A))/sum(diag(t(A)%*%A)))
}
methodNames <- c("APFactor", "LFM","GLLVM")
n_methods <- length(methodNames)
metricList <- list(H_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                   H_CCor=matrix(NA,N, n_methods), B_CCor=matrix(NA, N, n_methods),
                   Upsilon_tr=matrix(NA, N, n_methods),Upsilon_CCor=matrix(NA, N, n_methods),  
                   mu_norm=matrix(NA, N, n_methods), beta_norm = matrix(NA, N, n_methods),
                   timeMat = matrix(NA, N, n_methods))
metricList_no <- list(H_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                      H_CCor=matrix(NA,N, n_methods), B_CCor=matrix(NA, N, n_methods),
                      Upsilon_tr=matrix(NA, N, n_methods),Upsilon_CCor=matrix(NA, N, n_methods), 
                      mu_norm=matrix(NA, N, n_methods), timeMat = matrix(NA, N, n_methods))
library(parallel)
cl<- makeCluster(5)
registerDoParallel(cl)
result<-foreach(i=1:N,.combine="rbind",.packages = c("COAP","gllvm")) %dopar%
{
  #i <- 1
  hq <- q
  rankr<-rank0-3
  message('i = ', i)
  datList <- gendata_simu_lowrank(n=n, p=p, d=d,rank0 = rank0,q=q, rho=c(1, 2,10)/2, seed = i, sigma2_eps=1)
  X_count <- datList$X; Z <-datList$Z; bbeta0<-datList$bbeta0; H0 <- datList$H0; B0 <- datList$B0
  
  tic <- proc.time()
  reslist_no <- RR_COAP(datList$X_no, Z=matrix(1, n, 1), q= hq, rank_use =1, epsELBO = 1e-7)
  toc <- proc.time()
  time_apfactor_no <- toc[3] - tic[3]
  re1<-rep(NA,8)
  re1[1] <- time_apfactor_no
  re1[2]<- tr_sta(reslist_no$H, H0)
  re1[3] <- tr_sta(reslist_no$B, B0)
  re1[4]<- tr_sta(cbind(reslist_no$bbeta, reslist_no$B ), cbind(bbeta0[,1],B0))
  
  re1[5]<- MCCor(reslist_no$H, H0)
  re1[6] <- MCCor(reslist_no$B , B0)
  re1[7] <- MCCor(cbind(reslist_no$bbeta, reslist_no$B ), cbind(bbeta0[,1],B0))
  re1[8]<- norm_vec(reslist_no$bbeta- bbeta0[,1])
  
  
  tic <- proc.time()
  reslist <-  RR_COAP(X_count, Z=cbind(matrix(1, n, 1),Z), q= hq, rank_use = rankr, epsELBO = 1e-7)
  toc <- proc.time()
  time_apfactor <- toc[3] - tic[3]
  re2<-rep(NA,9)
  re2[1] <- time_apfactor
  re2[2]<- tr_sta(reslist$H, H0)
  re2[3] <- tr_sta(reslist$B, B0)
  re2[4]<- tr_sta(cbind(reslist$bbeta[,1], reslist$B ), cbind(bbeta0[,1],B0))
  
  re2[5] <- MCCor(reslist$H, H0)
  re2[6]<- MCCor(reslist$B , B0)
  re2[7] <- MCCor(cbind(reslist$bbeta[,1], reslist$B ), cbind(bbeta0[,1],B0))
  re2[8] <- norm_vec(reslist$bbeta[,1]- bbeta0[,1])
  re2[9] <- norm_vec(as.vector(reslist$bbeta) - as.vector(bbeta0)) 
  c(re1,re2)
  
}
















