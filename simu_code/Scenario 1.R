

rm(list=ls())
library(pracma)
source("simu_util_zqz.R")
source("simu_util_ai.R")
library(COAP)
n <- 200  ; p <- 100; d <- 50; q <- 5; rank0 <- 6  ;ai<-20   
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
hq <- q
for(i in 1: N){
  # i <- 1
  message('i = ', i)
  datList <- gendata_simu_lowrank(n=n, p=p, d=d,rank0 = rank0,q=q, rho=c(1, 3,6)/2, seed = i, sigma2_eps=1)
  X_count <- datList$X; Z <-datList$Z; bbeta0<-datList$bbeta0; H0 <- datList$H0; B0 <- datList$B0
 

  #######################################
  ### Scenario 1.1                                                           #
  #######################################
  tic <- proc.time()
  reslist <-  RR_COAP(X_count, Z=cbind(matrix(1, n, 1),Z), q= hq, rank_use = rank0, epsELBO = 1e-7)
  toc <- proc.time()
  time_apfactor <- toc[3] - tic[3]
  metricList$timeMat[i,1] <- time_apfactor
  metricList$H_tr[i,1] <- tr_sta(reslist$H, H0)
  metricList$B_tr[i,1] <- tr_sta(reslist$B, B0)
  metricList$Upsilon_tr[i,1] <- tr_sta(cbind(reslist$bbeta[,1], reslist$B ), cbind(bbeta0[,1],B0))
  
  metricList$H_CCor[i,1] <- MCCor(reslist$H, H0)
  metricList$B_CCor[i,1] <- MCCor(reslist$B , B0)
  metricList$Upsilon_CCor[i,1] <- MCCor(cbind(reslist$bbeta[,1], reslist$B ), cbind(bbeta0[,1],B0))
  metricList$mu_norm[i, 1] <- norm_vec(reslist$bbeta[,1]- bbeta0[,1])
  metricList$beta_norm[i, 1] <- norm_vec(as.vector(reslist$bbeta) - as.vector(bbeta0)) 



 #######################################
 ###  Scenario 1.2                                                             #
  #######################################
  tic <- proc.time()
  reslist_no <- RR_COAP(datList$X_no, Z=matrix(1, n, 1), q= hq, rank_use =1, epsELBO = 1e-7)
  toc <- proc.time()
  time_apfactor_no <- toc[3] - tic[3]
  metricList_no$timeMat[i,1] <- time_apfactor_no
  metricList_no$H_tr[i,1] <- tr_sta(reslist_no$H, H0)
  metricList_no$B_tr[i,1] <- tr_sta(reslist_no$B, B0)
  metricList_no$Upsilon_tr[i,1] <- tr_sta(cbind(reslist_no$bbeta, reslist_no$B ), cbind(bbeta0[,1],B0))
  
  metricList_no$H_CCor[i,1] <- MCCor(reslist_no$H, H0)
  metricList_no$B_CCor[i,1] <- MCCor(reslist_no$B , B0)
  metricList_no$Upsilon_CCor[i,1] <- MCCor(cbind(reslist_no$bbeta, reslist_no$B ), cbind(bbeta0[,1],B0))
  metricList_no$mu_norm[i, 1] <- norm_vec(reslist_no$bbeta- bbeta0[,1])
  

   #######################################
  ##### Scenario 1.3                                                             #
  #######################################
  datList <- gendata_simu_lowrankai(n=n, p=p, d=d,rank0 = rank0,q=q, rho=c(1, 2,10)/2, seed = i, sigma2_eps=1, ai=20)
  X_count <- datList$X; Z <-datList$Z; bbeta0<-datList$bbeta0; H0 <- datList$H0; B0 <- datList$B0
tic <- proc.time()
  reslist <-  RR_COAP(X_count, Z=cbind(matrix(1, n, 1),Z), q= hq, rank_use = rank0, epsELBO = 1e-7)
  toc <- proc.time()
  time_apfactor <- toc[3] - tic[3]
  metricList$timeMat[i,1] <- time_apfactor
  metricList$H_tr[i,1] <- tr_sta(reslist$H, H0)
  metricList$B_tr[i,1] <- tr_sta(reslist$B, B0)
  metricList$Upsilon_tr[i,1] <- tr_sta(cbind(reslist$bbeta, reslist$B ), cbind(bbeta0[,1],B0))
  
  metricList$H_CCor[i,1] <- MCCor(reslist$H, H0)
  metricList$B_CCor[i,1] <- MCCor(reslist$B , B0)
  metricList$Upsilon_CCor[i,1] <- MCCor(cbind(reslist$bbeta[,1], reslist$B ), cbind(bbeta0[,1],B0))
  metricList$mu_norm[i, 1] <- norm_vec(reslist$bbeta[,1]-log(20)- bbeta0[,1])
  c_coap<-reslist$bbeta
  c_coap[,1]<- c_coap[,1]-log(20)
  metricList$beta_norm[i, 1] <- norm_vec(as.vector(c_coap) - as.vector(bbeta0)) 
  
  }




