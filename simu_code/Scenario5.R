setwd("D:\\LearnFiles\\Research paper\\idea\\3. PoisFactor\\Simu\\simu_code")
source("./simu_util.R")


library(COAP)
n_set <- c(200, 500, 1000, 3000, 8000, 2e4)
p <- 1000; d= 100; 
q <- 5; rank0 <- 6
N <- 10
methodNames <- c("COAP", "COAP_nocovariates", "PLNPCA", "GFM", "MRRR", "MRRR2")
n_methods <- length(methodNames)
timeArray <- array(dim=c(N, n_methods, length(n_set)))
colnames(timeArray) <- methodNames
hq <- q;
for(i_n in 1:length(n_set)){ 
  # i_n <- 1
  
  
  for(i in 1: N){
    #i <- 1
    message('i = ', i, ", i_n = ", i_n)
    # i_n <- 1
    datList <- gendata_simu_lowrank_select(n=n_set[i_n], p=p, d=d, rank0 = rank0, q=q, rho=c(3,6)/2, seed = i)
    X_count <- datList$X; Z <-datList$Z;  H0 <- datList$H0; B0 <- datList$B0; bbeta0 <- datList$bbeta0
    summary(as.vector(X_count))
    
    tic <- proc.time() ## update beta is too slow
    reslist <- RR_COAP(datList$X, Z = datList$Z, q= hq, rank_use = rank0, epsELBO = 1e-20, maxIter = 30, fast_svd = TRUE)
    toc <- proc.time()
    time_apfactor <- toc[3] - tic[3]
    MCCor(reslist$H , H0)
    MCCor(reslist$B , B0)
    norm_vec(reslist$bbeta[,1]- bbeta0[,1])

    timeArray[i,1, i_n] <- time_apfactor

    tic <- proc.time() ## The case without covariates
    reslist <- RR_COAP(datList$X, q= hq, epsELBO = 1e-20, maxIter = 30, fast_svd = TRUE)
    toc <- proc.time()
    time_apfactor2 <- toc[3] - tic[3]

    timeArray[i,2, i_n] <- time_apfactor2



    ## GFM
    library(GFM)
    tic <- proc.time()
    res_gfm <- gfm(list(X_count), types=c("poisson"), q=hq, dc_eps=1e-20, maxIter = 30)
    toc <- proc.time()
    time_gfm <- toc[3] - tic[3]
    timeArray[i,4, i_n] <- time_gfm

    ## MRRR
    tic <- proc.time()
    res_mrrr <- mrrr_run(X_count, Z, rank0 = rank0,  epsilon = 1e-20, sv.tol = 1e-50, maxIter = 30)
    toc <- proc.time()
    time_mrrr <- toc[3] - tic[3]
    timeArray[i,5, i_n] <- time_mrrr

    ## MRRR variant
    tic <- proc.time()
    res_mrrr2 <- mrrr_run(X_count, Z, rank0 = rank0, q= q,  epsilon = 1e-20, sv.tol = 1e-20, maxIter = 30)
    toc <- proc.time()
    time_mrrr <- toc[3] - tic[3]
    timeArray[i,6, i_n] <- time_mrrr
    hbeta_mrrr2 <- res_mrrr2$coef[1:d,]
    hBH_mrrr2 <- res_mrrr2$coef[(d+1): (d+n_set[i_n]),]
    svdBH <- svd(hBH_mrrr2, nu=q, nv=q)

    MCCor(svdBH$u , H0)
    MCCor(svdBH$v %*% Diag(svdBH$d[1:q]), B0)


    ## PLNPCA 
    tic <- proc.time()
    res_plnpca <- PLNPCA_run(X_count, Z[,-1], hq, maxIter=30,ftol_rel=1e-20, xtol_rel= 1e-20)
    toc <- proc.time()
    time_plnpca <- toc[3] - tic[3]
    timeArray[i,3, i_n] <- time_plnpca
    
    
  }
  timeArray[,,i_n]
  
  save(timeArray, file=paste0("p",p, "d", d, "_varyn_timeComp_method_PLNPCA.rds"))
}

save(timeArray, file=paste0("p",p, "d", d, "_varyn_timeComp_method5.rds"))





# Compare COAP and GLLVM --------------------------------------------------



gllvm_run <- function(Y, X, q){
  library(gllvm)
  if(!is.null(X)){
    nc <- ncol(X)
    colnames(X) <-  paste0("V",1:nc)
  }
  
  tic <- proc.time()
  res_gllvm <- gllvm(y=Y, X=X, family=poisson(), num.lv= q, control = list(trace=T))
  # hH_gllvm <- predictLVs.gllvm(res_gllvm)
  toc <- proc.time()
  time_gllvm <- toc[3] - tic[3]
  res_gllvm$time_used <- time_gllvm
  return(res_gllvm)
}
library(COAP)
library(gllvm)

#### Fix p=50; change n
n_set <- c(80, 150, 200, 250, 400, 600)
p <- 50; d= 50; 
q <- 5; rank0 <- 6
N <- 5
methodNames <- c("COAP","GLLVM")
n_methods <- length(methodNames)
timeArray_gllvm <- array(dim=c(N, n_methods, length(n_set)))
colnames(timeArray_gllvm) <- methodNames
hq <- q;
for(i_n in 1:length(n_set)){ 
  # i_n <- 1
  
  
  for(i in 1: N){
    #i <- 1
    message('i = ', i, ", i_n = ", i_n)
    # i_n <- 1
    datList <- gendata_simu_lowrank_select(n=n_set[i_n], p=p, d=d, rank0 = rank0, q=q, rho=c(3,6)/2, seed = i)
    X_count <- datList$X; Z <-datList$Z;  H0 <- datList$H0; B0 <- datList$B0; bbeta0 <- datList$bbeta0
    summary(as.vector(X_count))
    
    tic <- proc.time() ## 
    reslist <- RR_COAP(datList$X, Z = datList$Z, q= hq, rank_use = rank0, epsELBO = 1e-20, maxIter = 30, fast_svd = TRUE)
    toc <- proc.time()
    time_apfactor <- toc[3] - tic[3]
    timeArray_gllvm[i,1, i_n] <- time_apfactor
    
    tic <- proc.time() ## 
    try({
      reslist <- gllvm_run(datList$X, X = datList$Z, q= hq)
    },
    silent=T 
    )
    
    toc <- proc.time()
    time_gllvm <- toc[3] - tic[3]
    # tic <- proc.time() ## 
    # reslist <- gllvm_run(datList$X, X=NULL, q= hq)
    # toc <- proc.time()
    # time_gllvm_nocov <- toc[3] - tic[3]
    timeArray_gllvm[i,2, i_n] <- time_gllvm
  }
  timeArray_gllvm[,,i_n]
  
  save(timeArray_gllvm, file=paste0("p",p, "d", d, "_varyn_timeComp_method_rep5_COAPvsGLLVM.rds"))
}
save(timeArray_gllvm, n_set, file=paste0("p",p, "d", d, "_varyn_timeComp_method_rep5_COAPvsGLLVM.rds"))

#### Fix n=50; change p
library(COAP)
library(gllvm)
p_set <-  c(80, 150, 200, 250, 400, 600)
n <- 80; d= 50; 
q <- 5; rank0 <- 6
N <- 5
methodNames <- c("COAP", "GLLVM")
n_methods <- length(methodNames)
timeArray_gllvm <- array(dim=c(N, n_methods, length(p_set)))
colnames(timeArray_gllvm) <- methodNames
hq <- q;
for(i_n in 1:length(p_set)){ 
  # i_n <- 1
  
  
  for(i in 1: N){
    #i <- 1
    message('i = ', i, ", i_n = ", i_n)
    # i_n <- 1
    datList <- gendata_simu_lowrank_select(n=n, p=p_set[i_n], d=d, rank0 = rank0, q=q, rho=c(3,6)/2, seed = i)
    X_count <- datList$X; Z <-datList$Z;  H0 <- datList$H0; B0 <- datList$B0; bbeta0 <- datList$bbeta0
    summary(as.vector(X_count))
    
    tic <- proc.time() ## update beta is too slow
    reslist <- RR_COAP(datList$X, Z = datList$Z, q= hq, rank_use = rank0, epsELBO = 1e-20, maxIter = 30, fast_svd = TRUE)
    toc <- proc.time()
    time_apfactor <- toc[3] - tic[3]
    MCCor(reslist$H , H0)
    MCCor(reslist$B , B0)
    norm_vec(reslist$bbeta[,1]- bbeta0[,1])
    
    timeArray_gllvm[i,1, i_n] <- time_apfactor
    
    tic <- proc.time() ## 
    try({
      reslist <- gllvm_run(datList$X, X = datList$Z, q= hq)
    },
    silent=T 
    )
    toc <- proc.time()
    time_gllvm <- toc[3] - tic[3]
    timeArray_gllvm[i,2, i_n] <- time_gllvm
    
    
  }
  timeArray_gllvm[,,i_n]
  
  save(timeArray_gllvm, file=paste0("n",n, "d", d, "_varyp_timeComp_COAPvsGLLVM.rds"))
}




