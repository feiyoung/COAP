# 
# datList <- gendata(seed=20, n=100, p=100)
# colSums(datList$X)
# X_count <- datList$X; Z <- datList$Z; H <- datList$H0; B <- datList$B0
# n <- nrow(X_count); p <- ncol(X_count); d <- ncol(Z); q <- ncol(H)
# 
# 
# 
# Mu_y_int = log(1+ X_count)#matrix(1, n, p);
# S_y_int = matrix(1, n, p);
# a <- rep(1, n)
# set.seed(1)
# fit_lfm <- GFM::Factorm(Mu_y_int, q=q)
# #B_int = matrix(rnorm(p*q), p, q)*0.1;
# # Mu_h_int <-  matrix(rnorm(n*q), n, q)#matrix(0, n, q)
# B_int <- fit_lfm$hB
# Mu_h_int <-  fit_lfm$hH
# bbeta_int <- matrix(0, p, d)
# S_h_int <- array(dim=c(q,q, n)) ## fast computation version that does not update Mu_y and S_y.
# for(i in 1:n) S_h_int[,,i] <- diag(rep(1, q))
# Sigma_h_int <- diag(rep(1, q))
# invLambda_int = rep(1, p); ## the default p is 100, the default q is 15
# tau_int = matrix(colMeans(Mu_y_int), p, 1);
# 
# 
# 
# # reslist <- APoisFactor(X_count, a, Z, Mu_y_int, S_y_int, invLambda_int, B_int,
# #                        bbeta_int, Mu_h_int, S_h_int, Sigma_h_int, epsELBO=1e-9,
# #                        maxIter=30, verbose=T, fast_version=1)
# #
# # reslist <- APoisFactor(X_count, a, Z=cbind(1, Z), Mu_y_int, S_y_int, invLambda_int, B_int,
# #                        cbind(tau_int,bbeta_int), Mu_h_int, S_h_int, Sigma_h_int, epsELBO=1e-9,
# #                        maxIter=30, verbose=T, fast_version=2)
# hq <- 6
# reslist <- APoissonFactor(X_count, Z=cbind(1,Z), q=hq)
# reslist <- APoissonFactor(X_count, Z=cbind(1,Z), q=hq, fast_version = "Rough_Approx")
# plot(reslist$ELBO_seq[-1], type='o')
# reslist$B[1:4,1:5]
# bbeta0 <- cbind( datList$mu0, datList$bbeta0)
# reslist$bbeta[1:5,]; bbeta0[1:5,]
# cancor(reslist$H, H)$cor
# cancor(reslist$B, B)$cor
# 
# str(reslist)
# 
# norm(reslist$bbeta-bbeta0, 'F')
# 
# reslist2 <-  APoissonFactor(X_count, Z=cbind(1,Z), q=hq)
# 
# norm(reslist2$bbeta-bbeta0, 'F')
# colMeans(reslist2$H)
# cancor(datList$H0, reslist2$H)$cor
# 
# fit_gfm <- GFM::gfm(X_count, group=rep(1,p), type='poisson', q= q, maxIter = 2)
# cancor(fit_gfm$hH, H)$cor
# fit_lfm <- GFM::Factorm(log(1+X_count), q=q)
# fit_lfm <- GFM::Factorm(X_count, q=q)
# cancor(fit_lfm$hH, H)$cor
# 
# 
# ### There is no covariates
# reslist3 <- APoissonFactor(X_count, Z=NULL, q=hq)
# colMeans(reslist3$H)
# norm(reslist3$tau - datList$mu0, 'F')
# cancor(datList$H0, reslist3$H)$cor
# 
# # Test running time for each version --------------------------------------
# datList <- gendata(seed=2, n=100, p=100)
# colSums(datList$X)
# X_count <- datList$X; Z <- datList$Z; H <- datList$H0; B <- datList$B0
# n <- nrow(X_count); p <- ncol(X_count); d <- ncol(Z); q <- ncol(H)
# version_vec <- c("Laplace_Taylor",  "Rough_Approx")
# time_vec <- rep(NA, length(version_vec))
# for(i in seq_along(version_vec)){
#   # i <- 1
#   cat(" version = ", version_vec[i], '\n')
#   tic <- proc.time()
#   reslist <- RR_COAP(X_count, Z=Z, q=q, fast_version = version_vec[i], maxIter = 20, epsELBO = 1e-20)
#   cat(cancor(reslist$H, H)$cor, '\n')
#   cat(cancor(reslist$B, B)$cor, '\n')
#   toc <- proc.time()
#   time_vec[i] <- toc[3] - tic[3]
# 
# }
# 
# reslist <- RR_COAP(X_count, q=q, fast_version = version_vec[i], maxIter = 20, epsELBO = 1e-20)
# 
# version_vec <- c('Factorm', "approxPCA")
# time_vec <- rep(NA, length(version_vec))
# for(i in seq_along(version_vec)){
#   # i <- 3
#   cat(" version = ", version_vec[i], '\n')
#   tic <- proc.time()
#   reslist <- APoissonFactor(X_count, Z=Z, q=q, initial = version_vec[i], maxIter = 20, epsELBO = 1e-20)
#   cat(cancor(reslist$H, H)$cor, '\n')
#   cat(cancor(reslist$B, B)$cor, '\n')
#   toc <- proc.time()
#   time_vec[i] <- toc[3] - tic[3]
# 
# }
# 
# 
# 
# # High-dimensional covariates ---------------------------------------------
# rank0 <- 6; q = 5; d= 50;
# datList <- gendata_simu_lowrank_select(seed = 3,n=1000, p=1000, d= d, rank0 = rank0, q= q, rho=c(3, 6)/5)
# X_count <- datList$X
# Z <- datList$Z
# 
# reslist <- RR_COAP(datList$X, Z = datList$Z, rank_use =  rank0, q= q, verbose = T)
# 
# 
# library(microbenchmark)
# mbm <- microbenchmark("fast" = { reslist <- RR_COAP(X_count, Z =  Z, rank_use =  rank0, fast_svd= T)
# },
# "slow" = {
#   reslist <- RR_COAP(X_count, Z = Z, rank_use = 3, fast_svd= F)
# },times = 2); mbm
# 
# 
# ## Compare with GFM
# reslist <- RR_COAP(X_count,  rank_use = 3, epsELBO = 1e-4, fast_svd= T)
# reslist <- RR_COAP(X_count,  rank_use = 3, q=1, epsELBO = 1e-4, fast_svd= T)
# 
# mbm <- microbenchmark("fast" = {
#   reslist <- RR_COAP(X_count,  rank_use = rank0, epsELBO = 1e-20, maxIter=30, fast_svd= T)
# },
# "slow" = {
#   reslist <-GFM::gfm(list(X_count), types="poisson", q=15, dc_eps = 1e-20, maxIter = 30)
# },times = 2); mbm
# 
# 
# ## Compare with MRRR
# mrrr_run <- function(Y, X, rank0, family=list(poisson()),
#                      familygroup=rep(1,ncol(Y)), epsilon = 1e-4, sv.tol = 1e-2, maxIter = 30, trace=TRUE){
#   require(rrpack)
# 
#   n <- nrow(Y); p <- ncol(Y)
# 
# 
#   svdX0d1 <- svd(X)$d[1]
#   init1 = list(kappaC0 = svdX0d1 * 5)
#   offset = NULL
#   control = list(epsilon = epsilon, sv.tol = sv.tol, maxit = maxIter,
#                  trace = trace, gammaC0 = 1.1, plot.cv = TRUE,
#                  conv.obj = TRUE)
#   fit.mrrr <- mrrr(Y=Y, X=X[,-1], family = family, familygroup = familygroup,
#                    penstr = list(penaltySVD = "rankCon", lambdaSVD = 0.1),
#                    control = control, init = init1, maxrank = rank0)
# 
#   return(fit.mrrr)
# }
# 
# mbm <- microbenchmark("fast" = {
#   reslist <- RR_COAP(X_count,  rank_use = rank0, q=q, epsELBO = 1e-4, fast_svd= T)
# },
# "slow" = {
#   res_mrrr <- mrrr_run(Y=X_count, X = Z, rank0 = rank0)
# },times = 2); mbm
# 
# 
# mbm <- microbenchmark("fast" = {
#   reslist <- RR_COAP(X_count,  rank_use = rank0, q=q, epsELBO = 1e-20, maxIter = 30 ,fast_svd= T)
# },
# "slow" = {
#   res_mrrr <- mrrr_run(Y=X_count, X = Z, rank0=10, epsilon = 1e-20, sv.tol=1e-20, maxIter = 30) # , rank0 = rank0
# },times = 2); mbm
# 
# 
