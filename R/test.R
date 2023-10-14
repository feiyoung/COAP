# 
# ## Check Cpp function ------------------------------------------------------
# 
# 
# # datList <- gendata(seed=20, n=100, p=100)
# # colSums(datList$X)
# # X_count <- datList$X; Z <- datList$Z; H <- datList$H0; B <- datList$B0
# # n <- nrow(X_count); p <- ncol(X_count); d <- ncol(Z); q <- ncol(H)
# #
# #
# #
# # Mu_y_int = log(1+ X_count)#matrix(1, n, p);
# # S_y_int = matrix(1, n, p);
# # a <- rep(1, n)
# # set.seed(1)
# # fit_lfm <- GFM::Factorm(Mu_y_int, q=q)
# # #B_int = matrix(rnorm(p*q), p, q)*0.1;
# # # Mu_h_int <-  matrix(rnorm(n*q), n, q)#matrix(0, n, q)
# # B_int <- fit_lfm$hB
# # Mu_h_int <-  fit_lfm$hH
# # bbeta_int <- matrix(0, p, d)
# # S_h_int <- array(dim=c(q,q, n)) ## fast computation version that does not update Mu_y and S_y.
# # for(i in 1:n) S_h_int[,,i] <- diag(rep(1, q))
# # Sigma_h_int <- diag(rep(1, q))
# # invLambda_int = rep(1, p); ## the default p is 100, the default q is 15
# # tau_int = matrix(colMeans(Mu_y_int), p, 1);
# #
# #
# #
# # # reslist <- APoisFactor(X_count, a, Z, Mu_y_int, S_y_int, invLambda_int, B_int,
# # #                        bbeta_int, Mu_h_int, S_h_int, Sigma_h_int, epsELBO=1e-9,
# # #                        maxIter=30, verbose=T, fast_version=1)
# # #
# # # reslist <- APoisFactor(X_count, a, Z=cbind(1, Z), Mu_y_int, S_y_int, invLambda_int, B_int,
# # #                        cbind(tau_int,bbeta_int), Mu_h_int, S_h_int, Sigma_h_int, epsELBO=1e-9,
# # #                        maxIter=30, verbose=T, fast_version=2)
# 
# 
# 
# # Compare with GFM --------------------------------------------------------
# rank0 <- 6; q = 5; d= 50
# datList <- gendata_simu(seed = 1, n=300, p=300, d= d, rank0 = rank0, q= q, rho=c(1, 4),
#                         sigma2_eps = 3)
# X_count <- datList$X; Z <- datList$Z
# H <- datList$H0; B <- datList$B0
# hq <- 5; hr <- 6
# system.time(
#   reslist <- RR_COAP(X_count, Z= Z, q=hq, rank_use= hr)
# )
# system.time(
#   reslist <- RR_COAP(X_count, Z= Z, q=hq, fast_svd = F)
# )
# 
# # reslist <- RR_COAP(X_count, Z=Z, q=hq, fast_version = "Rough_Approx")
# # plot(reslist$ELBO_seq[-1], type='o')
# # reslist$B[1:4,1:5]
# bbeta0 <- cbind( datList$mu0, datList$bbeta0)
# # reslist$bbeta[1:5,1:4]; bbeta0[1:5,1:4]
# GFM::measurefun(reslist$H, H)
# GFM::measurefun(reslist$B, B)
# 
# str(reslist)
# norm_vec <- function(x) sqrt(sum(x^2/ length(x)))
# norm_vec(reslist$bbeta-bbeta0)
# 
# 
# fit_gfm <- GFM::gfm(list(X_count),  type='poisson', q= q)
# GFM::measurefun(fit_gfm$hH, H)
# GFM::measurefun(fit_gfm$hB, B)
# norm_vec(fit_gfm$hmu- bbeta0[,1])
# fit_lfm <- GFM::Factorm(X_count, q=q)
# GFM::measurefun(fit_lfm$hH, H)
# GFM::measurefun(fit_lfm$hB, B)
# 
# 
# 
# 
# 
# ### There is no covariates
# reslist3 <- RR_COAP(X_count, q=hq)
# cancor(datList$H0, reslist3$H)$cor
# 
# 
# 
# 
# 
# 
# gendata_simu_lowrank_select <- function (seed = 1, n = 300, p = 50, d=20, q = 6, rank0=3, rho = c(1.5, 1), sigma2_eps=0.1){
#   require(MASS)
#   if(rank0<=1) stop("rank0 must be greater than 1!")
#   cor.mat<-function (p, rho, type = "toeplitz") {
#     if (p == 1)
#       return(matrix(1, 1, 1))
#     mat <- diag(p)
#     if (type == "toeplitz") {
#       for (i in 2:p) {
#         for (j in 1:i) {
#           mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
#         }
#       }
#     }
#     if (type == "identity") {
#       mat[mat == 0] <- rho
#     }
#     return(mat)
#   }
#   Diag<-function (vec){
#     q <- length(vec)
#     if (q > 1) {
#       y <- diag(vec)
#     }
#     else {
#       y <- matrix(vec, 1, 1)
#     }
#     return(y)
#   }
# 
#   if(length(rho)<2) stop("rho must be a numeric vector of length 2!")
# 
# 
#   factor_term <- rho[1]
#   factor_term_z <- rho[2]
#   set.seed(1) # Fixed bbeta0
#   #bbeta0 <- matrix(rnorm(p*d), p, d)
#   rank_true <- rank0 - 1
#   bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p)) / p *4 * factor_term_z
#   Ztmp <- matrix(rnorm(p * q), p, q)
#   B <- qr(Ztmp)
#   eigvs <- sqrt(sort(eigen(t(Ztmp) %*% Ztmp)$values, decreasing = T))
#   B1 <- qr.Q(B) %*% Diag(sqrt(eigvs))
#   B0 <- B1 %*% Diag(sign(B1[1, ])) ## Fixed B0 and mu0 for each repeat.
#   # mu0 <- rnorm(p) * factor_term
#   # Bm0 <- cbind(mu0, B0)
#   set.seed(seed)
# 
#   if(d<2) stop("d must be greater than 1!")
#   Z <- MASS::mvrnorm(n, mu=rep(0, d-1), Sigma = cor.mat(d-1, rho=0.5))
#   Z <- cbind(1, Z)
#   epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
#   H <- mvrnorm(n, mu = rep(0, q), cor.mat(q, 0.5))
#   H <- residuals(lm(H~Z))
#   svdH <- svd(cov(H))
#   H0 <- scale(H, scale = F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*%
#     svdH$v
# 
#   g1 <- 1:p
#   B0[g1, ] <- B0[g1, ]/max(B0[g1, ]) * factor_term ## scale
# 
# 
#   mu <- exp(Z %*% t(bbeta0) + H0 %*% t(B0) + epsi) # + matrix(mu0, n, p, byrow = T)
# 
#   X <- matrix(rpois(n * p, lambda = mu), n, p)
# 
#   return(list(X = X, Z=Z, bbeta0=bbeta0, B0 = B0, H0 = H0,  rank=rank0, q=q))
# }
# 
# rank0 <- 6; q = 5; d= 50
# datList <- gendata_simu_lowrank_select(seed = 1, n=300, p=300, d= d, rank0 = rank0, q= q, rho=c(3, 6)/2,
#                         sigma2_eps = 3)
# X_count <- datList$X; Z <- datList$Z
# H <- datList$H0; B <- datList$B0
# hq <- 5; hr <- 6
# system.time(
#   reslist <- RR_COAP(X_count, Z= Z, q=hq, rank_use= hr)
# )
# 
# GFM::measurefun(reslist$H, H)
# GFM::measurefun(reslist$B, B)
# 
# str(reslist)
# norm_vec <- function(x) sqrt(sum(x^2/ length(x)))
# norm_vec(reslist$bbeta-bbeta0)
# 
# 
# fit_gfm <- GFM::gfm(list(X_count),  type='poisson', q= q)
# GFM::measurefun(fit_gfm$hH, H)
# GFM::measurefun(fit_gfm$hB, B)
# fit_lfm <- GFM::Factorm(X_count, q=q)
# GFM::measurefun(fit_lfm$hH, H)
# GFM::measurefun(fit_lfm$hB, B)
# 
# # ## Tune the signals to make the rank and number of factors are i --------
# 
# 
# ## rho=c(3, 9) is good.
# datList <- gendata_simu(seed = 1,n=150, p=200, d=50, rank0 = 6, q=5, rho=c(3,5), sigma2_eps = 1)
# q_max <- 15
# d <- ncol(datList$Z)
# 
# reslist <- RR_COAP(datList$X, Z = datList$Z, rank_use = floor(d/2), q= q_max, verbose = T,joint_opt_beta=F)
# 
# threshold <- 0.1
# svalues <- svd(reslist$bbeta)$d
# # cumsum(svalues)/ sum(svalues)
# 
# par(mfrow=c(2,1))
# svalues <- svalues[svalues>threshold]
# ratio1 <- svalues[-length(svalues)] / svalues[-1]
# 
# dat1 <- data.frame(Ratio=ratio1, r=1:length(ratio1))
# library(ggplot2)
# p1 <- ggplot(data=dat1, aes(x=r, y=Ratio)) + geom_line(linewidth=0.8) +geom_point(size=1.8)+
#   scale_x_continuous(breaks=seq(3, length(ratio1), by=3)) +theme_bw(base_size = 20)
# plot(ratio1, type='o', xlab='r', ylab=paste0('ratio of eigenvalue of ',expression(beta)))
# abline(v=6, col='red')
# which.max(ratio1[-length(ratio1)])
# svalues <- svd(reslist$B)$d
# svalues <- svalues[svalues>1e-2]
# ratio_fac <- svalues[-length(svalues)] / svalues[-1]
# dat2 <- data.frame(Ratio=ratio_fac, q=1:length(ratio_fac))
# p2 <- ggplot(data=dat2, aes(x=q, y=Ratio)) + geom_line(linewidth=0.8) +geom_point(size=1.5) +
#   theme_bw(base_size=20)
# p12 <- SRTpipeline::drawFigs(pList=list(p1,p2), layout.dim = c(1,2), legend.position='none')
# p12
# 
# ### Use the functions
# res1 <- selectParams(X_count=datList$X, Z=datList$Z)
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
# 
# rank0 <- 6; q = 5; d= 50;
# datList <- gendata_simu(seed = 3,n=1000, p=1000, d= d, rank0 = rank0, q= q, rho=c(3, 6)/5)
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
