gendata<-function (seed = 1, n = 300, p = 50, d=3, q = 6, rho = 1, sigma2_eps=0.1){
  library(MASS)
  cor.mat<-function (p, rho, type = "toeplitz") {
    if (p == 1) 
      return(matrix(1, 1, 1))
    mat <- diag(p)
    if (type == "toeplitz") {
      for (i in 2:p) {
        for (j in 1:i) {
          mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
        }
      }
    }
    if (type == "identity") {
      mat[mat == 0] <- rho
    }
    return(mat)
  }
  Diag<-function (vec){
    q <- length(vec)
    if (q > 1) {
      y <- diag(vec)
    }
    else {
      y <- matrix(vec, 1, 1)
    }
    return(y)
  }
  if (length(rho) == 1) 
    rho <- c(rho, 1.5)
  factor_term <- rho[2]
  set.seed(1)
  bbeta0 <- matrix(rnorm(p*d), p, d) # Fixed bbeta0
  Ztmp <- matrix(rnorm(p * q), p, q)
  B <- qr(Ztmp)
  eigvs <- sqrt(sort(eigen(t(Ztmp) %*% Ztmp)$values, decreasing = T))
  B1 <- qr.Q(B) %*% Diag(sqrt(eigvs))
  B0 <- rho[1] * B1 %*% Diag(sign(B1[1, ])) ## Fixed B0 and mu0 for each repeat.
  mu0 <- 0.4 * rnorm(p)
  Bm0 <- cbind(mu0, B0)
  set.seed(seed)
  H <- mvrnorm(n, mu = rep(0, q), cor.mat(q, 0.5))
  svdH <- svd(cov(H))
  H0 <- scale(H, scale = F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*% 
    svdH$v
  
  g1 <- 1:p
  B0[g1, ] <- B0[g1, ]/max(B0[g1, ]) * factor_term ## scale
  
  Z <- MASS::mvrnorm(n, mu=rep(0, d), Sigma = cor.mat(d, rho=0.5))
  epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
  
  mu <- exp(Z %*% t(bbeta0) + H0 %*% t(B0) + matrix(mu0, n, p, byrow = T)+ epsi) 
  
  X <- matrix(rpois(n * p, lambda = mu), n, p)
  
  return(list(X = X, Z=Z, bbeta0=bbeta0, B0 = B0, H0 = H0, mu0 = mu0))
}


# gendata_simu_lowrank_select <-function (seed = 1, n = 300, p = 50, d=20, q = 6, rank0=3, rho = c(1.5, 1), sigma2_eps=0.1){
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


Diag <- function (vec) {
  q <- length(vec)
  if (q > 1) {
    y <- diag(vec)
  }
  else {
    y <- matrix(vec, 1, 1)
  }
  return(y)
}
add_identifiability <- function(H, B){
  q <- ncol(H); n <- nrow(H)
  svdHB <- svd(H %*% t(B), nu=q, nv = q)
  signB1 <- sign(svdHB$v[1,])
  H <- sqrt(n) * svdHB$u %*% Diag(signB1)
  
  B <- svdHB$v %*% Diag(svdHB$d[1:q]*signB1) / sqrt(n)
  
  return(list(H=H, B=B))
}

RR_COAP <- function(X_count, multiFac=rep(1, nrow(X_count)), Z=matrix(1, nrow(X_count),1),rank_use=5,
                    q=15, epsELBO=1e-5, 
                    maxIter=30, verbose=TRUE, fast_version=c("Laplace_Taylor", "Rough_Approx"),
                    joint_opt_beta=FALSE, fast_svd=TRUE){
  
  # Z=NULL; q=15; epsELBO=1e-6;  maxIter=10; verbose=TRUE
  
  get_initials <- function(X, q){
    require(irlba)
    n <- nrow(X); p <- ncol(X)
    mu <- colMeans(X)
    X <- X - matrix(mu, nrow=n, ncol=p, byrow=TRUE)
    svdX  <- irlba(A =X, nv = q)
    PCs <- sqrt(n) * svdX$u
    loadings <- svdX$v %*% Diag(svdX$d[1:q]) / sqrt(n)
    dX <- PCs %*% t(loadings) - X
    Lam_vec <- colSums(dX^2)/n
    return(list(hH = PCs, hB = loadings, hmu=mu,sigma2vec = Lam_vec))
    
  }
  
  fast_version <- match.arg(fast_version)
  fast_version <- switch (fast_version,
                          Laplace_Taylor = 1,
                          Rough_Approx = 0
  )
  
  
  
  message("Calculate initial values...")
  n <- nrow(X_count); p <- ncol(X_count); 
  if(any(Z[,1]!=1)) warning("The first column of covariates Z is not a full-one column vector, so it will fit a model without intercept!")
  Mu_y_int = log(1+ X_count)#matrix(1, n, p);
  S_y_int = matrix(1, n, p);
  a <- multiFac
  set.seed(1)
  fit_approxPCA <- get_initials(Mu_y_int, q=q)
  B_int <- fit_approxPCA$hB
  Mu_h_int <-  fit_approxPCA$hH
  
  
  S_h_int <- diag(rep(1, q))
  Sigma_h_int <- diag(rep(1, q))
  invLambda_int = rep(1, p); ## the default p is 100, the default q is 15
  
  d <- ncol(Z)
  rank_use <- min(rank_use, d)
  
  bbeta_int <- matrix(0, p, d)
  bbeta_int[,1] <- colMeans(Mu_y_int)
  
  reslist <- Reduced_Rank_COAP(X_count, a, Z, rank_use, Mu_y_int, S_y_int, invLambda_int, B_int, 
                               bbeta_int,  Mu_h_int, S_h_int, Sigma_h_int, epsELBO=epsELBO, 
                               maxIter=maxIter, verbose=verbose,  fast_version=fast_version, 
                               sep_opt_beta=!joint_opt_beta, fast_svd=fast_svd) 
  

  return(reslist)
  
}
