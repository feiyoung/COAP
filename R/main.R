# generate man files
# devtools::document()
# R CMD check --as-cran ProFAST_1.2.tar.gz
## usethis::use_data(dat_r2_mac)
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("COAPsimu")
# pkgdown::build_article("ProFASTdlpfc2")

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('html_document'))

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('pdf_document'), clean = F)

#' Generate simulated data
#' @description Generate simulated data from covariate-augmented Poisson factor models
#' @param seed a postive integer, the random seed for reproducibility of data generation process.
#' @param n a postive integer, specify the sample size. 
#' @param p a postive integer, specify the dimension of count variables.
#' @param d a postive integer,  specify the dimension of covariate matrix.
#' @param q a postive integer,  specify the number of factors.
#' @param rank0 a postive integer, specify the rank of the coefficient matrix.
#' @param rho a numeric vector with length 2 and positive elements, specify the signal strength of regression coefficient and loading matrix, respectively. 
#' @param sigma2_eps a positive real, the variance of overdispersion error.
#' @return return a list including the following components: (1) X, the high-dimensional count matrix; (2) Z, the high-dimensional covriate matrix; (3) bbeta0, the low-rank large coefficient matrix; (4) B0, the loading matrix; (5) H0, the factor matrix; (6) rank: the true rank of bbeta0; (7) q: the true number of factors.
#' @details None
#' @seealso \code{\link{RR_COAP}}
#' @references None
#' @export
#' @importFrom  MASS mvrnorm
#'
#'

gendata_simu <-function (seed = 1, n = 300, p = 50, d=20, q = 6, rank0=3, rho = c(1.5, 1), sigma2_eps=0.1){
  require(MASS)
  if(rank0<=1) stop("rank0 must be greater than 1!")
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
  
  if(length(rho)<2) stop("rho must be a numeric vector of length 2!")
  
  
  factor_term <- rho[1]
  factor_term_z <- rho[2]
  set.seed(1) # Fixed bbeta0
  #bbeta0 <- matrix(rnorm(p*d), p, d) 
  rank_true <- rank0 - 1
  bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p)) / p *4 * factor_term_z
  Ztmp <- matrix(rnorm(p * q), p, q)
  B <- qr(Ztmp)
  eigvs <- sqrt(sort(eigen(t(Ztmp) %*% Ztmp)$values, decreasing = T))
  B1 <- qr.Q(B) %*% Diag(sqrt(eigvs))
  B0 <- B1 %*% Diag(sign(B1[1, ])) ## Fixed B0 and mu0 for each repeat.
  # mu0 <- rnorm(p) * factor_term
  # Bm0 <- cbind(mu0, B0)
  set.seed(seed)
  
  if(d<2) stop("d must be greater than 1!")
  Z <- MASS::mvrnorm(n, mu=rep(0, d-1), Sigma = cor.mat(d-1, rho=0.5))
  Z <- cbind(1, Z)
  epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
  H <- mvrnorm(n, mu = rep(0, q), cor.mat(q, 0.5))
  H <- residuals(lm(H~Z))
  svdH <- svd(cov(H))
  H0 <- scale(H, scale = F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*% 
    svdH$v
  
  g1 <- 1:p
  B0[g1, ] <- B0[g1, ]/max(B0[g1, ]) * factor_term ## scale
  
  
  mu <- exp(Z %*% t(bbeta0) + H0 %*% t(B0) + epsi) # + matrix(mu0, n, p, byrow = T) 
  
  X <- matrix(rpois(n * p, lambda = mu), n, p)
  
  return(list(X = X, Z=Z, bbeta0=bbeta0, B0 = B0, H0 = H0,  rank=rank0, q=q))
}


# gendata_simu <-function (seed = 1, n = 300, p = 50, d=20, 
#                                         q = 6, rank0=3, rho = c(1.5, 1), sigma2_eps=0.1){
#   # require(MASS)
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
#   B0 <- B1 %*% Diag(sign(B1[1, ])) ## 
#   set.seed(seed)
# 
#   if(d<2) stop("d must be greater than 1!")
#   cov_S <-  cor.mat(d-1, rho=0.5)
#   zeroMat <- matrix(0, d-1, q)
#   cov_all <- rbind(cbind(cov_S, zeroMat), cbind(t(zeroMat), Diag(rep(1, q)) ) )
#   Z_all <- MASS::mvrnorm(n, mu=rep(0, q+d-1), Sigma = cov_all)
#   Z <- cbind(1, Z_all[,1:(d-1)])
#   epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
#   # H <- mvrnorm(n, mu = rep(0, q), cor.mat(q, 0.5))
#   # H <- residuals(lm(H~Z))
#   # svdH <- svd(cov(H))
#   # H0 <- scale(H, scale = F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*%svdH$v
#   H0 <- Z_all[, d: (q+d-1)]
#   # H0 <- scale(H0, scale=F)
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


#' Fit the COAP model
#' @description Fit the covariate-augmented overdispersed Poisson factor model
#' @param X_count a count matrix, the observed count matrix.
#' @param multiFac an optional vector, the normalization factor for each unit; default as full-one vector.
#' @param Z an optional matrix, the covariate matrix; default as a full-one column vector if there is no additional covariates.
#' @param rank_use an optional integer, specify the rank of the regression coefficient matrix; default as 5.
#' @param q an optional string, specify the number of factors; default as 15.
#' @param epsELBO  an optional positive vlaue, tolerance of relative variation rate of the envidence lower bound value, defualt as '1e-5'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param joint_opt_beta a logical value, whether use the joint optimization method to update bbeta. The default is \code{FALSE}, which means using the separate optimization method.
#' @param fast_svd a logical value, whether use the fast SVD algorithm in the update of bbeta; default is \code{TRUE}.
#' @return return a list including the following components: (1) H, the predicted factor matrix; (2) B, the estimated loading matrix; (3) bbeta, the estimated low-rank large coefficient matrix; (4) invLambda, the inverse of the estimated variances of error; (5) H0, the factor matrix; (6) ELBO: the ELBO value when algorithm stops; (7) ELBO_seq: the sequence of ELBO values.
#' @details None
#' @seealso \code{\link{RR_COAP}}
#' @references None
#' @export
#' @useDynLib COAP, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#'
RR_COAP <- function(X_count, multiFac=rep(1, nrow(X_count)), Z=matrix(1, nrow(X_count),1),
                    rank_use=5, q=15, epsELBO=1e-5, 
                    maxIter=30, verbose=TRUE,
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
  
  invLambda_int = rep(1, p); 
  
  d <- ncol(Z)
  rank_use <- min(rank_use, d)
  
  bbeta_int <- matrix(0, p, d)
  bbeta_int[,1] <- colMeans(Mu_y_int)
  
  reslist <- Reduced_Rank_COAP(X_count, a, Z, rank_use, Mu_y_int, S_y_int, invLambda_int, B_int, 
                               bbeta_int,  Mu_h_int, epsELBO=epsELBO, 
                               maxIter=maxIter, verbose=verbose, 
                               sep_opt_beta=!joint_opt_beta, fast_svd=fast_svd) 
  
  reslist$ELBO_seq <- reslist$ELBO_seq[-1]
  return(reslist)
}



#' Select the parameters in COAP models
#' @description Select the number of factors and the rank of coefficient matrix in the covariate-augmented overdispersed Poisson factor model
#' @param X_count a count matrix, the observed count matrix.
#' @param multiFac an optional vector, the normalization factor for each unit; default as full-one vector.
#' @param Z an optional matrix, the covariate matrix; default as a full-one column vector if there is no additional covariates.
#' @param q_max an optional string, specify the upper bound for the number of factors; default as 15.
#' @param r_max an optional integer, specify the upper bound for the rank of the regression coefficient matrix; default as 24.
#' @param threshold  an optional 2-dimensional positive vector, specify the the thresholds that filters the singular values of beta and B, respectively.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param joint_opt_beta a logical value, whether use the joint optimization method to update bbeta. The default is \code{FALSE}, which means using the separate optimization method.
#' @param fast_svd a logical value, whether use the fast SVD algorithm in the update of bbeta; default is \code{TRUE}.
#' @return return a named vector with names `hr` and `hq`, the estimated rank and number of factors.
#' @details The threshold is to filter the singular values with  low signal, to assist the identification of underlying model structure.
#' @seealso \code{\link{RR_COAP}}
#' @references None
#' @export
#'
#'

selectParams <- function(X_count, Z, multiFac=rep(1, nrow(X_count)), 
                         q_max=15, r_max=24,
                         threshold=c(1e-1, 1e-2),verbose=TRUE, ...){
  
  reslist <- RR_COAP(X_count, Z = datList$Z,multiFac=multiFac, rank_use = r_max, q= q_max,
                     verbose=verbose,...)
  
  thre1 <- threshold[1]
  beta_svalues <- svd(reslist$bbeta)$d
  beta_svalues <- beta_svalues[beta_svalues>thre1]
  ratio1 <- beta_svalues[-length(beta_svalues)] / beta_svalues[-1]
  hr <- which.max(ratio1[-length(ratio1)])
  
  
  thre2 <- threshold[2]
  B_svalues <- svd(reslist$B)$d
  B_svalues <- B_svalues[B_svalues>thre2]
  ratio_fac <- B_svalues[-length(B_svalues)] / B_svalues[-1]
  hq <- which.max(ratio_fac)
  
  return(c(hr=hr, hq=hq))
}

