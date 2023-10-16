gendata_simu_lowrankai <-function (seed = 1, n = 300, p = 50, d=20, q = 6, rank0=3, rho = 1, sigma2_eps=0.1,ai=20){
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
  if (length(rho) < 3) 
    rho <- c(rho, 1.5, 1)
  factor_term <- rho[2]
  factor_term_z <- rho[3]
  set.seed(1) # Fixed bbeta0
  #bbeta0 <- matrix(rnorm(p*d), p, d) 
  rank_true <- rank0 - 1
  bbeta0 <- t(matrix(rnorm(d*rank_true), d, rank_true) %*% matrix(rnorm(rank_true* p), rank_true, p)) / p * factor_term_z
  Ztmp <- matrix(rnorm(p * q), p, q)
  B <- qr(Ztmp)
  eigvs <- sqrt(sort(eigen(t(Ztmp) %*% Ztmp)$values, decreasing = T))
  B1 <- qr.Q(B) %*% Diag(sqrt(eigvs))
  B0 <- rho[1] * B1 %*% Diag(sign(B1[1, ])) ## Fixed B0 and mu0 for each repeat.
  mu0 <- rnorm(p) * factor_term
  Bm0 <- cbind(mu0, B0)
  set.seed(seed)
  
  Z <- MASS::mvrnorm(n, mu=rep(0, d), Sigma = cor.mat(d, rho=0.5))
  epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
  H <- mvrnorm(n, mu = rep(0, q), cor.mat(q, 0.5))
  H <- residuals(lm(H~Z))
  svdH <- svd(cov(H))
  H0 <- scale(H, scale = F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*% 
    svdH$v
  
  g1 <- 1:p
  B0[g1, ] <- B0[g1, ]/max(B0[g1, ]) * factor_term ## scale
  
  mu <- ai*exp(Z %*% t(bbeta0) + H0 %*% t(B0) + matrix(mu0, n, p, byrow = T) + epsi)
  X <- matrix(rpois(n * p, lambda = mu), n, p)

    mu_no <- ai*exp( H0 %*% t(B0) + matrix(mu0, n, p, byrow = T) + epsi)
  X_no <- matrix(rpois(n * p, lambda = mu_no), n, p)
  return(list(X = X,X_no=X_no, Z=Z, bbeta0=cbind(mu0, bbeta0), B0 = B0, H0 = H0,  rank=rank0, q=q))
}


norm_vec <- function(x) sqrt(sum(x^2/ length(x)))

gendata_simu_scnario1 <-function (seed = 1, n = 300, p = 50, q = 6, rho = 1, sigma2_eps=0.1){
  require(MASS)
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
  Ztmp <- matrix(rnorm(p * q), p, q)
  B <- qr(Ztmp)
  eigvs <- sqrt(sort(eigen(t(Ztmp) %*% Ztmp)$values, decreasing = T))
  B1 <- qr.Q(B) %*% Diag(sqrt(eigvs))
  B0 <- rho[1] * B1 %*% Diag(sign(B1[1, ])) ## Fixed B0 and mu0 for each repeat.
  mu0 <- 0.4 * rnorm(p)
  Bm0 <- cbind(mu0, B0)
  set.seed(seed)
  
  epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
  H <- mvrnorm(n, mu = rep(0, q), cor.mat(q, 0.5))
  svdH <- svd(cov(H))
  H0 <- scale(H, scale = F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*% 
    svdH$v
  
  g1 <- 1:p
  B0[g1, ] <- B0[g1, ]/max(B0[g1, ]) * factor_term ## scale
  
  
  
  
  mu <- ai*exp( H0 %*% t(B0) + matrix(mu0, n, p, byrow = T) + epsi)
  
  X <- matrix(rpois(n * p, lambda = mu), n, p)
  
  return(list(X = X, B0 = B0, H0 = H0, mu0 = mu0))
}

MCCor <- function(H, H0){
  median(cancor(H, H0)$cor)
}
trace_statistic_fun <- function(H, H0){
  
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% qr.solve(t(H) %*% H) %*% t(H) %*% H0
  
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
  
}

factorm <- function(X, q=NULL){
  
  signrevise <- GFM:::signrevise
  if ((!is.null(q)) && (q < 1)) 
    stop("q must be NULL or other positive integer!")
  if (!is.matrix(X)) 
    stop("X must be a matrix.")
  mu <- colMeans(X)
  X <- scale(X, scale = FALSE)
  
  n <- nrow(X)
  p <- ncol(X)
  if (p > n) {
    svdX <- eigen(X %*% t(X))
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hatF <- as.matrix(svdX$vector[, 1:q] * sqrt(n))
    B2 <- n^(-1) * t(X) %*% hatF
    sB <- sign(B2[1, ])
    hB <- B2 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(k) hatF[, k] * sign(B2[1, 
                                                      ])[k])
  }
  else {
    svdX <- eigen(t(X) %*% X)
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hB1 <- as.matrix(svdX$vector[, 1:q])
    hH1 <- n^(-1) * X %*% hB1
    svdH <- svd(hH1)
    hH2 <- signrevise(svdH$u * sqrt(n), hH1)
    if (q == 1) {
      hB1 <- hB1 %*% svdH$d[1:q] * sqrt(n)
    }
    else {
      hB1 <- hB1 %*% diag(svdH$d[1:q]) * sqrt(n)
    }
    sB <- sign(hB1[1, ])
    hB <- hB1 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(j) hH2[, j] * sB[j])
  }
  sigma2vec <- colMeans((X - hH %*% t(hB))^2)
  res <- list()
  res$hH <- hH
  res$hB <- hB
  res$mu <- mu
  res$q <- q
  res$sigma2vec <- sigma2vec
  res$propvar <- sum(evalues[1:q])/sum(evalues)
  res$egvalues <- evalues
  attr(res, "class") <- "fac"
  return(res)
}


PLNPCA_run <- function(X_count, covariates, q,  Offset=rep(ai, nrow(X_count)), workers=NULL){
  require(PLNmodels)
  if(!is.null(workers)){
    future::plan("multisession", workers = workers)
  }
  d <- ncol(covariates)
  dat_plnpca <- prepare_data(X_count, covariates, Offset)
  formu <- paste0("Abundance ~ 1 + offset(log(Offset))+",paste(paste0("V",1:d), collapse = '+'))
  myPCA <- PLNPCA(as.formula(formu), data = dat_plnpca, ranks = q)
  myPCA1 <- getBestModel(myPCA)
  res_plnpca <- list(PCs= myPCA1$scores, bbeta= myPCA1$model_par$Theta, 
                     loadings=myPCA1$model_par$C)
  return(res_plnpca)
}


PLNPCA_runno <- function(X_count, covariates, q,  Offset=rep(ai, nrow(X_count)), workers=NULL){
  require(PLNmodels)
  if(!is.null(workers)){
    future::plan("multisession", workers = workers)
  }
  dat_plnpca <- prepare_data(X_count, covariates, Offset)
  formu <- paste0("Abundance ~ 1 + offset(log(Offset))")
  myPCA <- PLNPCA(as.formula(formu), data = dat_plnpca, ranks = q)
  myPCA1 <- getBestModel(myPCA)
  res_plnpca <- list(PCs= myPCA1$scores, bbeta= myPCA1$model_par$Theta, 
                     loadings=myPCA1$model_par$C)
  
  return(res_plnpca)
}


mrrr_run <- function(Y, X, ai,rank0, q, family=list(poisson()),
                     familygroup=rep(1,ncol(Y)), epsilon = 1e-4, sv.tol = 1e-2, maxIter = 30, trace=TRUE){
  # epsilon = 1e-4; sv.tol = 1e-2; maxIter = 30; trace=TRUE
  n <- nrow(Y); p <- ncol(Y)
   X <- cbind(X, diag(n))
  
  svdX0d1 <- svd(X)$d[1]
  init1 = list(kappaC0 = svdX0d1 * 5)
  offset = matrix(log(ai), n, p)
  control = list(epsilon = epsilon, sv.tol = sv.tol, maxit = maxIter,
                 trace = trace, gammaC0 = 1.1, plot.cv = TRUE,
                 conv.obj = TRUE)
  fit.mrrr <- mrrr(Y=Y, X=X,  offset= offset,family = family, familygroup = familygroup,
                   penstr = list(penaltySVD = "rankCon", lambdaSVD = 0.1),
                   control = control, init = init1, maxrank = rank0+q)
  
  return(fit.mrrr)
}



