// This script implement covariate-augumented Poisson factor model.
// Date: 2022-12-27

// Revised log:
// 2023-07-16: replace the for loop with matrix operation in updating variational parameters of Z of E-step!
// 2023-07-31: Regard H as an infinite parameter and estimate it in the algorithm.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 

// diag(W0^t* Cki * W0)
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  


List irlbaCpp(const mat& X, const int& q){
  Rcpp::Environment irlba("package:irlba");
  Rcpp::Function f = irlba["irlba"];
  
  return f(Named("A") = X, Named("nv") = q);
  
}  

double calELBO(const mat& X_count, const vec& a, const mat& Z, const mat& Mu_y, const mat& S_y,
               const vec& invLambda, const mat& B, const mat& bbeta, const mat& H){
  
  int i, n = X_count.n_rows, p = X_count.n_cols, q= B.n_cols;
  double pois_term1, pois_term2, logS, entropy=0.0, ELBO=0.0;
  pois_term1 = accu(X_count % Mu_y); 
  pois_term2 = -accu(repmat(a, 1, p) % exp(Mu_y + S_y/2));
  
  
  mat dX = (Mu_y - Z*bbeta.t() - H * B.t()) % repmat(sqrt(invLambda.t()), n, 1);
  mat LamB = B % repmat(sqrt(invLambda), 1, q); 
  double dimreduction_term1 = -0.5*(accu(dX % dX)+ // trace(LamB * S_bar * LamB.t())
                                    accu(S_y % repmat(invLambda.t(), n, 1)) - n* accu(log(invLambda)) );
  
  entropy = 0.5*accu(log(S_y));
  
  ELBO = pois_term1 + pois_term2 + dimreduction_term1 + entropy;
  return ELBO;
}

// update bbeta
// separate method
mat update_bbeta_sep(const mat& Z, const mat&Y, const vec& invLambda, const int& rank_use, const bool& fast_svd=true){ // can be speeded up !!!!
  
  // Perform singular value decomposition of X
  int n = Y.n_rows, d = Z.n_cols;
  mat C_ls = pinv(Z.t() * Z) * Z.t()*Y; // d*p
  rowvec sqrt_sigma_inv = sqrt(invLambda.t());
  //Rprintf("good2\n");
  mat ZC = Z * (C_ls % repmat(sqrt_sigma_inv, d, 1)); // n*p
  mat V;
  if(fast_svd){
    Rcpp::List svdX = irlbaCpp(ZC, rank_use);
    mat V1 = svdX["v"];
    V = V1;
    V1.reset();
  }else{
    mat U, V1;
    vec s;
    svd(U, s, V1, ZC);  // how to speed up using approximate SVD
    V = V1;
    U.reset();
    V1.reset();
  }
  ZC.reset();
  mat C_rr = (C_ls % repmat(sqrt_sigma_inv, d, 1))* (V.cols(0, rank_use-1)) * (trans(V.cols(0, rank_use-1)) % repmat(1/sqrt_sigma_inv, rank_use, 1)) ; // d*p
  //Rprintf("good2\n");
  
  return C_rr.t();
  
}
// joint method
mat update_bbeta_joint(const mat& Z, const mat&Y, const mat& B, const mat& S_y, const int& rank_use, const bool& fast_svd=true){ // can be speeded up !!!!
  
  // Perform singular value decomposition of X
  int n = Y.n_rows, d = Z.n_cols;
  mat C_ls = pinv(Z.t() * Z) * Z.t()* Y; // d*p
  mat dX = Y-Z*C_ls;
  vec lambdaVec =  trans(mean(dX % dX + S_y)); // tSigma_eps
  dX.reset();
  
  
  // mat C_ls = pinv(Z.t() * Z) * Z.t()*Y; // d*p
  rowvec sqrt_sigma_inv = sqrt(1/lambdaVec.t());
  //Rprintf("good2\n");
  mat ZC = Z * (C_ls % repmat(sqrt_sigma_inv, d, 1)); // n*p
  //Rprintf("good1\n");
  mat V;
  if(fast_svd){
    Rcpp::List svdX = irlbaCpp(ZC, rank_use);
    mat V1 = svdX["v"];
    V = V1;
    V1.reset();
  }else{
    mat U, V1;
    vec s;
    svd(U, s, V1, ZC);  // how to speed up using approximate SVD
    V = V1;
    U.reset();
    V1.reset();
  }
  ZC.reset();
  mat C_rr = (C_ls % repmat(sqrt_sigma_inv, d, 1))* (V.cols(0, rank_use-1)) * (trans(V.cols(0, rank_use-1)) % repmat(1/sqrt_sigma_inv, rank_use, 1)) ; // d*p
  
  return C_rr.t();
  
}

void VB_Estep(const mat& X_count, const vec& a, const mat& Z, mat& Mu_y,  mat& S_y,
              const vec& invLambda, const mat& B, const mat& bbeta, mat& H){
  
  int i, j, n = X_count.n_rows, p = X_count.n_cols, q= B.n_cols;
  //  ## VB E-step
  // update posterior variance of y: S_y
  // update posterior mean of y: Mu_y
  // double elbo0 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H, S_h, Sigma_h);
  // Laplace approximation plus Taylor approximation.
  // Method 2: Laplace approximation + Tylor approximation
  // double to_ij =0.0, tz_ij;
  // for(i = 0; i<n ; ++i){
  //   for(j = 0; j<p; ++j){
  //     to_ij = a(i) * exp(Mu_y(i,j));
  //     tz_ij = tau(j)+ as_scalar(Z.row(i) * bbeta.row(j).t()+ H.row(i) * B.row(j).t());
  //     Mu_y(i,j) = (X_count(i,j)- to_ij*(1-Mu_y(i,j)) + invLambda(j)* tz_ij) / (invLambda(j) + to_ij);
  //     
  //     S_y(i,j) = 1.0/(a(i)*exp(Mu_y(i,j)) + invLambda(j));
  //   }
  // }
  // Matrix operation for this method:
  mat tmp_mat = Z * bbeta.t() + H * B.t();
  Mu_y = (X_count - repmat(a, 1, p) % exp(Mu_y) % (1- Mu_y) + repmat(invLambda.t(), n,1) % tmp_mat) / 
    (repmat(a, 1, p) % exp(Mu_y) + repmat(invLambda.t(), n,1) );
  S_y = 1.0 / (repmat(a, 1, p) % exp(Mu_y) + repmat(invLambda.t(), n,1) );
  //}
  // double elbo1 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H, S_h, Sigma_h);
  // Rprintf("VB dY= %4f\n", elbo1 - elbo0);
  
  // ## update posterior variance of h_i, S_i
  // ## update posterior mean of h_i, H
  // mat Si_inv;
  // Si_inv = B.t() *  (B % repmat(invLambda, 1, q)) + inv_sympd(Sigma_h); // diagmat(invLambda) * B
  // S_h = inv_sympd(Si_inv);
  // H = (Mu_y - Z*bbeta.t()) * (B % repmat(invLambda, 1, q)) * S_h;
  
  // double elbo2 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H, S_h, Sigma_h);
  // Rprintf("VB dH= %4f\n", elbo2 - elbo1);
}




// Augumented Poisson factor model
// [[Rcpp::export]]
Rcpp::List Reduced_Rank_COAP(const arma::mat& X_count, const arma::vec& a, const arma::mat& Z, const int& rank_use,
                              const arma::mat& Mu_y_int, 
                              const arma::mat& S_y_int,
                              const arma::vec& invLambda_int, const arma::mat& B_int, const arma::mat& bbeta_int,
                              const arma::mat& H_int, 
                              const double& epsELBO, const int& maxIter, const bool& verbose, 
                              const bool& sep_opt_beta, const bool& fast_svd=true){
  
  int i, n = X_count.n_rows, p = X_count.n_cols, q= B_int.n_cols;
  
  //bool sep_opt_beta = true;
  // Initialize
  mat Mu_y(Mu_y_int), S_y(S_y_int), B(B_int), bbeta(bbeta_int), H(H_int);
  vec invLambda(invLambda_int);
  
  vec ELBO_vec(maxIter), bsb, Lambda;;
  ELBO_vec(0) = INT_MIN;
  mat S_bar, dX, tY;
  int iter;
  
  
  for(iter = 1; iter < maxIter; ++iter){
    
    
    // Rprintf("E step starting!\n");
    // VB E-step
    VB_Estep(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H);
    
    //VB M-step
    //double elbo1 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H);
    //update B
    B = trans(Mu_y - Z * bbeta.t()) * H * inv_sympd(H.t() * H);
    // double elbo2 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H);
    // Rprintf("dB= %4f \n", elbo2 - elbo1);
    
    // update H
    // H = (Mu_y - Z * bbeta.t()) * B * inv_sympd(B.t() * B);
    H = (Mu_y - Z * bbeta.t()) *(B % repmat(invLambda, 1, q)) * inv(B.t()*  (B % repmat(invLambda, 1, q)));
    // double elbo3 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H);
    // Rprintf("dH= %4f \n", elbo3 - elbo2);
    // update bbeta
    tY = Mu_y - H * B.t();
    if(sep_opt_beta){
      bbeta = update_bbeta_sep(Z, tY, invLambda, rank_use, fast_svd); // Fixed-rank regression.
    }else{
      bbeta = update_bbeta_joint( Z, tY,  B, S_y, rank_use, fast_svd);
    }
    //bbeta = trans(Mu_y - H * B.t()) * Z * inv_sympd(Z.t()*Z);
    
    // double elbo4 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H);
    // Rprintf("dbbeta= %4f \n", elbo4 - elbo3);
    
    // update Lambda
    dX = Mu_y -  Z * bbeta.t()- H * B.t();
    Lambda =  trans(mean(dX % dX + S_y));
    Lambda.elem( find(Lambda < 1e-4) ).fill(1e-4); // increase the numerical stability!
    invLambda = 1.0 / Lambda;
    // double elbo5 = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H);
    // Rprintf("dLambda= %4f\n", elbo5 - elbo4);
    
    
    ELBO_vec(iter) = calELBO(X_count, a, Z, Mu_y, S_y, invLambda, B, bbeta, H);
    
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n", 
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
  }
  
  // Add identifiability condition for H and B;
  mat alpha = pinv(Z.t() * Z) * Z.t()* H; // d*q
  mat H_new = H - Z * alpha;
  vec mu = bbeta.col(0) + B * alpha.row(0).t();
  // mat U1, V1;
  // vec s1;
  // svd(U1, s1, V1, H_new*B.t()); // speed up using approxPCA.
  // mat V2 = V1.cols(0, q-1), U2 = U1.cols(0, q-1);
  Rcpp::List svdX = irlbaCpp(H_new*B.t(), q);
  mat V2 = svdX["v"];
  mat U2 = svdX["u"];
  vec s1 = svdX["d"];
  rowvec signB1 = sign(V2.row(0));
  H_new = sqrt(n) * U2 * diagmat(signB1.t());
  mat B_new = V2 * diagmat(s1.subvec(0, q-1) % signB1.t()) / sqrt(n);
  // then re-estimate beta to ensure beta0 is rank-fixed
  //bbeta = trans(Mu_y - H * B.t()) * Z * inv_sympd(Z.t()*Z);
  bbeta = update_bbeta_sep(Z, Mu_y - H_new * B_new.t(), invLambda, rank_use); // Fixed-rank regression.
  mat bbeta_new = bbeta;
  bbeta_new.col(0) = mu;
  
  // output return value
  List resList = List::create(
    // Rcpp::Named("H") = H,
    // Rcpp::Named("B") = B,
    // Rcpp::Named("bbeta") = bbeta,
    Rcpp::Named("H") = H_new,
    Rcpp::Named("B") = B_new,
    Rcpp::Named("bbeta") = bbeta_new,
    Rcpp::Named("invLambda") = invLambda,
    Rcpp::Named("Mu_y") = Mu_y,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);
  
}