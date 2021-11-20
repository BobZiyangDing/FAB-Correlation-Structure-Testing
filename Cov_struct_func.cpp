#include <RcppArmadillo.h>
#include <cstdlib>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat getFABnUMPUpvals(arma::mat& Cor_F_r, 
                           arma::vec& F_r,
                           int n,
                           int p) {
  
  arma::mat one = ones(p, 1);
  arma::mat result = arma::mat(2, p);
  
  for(int j=0; j<p; j++){
    mat G_j = null(Cor_F_r.row(j));
    
    mat tG_j = G_j.t();
    mat tone = one.t();
    
    vec indirect_info = tG_j *F_r;
    vec eta_MLE = inv(tone*G_j*tG_j*one)*tone*G_j*indirect_info;
    vec mu_MLE = one * eta_MLE;
    double psi2_MLE = var(indirect_info) - 1/(n-3);
    
    mat temp = G_j*inv(tG_j *Cor_F_r*G_j/(n-3))*tG_j;
    mat V = inv(eye(p,p)/psi2_MLE + temp);
    vec m = V*((eye(p,p)/psi2_MLE)*mu_MLE+temp*F_r);
    
    double tilde_m = m(j);
    double tilde_v = V(j,j);
    double umpu_p = 1 - abs(normcdf(F_r(j)*sqrt(n-3)) - normcdf(-F_r(j)*sqrt(n-3)));
    double fab_p = 1 - abs(normcdf(F_r(j)*sqrt(n-3)+2*tilde_m*sqrt(n-3)/tilde_v) - normcdf(-F_r(j)*sqrt(n-3)));
    
    result(0, j) = umpu_p;
    result(1, j) = fab_p;
    // result(2, j) = tilde_m;
    // result(3, j) = tilde_v;
    
  }
  return result;
}


// [[Rcpp::export]]
arma::mat getGeneralFABnUMPUpvals(arma::mat& Cor_F_r, 
                                  arma::vec& F_r,
                                  arma::mat& one,
                                  double lambda,
                                  int n,
                                  int p) {
  
  // arma::mat one = ones(p, 1);
  arma::mat result = arma::mat(2, p);
  
  for(int j=0; j<p; j++){
    mat G_j = null(Cor_F_r.row(j));
    
    mat tG_j = G_j.t();
    mat tone = one.t();
    
    vec indirect_info = tG_j *F_r;
    vec eta_MLE = inv(tone*G_j*tG_j*one + lambda * eye(size(one)[1] ,size(one)[1]) )*tone*G_j*indirect_info;
    vec mu_MLE = one * eta_MLE;
    double psi2_MLE = var(indirect_info) - 1/(n-3);
    
    mat temp = G_j*inv(tG_j *Cor_F_r*G_j/(n-3))*tG_j;
    mat V = inv(eye(p,p)/psi2_MLE + temp);
    vec m = V*((eye(p,p)/psi2_MLE)*mu_MLE+temp*F_r);
    
    double tilde_m = m(j);
    double tilde_v = V(j,j);
    double umpu_p = 1 - abs(normcdf(F_r(j)*sqrt(n-3)) - normcdf(-F_r(j)*sqrt(n-3)));
    double fab_p = 1 - abs(normcdf(F_r(j)*sqrt(n-3)+2*tilde_m*sqrt(n-3)/tilde_v) - normcdf(-F_r(j)*sqrt(n-3)));
    
    result(0, j) = umpu_p;
    result(1, j) = fab_p;
    // result(2, j) = tilde_m;
    // result(3, j) = tilde_v;
    
  }
  return result;
}