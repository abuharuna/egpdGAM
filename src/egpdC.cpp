// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps  = -11.51293;

// //' extended Generalized Pareto distribution (eGPD) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each eGPD parameter
// //' @param X1 a design matrix for the eGPD log scale parameter
// //' @param X2 a design matrix for the eGPD log shape parameter
// //' @param X3 a design matrix for the eGPD kappa parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return egpdd0 a scalar, the negative log-likelihood
// //' @return egpdd12 a matrix, first then second derivatives w.r.t. eGPD parameters
// //' @return egpdd34 a matrix, third then fourth derivatives w.r.t. eGPD parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double egpdd0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
arma::vec lpkvec = X1 * Rcpp::as<arma::vec>(pars[0]);    
arma::vec lpvivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec lpxivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

if (dcate == 1) {
    lpkvec = lpkvec.elem(dupid);
    lpvivec = lpvivec.elem(dupid);
    lpxivec = lpxivec.elem(dupid);
}

double y, lpvi, lpxi, lpk;
double ee1, ee0;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpk = lpkvec[j];
lpvi = lpvivec[j];
lpxi = lpxivec[j];


ee0 = exp(-lpxi);
ee1 =y*(exp(lpxi-lpvi) + exp(2*lpxi-lpvi));

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {
  nllh += lpvi - log(1.0 + exp(lpxi)) + (1.0 + ee0) * log(1.0 + ee1) - (exp(lpk) -  1.0) * log(1 - pow((1.0 + ee1), (-ee0))) - lpk;
}




}

return(nllh);

}
// //' @rdname egpd0
// [[Rcpp::export]]
arma::mat  egpdd12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
{
  
  arma::vec lpkvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lpvivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lpxivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9);
  
  if (dcate == 1) {
    lpkvec = lpkvec.elem(dupid);
    lpvivec = lpvivec.elem(dupid);
    lpxivec = lpxivec.elem(dupid);
  }
  
  double y, lpk, lpvi, lpxi;
  double e3, e4, e5, e6,  e8, e9, e10;  
  double e11, e12, e13,  e15, e16, e17, e18,   e19;  
  double e20, e21, e22, e23, e24, e25, e26,  e28, e29;  
  double e30,  e32, e33, e36, e37, e38, e39; 
  double e40,  e42, e43, e44, e45,    e48, e49; 
  double e53; 
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    lpk = lpkvec[j];
    lpvi = lpvivec[j];
    lpxi = lpxivec[j];
    // 1 = k, 2 = log(sigma), 3 = log(xi)
    //order: 1, 2, 3, 11, 12, 13, 22, 23, 33
    e3 = exp(2 * lpxi - lpvi);
    e4 = exp(lpxi - lpvi);
    e5 = e3 + e4;
    e6 = y * e5;
    e8 = exp(-lpxi);
    e9 = 1 + e6;
    e10 = pow(e9, e8);
    e11 = 1 + e8;
    e12 = 1/e10;
    e13 = 1 - e12;
    e15 = 2 * e3 + e4;
    e16 = pow(e9, e11);
    e17 = exp(lpk);
    e18 = log1p(e6);
    e19 = e17 - 1;
    e20 = e13 * e10;
    e21 = exp(lpxi);
    e22 = e18/e10;
    e23 = y * e15;
    e24 = pow(e9, (e8 - 1));
    e25 = 1/e16;
    e26 = e22 - e23/e16;
    e28 = (e19/e20 - 1) * e8 - 1;
    e29 = e13 * e16;
    e30 = pow(e9, (2 * e8));
    e32 = 1 + e21;
    e33 = y * e11;
    e36 = y * e24 * e15 - e10 * e18;
    e37 = pow(e20, 2);
    e38 = e13 * e30;
    e39 = 1/e20;
    e40 = 2 * e11;
    e42 = 4 * e3 + e4;
    e43 = e8 * e17;
    e44 = e17 * log(e13);
    e45 = e21/e32;
    e48 = y * e28 * e5/e9;
    e49 = e6/e9;
    e53 = y * e8 * e5 * e17/e29;
    
    out(j, 0) = lpk = -(1 + e44);
    out(j, 1) = 1 + e48;
    out(j, 2) = (e19 * e26/e13 - e18) * e8 + e33 * e15/e9 - e45;
    out(j, 3) = -e44;
    out(j, 4) = e53;
    out(j, 5) = e43 * e26/e13;
    out(j, 6) = y * ((e19 * (y * (e13 * e24 + 1/e9) * e8 * e5/e37 -  e39) + 1) * e8 + 1 + e48) * e5/e9;
    out(j, 7) = y * (((((e26/e29 + e18/e16) * e8 - e25) * e5 -  e15 * (e33 * pow(e9, (e8 - e40)) * e5 - e25)) * e19/e13 + e5/e9) * e8 + e11 * e15 * (e49 - 1)/e9);
    out(j, 8) = (e19 * (y * (((e33 * e10 * e15 - e16 * e8 * e18)/pow(e9, e40) + e25 + e25) * e15 - e42/e16) - ((e26/e38 + e18/e30) *  e8 * e36 + e22))/e13 + e18 - e23/e9) *
      e8 + y * (e11 *  (e42 - y * pow(e15, 2)/e9) - e15 * e8)/e9 - (1 - e45) * e21/e32;
    
    
    
  }
  
  return out;
}
