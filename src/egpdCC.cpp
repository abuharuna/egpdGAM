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
double egpdcd0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate)
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

double cens, lpvi, lpxi, lpk;
double ee1, ee0;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

cens = yvec[j];
lpk = lpkvec[j];
lpvi = lpvivec[j];
lpxi = lpxivec[j];


ee0 = exp(-lpxi);
ee1 = cens*(exp(lpxi-lpvi) + exp(2*lpxi-lpvi));

if (ee1 <= -1.0) {
  nllh = 1e20;
  break;
} else {
  nllh += -log(pow((1 - pow((1 + cens*(exp(lpxi-lpvi) + exp(2*lpxi-lpvi))), (-ee0))), (exp(lpk)))) ;
  
}




}

return(nllh);

}
// //' @rdname egpd0
// [[Rcpp::export]]
arma::mat  egpdcd12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
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
  
  double cens, lpk, lpvi, lpxi;
  double e3, e4,  e6, e7,  e8, e9, e10;  
  double e11, e12, e14, e13,  e16, e17, e18,   e19;  
  double e22,  e25,e28;  
  double e30, e31,  e32, e35,  e37, e38; 
  double e41,  e46;
  
  for (int j=0; j < nobs; j++) {
    
    cens = yvec[j];
    lpk = lpkvec[j];
    lpvi = lpvivec[j];
    lpxi = lpxivec[j];
    // 1 = k, 2 = log(sigma), 3 = log(xi)
    //order: 1, 2, 3, 11, 12, 13, 22, 23, 33
    
    e3 = exp(2 * lpxi - lpvi);
    e4 = exp(lpxi - lpvi);
    e6 = exp(-lpxi);
    e7 = e3 + e4;
    e8 = cens * e7;
    e9 = 1 + e8;
    e10 = pow(e9, e6);
    e11 = 1 + e6;
    e12 = pow(e9, e11);
    e13 = 1/e10;
    e14 = 1 - e13;
    e16 = 2 * e3 + e4;
    e17 = exp(lpk);
    e18 = log1p(e8);
    e19 = e14 * e12;
    e22 = 1/e12;
    e25 = cens * pow(e9, (e6 - 1)) * e16 - e10 * e18;
    e28 = cens * e6;
    e30 = e18/e10 - cens * e16/e12;
    e31 = pow(e9, (2 * e6));
    e32 = e25 * e6;
    e35 = e28 * e7 * e17/e19;
    e37 = pow(e19, 2);
    e38 = e14 * e31;
    e41 = 2 * e11;
    e46 = cens * e10 * e11 * e16 - e12 * e6 * e18;
    
    out(j, 0) = -(e17 * log(e14));
    out(j, 1) = e35;
    out(j, 2) = e6 * e17 * e30/e14;
    out(j, 3) = -(e17 * log(e14));
    out(j, 4) = e35;
    out(j, 5) = e6 * e17 * e30/e14;
    out(j, 6) = cens * (cens * (e14 * e10 * e11 + e6) * e7/e37 - 1/e19) * e6 * e7 * e17;
    out(j, 7) = cens * (((e30/e19 + e18/e12) * e6 - e22) * e7 - e16 * (cens * pow(e9, (e6 - e41)) * e11 * e7 - e22)) * 
      e6 * e17/e14;
    out(j, 8) = (cens * ((e46/pow(e9, e41) + e22 + e22) * e16 -  (4 * e3 + e4)/e12) - ((e32/e31 + e13) * e18 + e32 * e30/e38)) * 
      e6 * e17/e14;
    
    
  }
  
  return out;
}
