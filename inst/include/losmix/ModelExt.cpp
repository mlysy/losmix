
#include <TMB.hpp> // _must_ put this first
// include losmix library
#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace losmix;
  // --- INPUT VARIABLES ---
  // data input: output of the call losmix::format_data
  DATA_MATRIX(y); // response vector y (as a matrix)
  DATA_MATRIX(Xtr); // _transpose_ of covariate matrix X
  DATA_IVECTOR(iStart); // start of each subject
  DATA_IVECTOR(nObs); // number of observations per subject
  // parameter input (order is arbitrary)
  PARAMETER(gamma);
  PARAMETER_VECTOR(lambda);
  PARAMETER_VECTOR(logC_Omega); // cholesky factor of Omega
  PARAMETER(log_nu);
  PARAMETER(log_tau);
  
  // --- INTERMEDIATE VARIABLES --- 
  int nSub = nObs.size(); // number of subjects
  int p = lambda.size(); // number of regression coefficients
  // hyperparameter conversions of Omega, nu, tau
  matrix<Type> Omega(p,p);
  // lchol2var expects one-column matrix as second input
  utils<Type>::lchol2var(Omega, logC_Omega.matrix()); 
  Type nu = exp(log_nu);
  Type tau = exp(log_tau);
  matrix<Type> XtrGamma(p, Xtr.cols()); // hyperparameter-dependent covariate matrix
  XtrGamma.row(0).setOnes();
  XtrGamma.row(1) = (-gamma * Xtr.row(0)).array().exp();
  Type nuX; // covariate-dependent degrees-of-freedom
  mNIX<Type> mnix(p); // mNIX distribution object

  // --- OUTPUT VARIABLES ---
  Type mll = 0; // marginal log-likelihood
  // random effects variables
  matrix<Type> beta(p, nSub); // regression coefficients
  vector<Type> sigma(nSub); // error standard deviations

  // --- CALCULATIONS ---
  for(int ii=0; ii<nSub; ii++) {
    nuX = nu * Xtr(1,iStart(ii)); // covariate-dependent degrees of freedom
    // posterior calculations
    mnix.set_prior(lambda, Omega, nuX, tau);
    mnix.set_suff(y.block(iStart[ii],0,nObs[ii],1),
		  XtrGamma.block(0,iStart[ii],p,nObs[ii]));
    mnix.calc_post();
    mll += mnix.log_marg(); // increment marginal likelihood
    SIMULATE {
      // simulate random effect for each subject
      // enclose in SIMULATE block to avoid AD derivative calculations,
      // i.e., faster
      mnix.simulate(beta.col(ii), sigma(ii));
    }
  }
  // hyperparameter prior
  // the default prior on unconstrained scalars is "uniform(-Inf,Inf)"
  // the default prior on the log-Cholesky decomposition of a variance matrix
  // is the 'lchol_prior' below.
  mll += utils<Type>::lchol_prior(logC_Omega);
  // return random effects
  SIMULATE {
    REPORT(beta);
    REPORT(sigma);
  }
  return -mll; // TMB expects a _negative_ log marginal posterior
}
