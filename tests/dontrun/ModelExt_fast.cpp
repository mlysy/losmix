
#include <TMB.hpp>
#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type ModelExt_marg(objective_function<Type>* obj) {
  using namespace losmix;
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

  // intermediate variables
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
  vector<Type> lm(nSub);
  Type lpi;
  mNIX<Type> mnix(p); // mNIX distribution object

  // marginal log-likelihood
  Type mll = 0.0; 
  for(int ii=0; ii<nSub; ii++) {
    nuX = nu * Xtr(1,iStart(ii)); // covariate-dependent degrees of freedom
    // posterior calculations
    mnix.set_prior(lambda, Omega, nuX, tau);
    mnix.set_suff(y.block(iStart[ii],0,nObs[ii],1),
		  XtrGamma.block(0,iStart[ii],p,nObs[ii]));
    mnix.calc_post();
    lm(ii) = mnix.log_marg();
    mll += lm(ii);
    // mll += mnix.log_marg();
  }
  // hyperparameter prior
  // the default prior on unconstrained scalars is "uniform(-Inf,Inf)"
  // the default prior on the log-Cholesky decomposition of a variance matrix
  // is the 'lchol_prior' below.
  lpi = utils<Type>::lchol_prior(logC_Omega);
  mll += lpi;
  // mll += utils<Type>::lchol_prior(logC_Omega);
  SIMULATE {
    matrix<Type> beta(p,1);
    Type sigma;
    mnix.simulate(beta, sigma);
    REPORT(beta);
    REPORT(sigma);
  }
  REPORT(lm);
  REPORT(lpi);
  return -mll; // TMB expects a _negative_ log marginal posterior
}

template<class Type>
Type ModelExt_sim(objective_function<Type>* obj) {
  using namespace losmix;
  // data input: output of the call losmix::format_data
  DATA_MATRIX(y); // response vector y (as a matrix)
  DATA_MATRIX(Xtr); // _transpose_ of covariate matrix X
  DATA_IVECTOR(iStart); // dummy variable = 0 (only one subject)
  DATA_IVECTOR(nObs); // nObs(0) = number of observations for subject

  // hyperparameter inputs (as data since they don't change here)
  DATA_VECTOR(gamma); // vector of length `nSim` (number of simulations)
  DATA_MATRIX(lambda); // matrix of size `p x nSim`
  DATA_MATRIX(logC_Omega); // matrix of size `p(p+1)/2 x nSim`
  DATA_VECTOR(log_nu); // vector of length `nSim`
  DATA_VECTOR(log_tau); // vector of length `nSim`
  PARAMETER(theta); // dummy parameter for model to compile

  // simulation block prevents derivative calculations, so faster
  SIMULATE {
    // intermediate variables
    int nSim = gamma.size(); // number of subjects
    int p = lambda.rows(); // number of regression coefficients
    // hyperparameter conversion of Omega, nu, and tau
    matrix<Type> Omega(p, p); 
    Type nuX; // parameter-dependent degrees-of-freedom
    Type tau;
    matrix<Type> XtrGamma(p, Xtr.cols()); // hyperparameter-dependent covariate matrix
    XtrGamma.row(0).setOnes();
    mNIX<Type> mnix(p); // mNIX distribution object
    // output variables
    matrix<Type> beta(p, nSim);
    vector<Type> sigma(nSim);

    // simulation
    mnix.set_suff(y, Xtr);
    for(int ii=0; ii<nSim; ii++) {
      // parameter dependent covariate matrix
      XtrGamma.row(1) = (-gamma(ii) * Xtr.row(0)).array().exp();
      nuX = exp(log_nu(ii)) * Xtr(1,0); // covariate-dependent degrees of freedom
      tau = exp(log_tau(ii));
      utils<Type>::lchol2var(Omega, logC_Omega.col(ii));
      // set up conditional mNIX distribution
      mnix.set_prior(lambda.col(ii), Omega, nuX, tau);
      mnix.set_suff(y, XtrGamma);
      mnix.calc_post();
      // random draw
      mnix.simulate(beta.col(ii), sigma(ii));
    }
    REPORT(beta);
    REPORT(sigma);
  }
  return 0; // dummy return value
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model_name);
  if(model_name == "ModelExt_marg") {
    return ModelExt_marg(this);
  } else if(model_name == "ModelExt_sim") {
    return ModelExt_sim(this);
  } else {
    error("Unknown model_name.");
  }
  return 0;
}
