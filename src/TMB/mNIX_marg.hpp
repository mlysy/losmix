/// @file mNIX_marg.hpp

#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type mNIX_marg(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  DATA_IVECTOR(iStart);
  DATA_IVECTOR(nObs);
  // parameters
  PARAMETER_VECTOR(lambda);
  PARAMETER_VECTOR(logC_Omega);
  PARAMETER(log_nu);
  PARAMETER(log_tau);
  // internal variables
  int nSub = nObs.size(); // number of subjects
  int p = X.cols(); // number of covariates
  matrix<Type> Omega = utils<Type>::lchol2var(logC_Omega);
  Type nu = exp(log_nu);
  Type tau = exp(log_tau);
  mNIX<Type> Phi(p);
  // negative log-marginal
  Type nll = nSub * mNIX<Type>::zeta(Omega, nu, tau);
  Phi.set_prior(lambda, Omega, nu, tau); // conjugate prior hyperparameters
  for(int ii=0; ii<nSub; ii++) {
    Phi.set_suff(y.segment(iStart[ii],nObs[ii]).matrix(),
		 X.block(iStart[ii],0,nObs[ii],p));
    Phi.calc_post(); // conjugate posterior hyperparameters
    nll -= Phi.zeta();
  }
  nll *= .5;
  // change-of-variables correction
  nll -= utils<Type>::lchol_prior(logC_Omega);
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
