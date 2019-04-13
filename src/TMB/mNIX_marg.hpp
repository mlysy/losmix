/// @file mNIX_marg.hpp

#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type mNIX_marg(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(Xtr);
  DATA_MATRIX(y);
  DATA_IVECTOR(iStart);
  DATA_IVECTOR(nObs);
  // parameters
  PARAMETER_VECTOR(lambda);
  PARAMETER_VECTOR(logC_Omega);
  PARAMETER(log_nu);
  PARAMETER(log_tau);
  // internal variables
  int nSub = nObs.size(); // number of subjects
  int p = Xtr.rows(); // number of covariates
  matrix<Type> Omega(p,p);
  utils<Type>::lchol2var(Omega, logC_Omega.matrix());
  Type nu = exp(log_nu);
  Type tau = exp(log_tau);
  mNIX<Type> mnix(p);
  // negative log-marginal
  Type nll = nSub * mnix.zeta(Omega, nu, tau);
  mnix.set_prior(lambda, Omega, nu, tau); // conjugate prior hyperparameters
  for(int ii=0; ii<nSub; ii++) {
    mnix.set_suff(y.block(iStart[ii],0,nObs[ii],1),
		  Xtr.block(0,iStart[ii],p,nObs[ii]));
    mnix.calc_post(); // conjugate posterior hyperparameters
    nll -= mnix.zeta();
  }
  // change-of-variables correction
  nll -= utils<Type>::lchol_prior(logC_Omega);
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
