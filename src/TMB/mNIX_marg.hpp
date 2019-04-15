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
  // output variables
  Type nll = 0; // negative log-marginal
  // random effects
  matrix<Type> beta(p, nSub); // regression coefficients
  vector<Type> sigma(nSub); // error standard deviations
  mnix.set_prior(lambda, Omega, nu, tau); // conjugate prior hyperparameters
  nll = nSub * mnix.zeta(Omega, nu, tau); // normalizing constant: prior half
  for(int ii=0; ii<nSub; ii++) {
    mnix.set_suff(y.block(iStart[ii],0,nObs[ii],1),
		  Xtr.block(0,iStart[ii],p,nObs[ii]));
    mnix.calc_post(); // conjugate posterior hyperparameters
    nll -= mnix.zeta(); // normalizing constant: posterior half
    SIMULATE {
      // random effects posterior draw
      mnix.simulate(beta.col(ii), sigma(ii));
    }
  }
  // change-of-variables correction
  nll -= utils<Type>::lchol_prior(logC_Omega);
  SIMULATE {
    REPORT(beta);
    REPORT(sigma);
  }
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
