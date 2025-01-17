/// @file mNIX_NLL.hpp
///
/// Negative loglikelihood for the LOSMIX model.
///
/// @note Depends on `TMB.hpp` which is *not* header-guarded, so don't include it here.

#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type mNIX_NLL(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(Xtr);
  DATA_MATRIX(y);
  // hyperparameters
  PARAMETER_VECTOR(lambda);
  PARAMETER_VECTOR(logC_Omega); // Omega on the log-cholesky scale
  PARAMETER(log_nu);
  PARAMETER(log_tau);
  // parameters of mNIX posterior
  int N = Xtr.cols();
  int p = Xtr.rows();
  matrix<Type> Omega(p,p);
  utils<Type>::lchol2var(Omega, logC_Omega.matrix());
  Type nu = exp(log_nu);
  Type tau = exp(log_tau);
  matrix<Type> lambda_hat(p,1);
  matrix<Type> Omega_hat(p,p);
  Type nu_hat;
  Type tau_hat;
  // conversion class
  mNIX<Type> mnix(p);
  mnix.set_suff(y, Xtr);
  mnix.set_prior(lambda.matrix(), Omega, nu, tau);
  mnix.get_post(lambda_hat, Omega_hat, nu_hat, tau_hat);
  SIMULATE {
    // sufficient statistic outputs
    Type yy;
    matrix<Type> Xy(p,1);
    matrix<Type> XX(p,p);
    mnix.get_suff(yy, Xy, XX);
    REPORT(yy);
    REPORT(Xy);
    REPORT(XX);
    REPORT(N);
    // hyperparameters of conjugate posterior
    REPORT(lambda_hat);
    Omega_hat = XX + Omega;
    REPORT(Omega_hat);
    REPORT(nu_hat);
    REPORT(tau_hat);
  }
  // normalizing constant
  mnix.calc_post();
  return - mnix.log_marg();
  // return .5 * (mnix.zeta(Omega, nu, tau) -
  // 	       mnix.zeta(Omega_hat, nu_hat, tau_hat));
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
