/// @file mNIX_NLL.hpp
///
/// Negative loglikelihood for the LOSMIX model.
///
/// @note Depends on `TMB.hpp` which is *not* header-guarded, so don't include it here.
///
/// Tests the following:
///
/// - `template <class Type>` functions with `void` return type, for multi-argument assignment.
/// - `LDLT` through the AD mechanism.
///
/// Conclusions:
///
/// - `TMB` uses `vector<CppAD>` for `vector<Type>`, but `Eigen::matrix` for `matrix<Type>`.  Therefore, matrix-vector operations typically don't work as we would expect.  Current approach is to explicitly convert vectors to one-column matrices via `.matrix()`, and use temporary matrices for intermediate calculations when needed.
/// - Can't figure out how to initialize a defined matrix/vector<Type> of unknown size...

#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type mNIX_NLL(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  // DATA_MATRIX(lambda_hat);
  // DATA_MATRIX(Omega_hat);
  // DATA_SCALAR(tau_hat);
  // DATA_SCALAR(nu_hat);
  // hyperparameters
  PARAMETER_VECTOR(lambda);
  PARAMETER_VECTOR(Omega_LC); // Omega on the log-cholesky scale
  PARAMETER(nu);
  PARAMETER(tau);
  matrix<Type> Omega = utils<Type>::lchol2var(Omega_LC);
  // parameters of mNIX posterior
  int N = X.rows();
  int p = X.cols();
  matrix<Type> lambda_hat(p,1);
  matrix<Type> Omega_hat(p,p);
  Type nu_hat;
  Type tau_hat;
  // conversion class
  mNIX<Type> Phi;
  Phi.set_suff(y.matrix(), X);
  Phi.set_prior(lambda.matrix(), Omega, nu, tau);
  Phi.get_post(lambda_hat, Omega_hat, nu_hat, tau_hat);
  SIMULATE {
    // sufficient statistic outputs
    Type yy;
    matrix<Type> Xy(p,1);
    matrix<Type> XX(p,p);
    Phi.get_suff(yy, Xy, XX);
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
  return .5 * (mNIX<Type>::zeta(Omega, nu, tau) - mNIX<Type>::zeta(Omega_hat, nu_hat, tau_hat));
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
