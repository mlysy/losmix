/// @file mNIX_post.hpp

#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type mNIX_post(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  DATA_IVECTOR(iStart);
  DATA_IVECTOR(nObs);
  // parameters (input as data because can't resize anyways)
  DATA_MATRIX(lambda);
  DATA_MATRIX(Omega);
  DATA_VECTOR(nu);
  DATA_VECTOR(tau);
  DATA_INTEGER(nPost);
  PARAMETER(theta); // dummy parameter
  SIMULATE {
    int p = X.cols();
    // output dimensions
    matrix<Type> lambda_hat(p,nPost);
    matrix<Type> Omega_hat(p,p*nPost);
    vector<Type> nu_hat(nPost);
    vector<Type> tau_hat(nPost);
    // determine whether each input is single or vectorized
    bool vData = nObs.size() > 1;
    bool vLambda = lambda.cols() > 1;
    bool vOmega = Omega.cols() > p;
    bool vNu = nu.size() > 1;
    bool vTau = tau.size() > 1;
    bool vPars = vLambda || vOmega || vNu || vTau;  
    mNIX<Type> Phi(p);
    for(int ii=0; ii<nPost; ii++) {
      if(vData || ii == 0) {
	Phi.set_suff(y.segment(iStart[ii],nObs[ii]).matrix(),
		     X.block(iStart[ii],0,nObs[ii],p));
      }
      if(vPars || ii == 0) {
	Phi.set_prior(lambda.col(vLambda*ii),Omega.block(0,vOmega*p*ii,p,p),
		      nu(vNu*ii), tau(vTau*ii));
      }
      Phi.get_post(lambda_hat.col(ii), Omega_hat.block(0,p*ii,p,p),
		   nu_hat(ii), tau_hat(ii));
    }
    REPORT(lambda_hat);
    REPORT(Omega_hat);
    REPORT(nu_hat);
    REPORT(tau_hat);
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
