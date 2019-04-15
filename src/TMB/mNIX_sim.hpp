/// @file mNIX_sim.hpp

#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type mNIX_sim(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(Xtr);
  DATA_MATRIX(y);
  DATA_IVECTOR(iStart);
  DATA_IVECTOR(nObs);
  // parameters (input as data because can't resize anyways)
  DATA_MATRIX(lambda);
  DATA_MATRIX(Omega);
  DATA_VECTOR(nu);
  DATA_VECTOR(tau);
  DATA_INTEGER(nOut);
  DATA_INTEGER(doPost);
  PARAMETER(theta); // dummy parameter
  SIMULATE {
    int p = Omega.rows();
    // output dimensions
    matrix<Type> beta(p,nOut);
    vector<Type> sigma(nOut);
    // internal variables
    mNIX<Type> mnix(p);
    Eigen::LLT<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > llt(p);
    // determine whether each input is single or vectorized
    bool vData = nObs.size() > 1;
    bool vLambda = lambda.cols() > 1;
    bool vOmega = Omega.cols() > p;
    bool vNu = nu.size() > 1;
    bool vTau = tau.size() > 1;
    bool vPars = vLambda || vOmega || vNu || vTau;
    if(doPost == 0) {
      // arbitrary mNIX distribution
      for(int ii=0; ii<nOut; ii++) {
	if(vOmega || ii == 0) {
	  llt.compute(Omega.block(0,p*ii,p,p));
	}
	mnix.simulate(beta.col(ii), sigma(ii),
		     lambda.col(vLambda*ii), llt, nu(vNu*ii), tau(vTau*ii));
      }
    } else {
      // posterior mNIX simulation
      for(int ii=0; ii<nOut; ii++) {
	if(vData || ii == 0) {
	  mnix.set_suff(y.block(iStart[ii],0,nObs[ii],1),
			Xtr.block(0,iStart[ii],p,nObs[ii]));
	}
	if(vPars || ii == 0) {
	  mnix.set_prior(lambda.col(vLambda*ii),Omega.block(0,vOmega*p*ii,p,p),
			nu(vNu*ii), tau(vTau*ii));
	}
	if(vPars || vData || ii == 0) {
	  mnix.calc_post();
	}
	mnix.simulate(beta.col(ii), sigma(ii));
      }
    }
    REPORT(beta);
    REPORT(sigma);
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
