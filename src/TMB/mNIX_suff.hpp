/// @file mNIX_suff.hpp

#include "losmix/utils.hpp"
#include "losmix/mNIX.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type mNIX_suff(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  DATA_IVECTOR(iStart);
  DATA_IVECTOR(nObs);
  PARAMETER(theta); // dummy parameter
  SIMULATE {
    int p = X.cols();
    int nSub = nObs.size();
    mNIX<Type> Phi(p);
    vector<Type> yy(nSub);
    matrix<Type> Xy(p, nSub);
    matrix<Type> XX(p, p*nSub);
    for(int ii=0; ii<nSub; ii++) {
      Phi.set_suff(y.segment(iStart[ii],nObs[ii]).matrix(),
		   X.block(iStart[ii],0,nObs[ii],p));
      Phi.get_suff(yy(ii), Xy.col(ii), XX.block(0,p*ii,p,p));
    }
    REPORT(yy);
    REPORT(Xy);
    REPORT(XX);
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
