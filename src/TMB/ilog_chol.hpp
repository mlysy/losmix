/// @file ilog_chol.hpp

#include "losmix/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type ilog_chol(objective_function<Type>* obj) {
  using namespace losmix;
  // data
  DATA_MATRIX(logC);
  PARAMETER(theta); // dummy parameter
  SIMULATE {
    // size of output
    int nOut = logC.cols();
    int p = (-1 + sqrt(1 + 8*logC.rows()))/2;
    matrix<Type> V(p, nOut*p);
    for(int ii=0; ii<nOut; ii++) {
      utils<Type>::lchol2var(V.block(0,ii*p,p,p), logC.col(ii));
    }
    REPORT(V);
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
