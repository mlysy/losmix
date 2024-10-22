// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_losmix_TMBExports
#include <TMB.hpp>
#include "ilog_chol.hpp"
#include "mNIX_marg.hpp"
#include "mNIX_NLL.hpp"
#include "mNIX_post.hpp"
#include "mNIX_sim.hpp"
#include "mNIX_suff.hpp"
#include "ModelExt.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "ilog_chol") {
    return ilog_chol(this);
  } else if(model == "mNIX_marg") {
    return mNIX_marg(this);
  } else if(model == "mNIX_NLL") {
    return mNIX_NLL(this);
  } else if(model == "mNIX_post") {
    return mNIX_post(this);
  } else if(model == "mNIX_sim") {
    return mNIX_sim(this);
  } else if(model == "mNIX_suff") {
    return mNIX_suff(this);
  } else if(model == "ModelExt") {
    return ModelExt(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
