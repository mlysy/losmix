/// @file utils.hpp
///
/// Utility functions.
///
/// @note Depends on `TMB.hpp` which is *not* header-guarded, so don't include it here.

#ifndef losmix_utils_hpp
#define losmix_utils_hpp 1

namespace losmix {

  template <class Type>
  struct utils {
  public:
    /// Typedef equivalent to `MatrixXd<Type>`.
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_t;
    /// Typedef equivalent to `Ref <MatrixXd<Type> >`
    typedef Eigen::Ref <MatrixXd_t> RefMatrix_t;
    /// Typedef equivalent to `const Ref <const MatrixXd<Type> >`
    typedef const Eigen::Ref <const MatrixXd_t> cRefMatrix_t;
    
    /// Dot product for one-column matrices.
    ///
    /// @param[in] x First column matrix.
    /// @param[in] y Second column matrix.
    /// @return The scalar value of `x.transpose() * y`.
    static Type dot_product(cRefMatrix_t& x, cRefMatrix_t& y) {
      return (x.array() * y.array()).sum();
    }

    /// Convert log-Cholesky decomposition to variance matrix.
    ///
    /// For a `p x p` variance matrix `V`, the log-Cholesky decomposition is defined as its upper-Cholesky factor `V = U.transpose() * U`, but with logs of the diagonal elements, and concatenated in column-major order into a vector of length `p*(p+1)/2`.
    ///
    /// @param[in] logC The log-Cholesky factor of `V`.
    /// @return If `logC` is a vector of length `n`, the `p x p` variance matrix `V`, where `n = 1/2 * (-1 + sqrt(1 + 8n))`.
    static matrix<Type> lchol2var(vector<Type>& logC) {
      // determine size of matrix
      int p = (-1 + sqrt(1 + 8*logC.size()))/2;
      matrix<Type> C(p,p);
      C.setZero();
      int kk=0;
      // convert cholesky decomposition from vector to matrix
      for(int ii=0; ii<p; ii++) {
	for(int jj=0; jj<=ii; jj++) {
	  C(jj,ii) = logC[kk++];
	  if(ii==jj) C(ii,ii) = exp(C(ii,ii));
	}
      }
      return C.transpose() * C;
    }

    /// Compute the variance flat prior on the log-Cholesky scale.
    ///
    /// If `V` is a `p x p` variance matrix is `logC` is its log-Cholesky factor, computes `log pi(logC)` corresponding to `pi(V) ~ 1`.
    ///
    /// @param[in] logC The log-Cholesky factor of `V`.
    /// @return Scalar value of `log pi(logC)`.
    static Type lchol_prior(cRefMatrix_t logC) {
      // determine size of matrix
      int p = (-1 + sqrt(1 + 8*logC.size()))/2;
      Type lpi = Type(0.0);
      for(int ii=0; ii<p-1; ii++) {
	lpi += (p-ii-1) * logC[p*ii+ii];
      }
      return lpi;
    }
  };

} // namespace losmix
 
#endif

