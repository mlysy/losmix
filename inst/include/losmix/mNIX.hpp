/// @file mNIX.hpp
///
/// TODO:
///
/// - reorder lambda, Omega, tau, nu??
///
/// @note Depends on `TMB.hpp` which is *not* header-guarded, so don't include it here.

#ifndef losmix_mNIX_hpp
#define losmix_mNIX_hpp 1

#include "utils.hpp"

namespace losmix {

  /// The multivariate Normal-Inverse-Chi-Square (mNIX) distribution.
  ///
  /// The multivariate normal-inverse-chi-squared (mNIX) distribution on a \f$p\f$-dimensional vector \f$x\$ and scalar \f$v\f$ is defined as
  /// \f[
  /// (x, v) \sim \mathrm{mNIX}(\lambda, \Omega, \nu, \tau)
  /// \qquad \iff \qquad
  /// \begin{aligned}
  /// 1/v & \sim \mathrm{Gamma}(\nu/2, \nu\tau/2) \\
  /// x \mid v & \sim \mathcal{N}(\lambda, v \Omega^{-1}).
  /// \end{aligned}
  /// \f]
  /// For the linear regression model
  /// /f[
  /// y_i \mid \beta, \sigma & \stackrel{iid}{\sim} \mathcal N(x_i'\beta, \sigma^2),
  /// /f]
  /// the mNIX is a conjugate prior distribution, in the sense that if the prior is given by \f$(\beta, \sigma^2) \sim \mathrm{mNIX}(\lambda, \Omega, \nu, \tau)\f$, then we have \f$(\beta, \sigma^2) \mid \mathrm{mNIX}(\hat \lambda, \hat \Omega, \hat \nu, \hat \tau)\f$.
  ///
  /// @note Due to how `TMB` treats `vector<Type>` objects, this class uses (one column) matrix inputs for vectors to do the linear algebra with minimal temporary assignments.
  template <class Type>
  class mNIX {
  public:
    //typedefs
    /// Typedef equivalent to `matrix<Type>`.
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_t;
    /// Typedef equivalent to `Ref <matrix<Type> >`.
    typedef Eigen::Ref <MatrixXd_t> RefMatrix_t;
    /// Typedef equivalent to `const Ref <const matrix<Type> >`.
    typedef const Eigen::Ref <const MatrixXd_t> cRefMatrix_t;
  private:
    int N_; // number of observations
    int p_; // number of covariates
    // sufficient statistics (data)
    matrix<Type> XX_;
    matrix<Type> Xy_;
    Type yy_;
    // prior hyperparameters
    matrix<Type> lambda_;
    matrix<Type> Omega_;
    Type nu_;
    Type tau_;
    // posterior hyperparameters
    matrix<Type> lambda_hat_;
    matrix<Type> Omega_hat_;
    Type nu_hat_;
    Type tau_hat_;
    // prior precomputations
    matrix<Type> Ol_;
    Type lOl_;
    // posterior storage
    matrix<Type> Ol_hat_;
    // Cholesky solver
    Eigen::LDLT<MatrixXd_t> chol_Ohat_;
  public:
    /// Set prior parameters.
    void set_prior(cRefMatrix_t& lambda,
		   cRefMatrix_t& Omega,
		   const Type& nu, const Type& tau);
    /// Set sufficient statistics.
    void set_suff(cRefMatrix_t& y, cRefMatrix_t& Xtr);
    /// Get sufficient statistics.
    void get_suff(Type& yy, RefMatrix_t Xy, RefMatrix_t XX);
    /// Get posterior parameters.
    void get_post(RefMatrix_t lambda_hat,
		  RefMatrix_t Omega_hat,
		  Type& nu_hat, Type& tau_hat);
    /// Calculate posterior parameters to internal values.
    ///
    /// @warning Must be run after a call to `set_suff` and `set_prior`.
    void calc_post() {
      get_post(lambda_hat_, Omega_hat_, nu_hat_, tau_hat_);
      return;
    }
    /// Normalizing constant for mNIX distribution.
    Type zeta();
    /// Constructor
    mNIX(int p);
    
    /// Normalizing constant for mNIX distribution.
    ///
    /// @param[in] Omega mNIX precision matrix.
    /// @param[in] nu mNIX shape parameter.
    /// @param[in] tau mNIX scale parameter.
    /// @return The mNIX normalizing constant, defined as
    /// \f[
    /// 2 \log \Gamma(\nu/2) - nu \log(\tau\nu/2) - \log |\Omega|.
    /// \f]
    static Type zeta(cRefMatrix_t& Omega,
		     const Type& nu, const Type& tau) {
      matrix<Type> Omega_ = Omega;
      return 2.0 * lgamma(.5 * nu) -
	nu * log(.5 * tau*nu) - atomic::logdet(Omega_);
    }
  };

  /// @param[in] p Integer number of covariates.
  ///
  /// @note Explicit memory allocation is required for `calc_post` to function.
  template <class Type>
  inline mNIX<Type>::mNIX(int p) {
    p_ = p;
    // allocate memory
    // sufficient statistics
    Xy_ = utils<Type>::zero_matrix(p_,1);
    XX_ = utils<Type>::zero_matrix(p_,p_);
    // prior parameters
    Omega_ = utils<Type>::zero_matrix(p_,p_);
    lambda_ = utils<Type>::zero_matrix(p_,1);
    // posterior parameters
    lambda_hat_ = utils<Type>::zero_matrix(p_,1);
    Omega_hat_ = utils<Type>::zero_matrix(p_, p_);
  }

  /// @param[in] lambda Prior mean vector.
  /// @param[in] Omega Prior precision matrix.
  /// @param[in] nu Prior shape parameter.
  /// @param[in] tau Prior scale parameter.
  template <class Type>
  inline void mNIX<Type>::set_prior(cRefMatrix_t& lambda,
				    cRefMatrix_t& Omega,
				    const Type& nu, const Type& tau) {
    lambda_ = lambda;
    Omega_ = Omega;
    nu_ = nu;
    tau_ = tau;
    // precomputations
    Ol_ = Omega_ * lambda_;
    lOl_ = utils<Type>::dot_product(lambda_, Ol_);
    return;
  }


  /// @param[in] y Response vector of length `N`.
  /// @param[in] Xtr Transpose of design matrix, having size `p x N`.
  template <class Type>
  inline void mNIX<Type>::set_suff(cRefMatrix_t& y,
				   cRefMatrix_t& Xtr) {
    XX_ = Xtr * Xtr.transpose(); // can't use .noalias()
    Xy_ = Xtr * y;
    yy_ = utils<Type>::dot_product(y, y);
    N_ = Xtr.cols();
    return;
  }

  /// For a regression model of the form
  /// \f[
  /// y_i \stackrel{iid}{\sim} \mathcal N(x_i'\beta, \sigma^2),
  /// \f]
  /// the sufficient statistics are `XX = Xt * X`, `Xy = Xt * y`, and `yy = yt * y`.
  ///
  /// @param[out] XX Transpose-product of design matrix with itself.
  /// @param[out] Xy Transpose-product of the design matrix `X` and response vector `y`.
  /// @param[out] yy Dot-product of response vector with itself.
  template <class Type>
  inline void mNIX<Type>::get_suff(Type& yy,
				   RefMatrix_t Xy, RefMatrix_t XX) {
    XX = XX_;
    Xy = Xy_;
    yy = yy_;
    return;
  }

  /// @param[out] lambda_hat Posterior mean vector.
  /// @param[out] Omega_hat Posterior precision matrix.
  /// @param[out] nu_hat Posterior shape parameter.
  /// @param[out] tau_hat Posterior scale parameter.
  ///
  /// @warning Must be run after a call to `set_suff` and `set_prior`.
  template <class Type>
  inline void mNIX<Type>::get_post(RefMatrix_t lambda_hat,
				   RefMatrix_t Omega_hat,
				   Type& nu_hat, Type& tau_hat) {
    Omega_hat = XX_ + Omega_;
    chol_Ohat_.compute(Omega_hat);
    lambda_hat = chol_Ohat_.solve(Ol_ + Xy_);
    nu_hat = nu_ + N_;
    Ol_hat_ = Omega_hat * lambda_hat;
    tau_hat = yy_ - utils<Type>::dot_product(lambda_hat, Ol_hat_) + lOl_ + nu_*tau_;
    tau_hat /= nu_hat;
    return;
  }

  /// @return The mNIX normalizing constant, defined as
  /// \f[
  /// 2 \log \Gamma(\nu/2) - nu \log(\tau\nu/2) - \log |\Omega|,
  /// \f]
  /// applied to the internal hyperparameter values of the conjugate posterior.
  ///
  /// @warning Must be run after a call to `calc_post`.
  template <class Type>
  inline Type mNIX<Type>::zeta() {
    // log-determinant of chol(Omega_hat)
    Type ldC = 0.0;
    for(int ii=0; ii<p_; ii++) {
      ldC += log(chol_Ohat_.vectorD()(ii));
    }
    return 2.0 * lgamma(.5 * nu_hat_) - ldC -
      nu_hat_ * log(.5 * tau_hat_*nu_hat_);
  } 

} // namespace losmix

#endif
