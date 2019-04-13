/// @file mNIX.hpp
///
/// TODO:
///
/// - reorder lambda, Omega, tau, nu?? (Probably not.)
/// - more efficient computation of normalizing constant. 
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
    // Cholesky solvers
    Eigen::LLT<MatrixXd_t> llt_; // for Omega_hat
    Eigen::LLT<MatrixXd_t> lltx_; // for external Omega
    // simulation
    matrix<Type> z_;
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
    void calc_post();
    /// Normalizing constant for mNIX distribution.
    Type zeta(const Eigen::LLT<MatrixXd_t>& llt,
	      const Type& nu, const Type& tau);
    /// Normalizing constant for mNIX distribution.
    Type zeta(cRefMatrix_t& Omega,
	      const Type& nu, const Type& tau);
    /// Normalizing constant for mNIX distribution.
    Type zeta();
    /// Marginal likelihood for mNIX distribution.
    Type log_marg();
    /// Random draw from the mNIX distribution.
    void simulate(RefMatrix_t beta, Type& sigma,
		  cRefMatrix_t& lambda,
		  const Eigen::LLT<MatrixXd_t>& llt,
		  const Type& nu, const Type& tau);
    /// Random draw from the mNIX distribution.
    void simulate(RefMatrix_t beta, Type& sigma,
		  cRefMatrix_t& lambda,
		  cRefMatrix_t& Omega,
		  const Type& nu, const Type& tau);
    /// Random draw from the mNIX distribution.
    void simulate(RefMatrix_t beta, Type& sigma);
    /// Constructor.
    mNIX(int p);    
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
    Ol_ = utils<Type>::zero_matrix(p_,1);
    Ol_hat_ = utils<Type>::zero_matrix(p_,1);
    // cholesky solvers
    llt_.compute(utils<Type>::identity_matrix(p_,p_));
    lltx_.compute(utils<Type>::identity_matrix(p_,p_));
    // simulation
    z_ = utils<Type>::zero_matrix(p_,1);
  }

  /// @param[in] lambda Prior mean vector, as a `p x 1` matrix.
  /// @param[in] Omega Prior precision matrix of size `p x p`.
  /// @param[in] nu Prior shape parameter (scalar).
  /// @param[in] tau Prior scale parameter (scalar).
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
    llt_.compute(Omega_hat);
    lambda_hat = llt_.solve(Ol_ + Xy_);
    nu_hat = nu_ + N_;
    Ol_hat_ = Omega_hat * lambda_hat;
    tau_hat = yy_ - utils<Type>::dot_product(lambda_hat, Ol_hat_) + lOl_ + nu_*tau_;
    tau_hat /= nu_hat;
    return;
  }

  /// @warning Must be run after a call to `set_suff` and `set_prior`.
  template <class Type>
  inline void mNIX<Type>::calc_post() {
    get_post(lambda_hat_, Omega_hat_, nu_hat_, tau_hat_);
    return;
  }

  /// @param[out] beta Regression coefficients, as a `p x 1` matrix.
  /// @param[out] sigma Standard deviation of errors, as a positive scalar.
  /// @param[in] lambda Mean vector, as a `p x 1` matrix.
  /// @param[in] llt Cholesky decomposition of the precision matrix of size `p x p`.  Must be precomputed.
  /// @param[in] nu Shape parameter (positive scalar).
  /// @param[in] tau Scale parameter (positive scalar).
  template <class Type>
  inline void mNIX<Type>::simulate(RefMatrix_t beta, Type& sigma,
  				   cRefMatrix_t& lambda,
  				   const Eigen::LLT<MatrixXd_t>& llt,
  				   const Type& nu, const Type& tau) {
    sigma = 1.0/sqrt(rgamma(.5*nu, 2.0/(nu*tau)));
    z_ = rnorm(p_, Type(0),sigma).matrix();
    beta = llt.matrixU().solve(z_) + lambda;
    return;
  }

  /// @param[out] beta Regression coefficients, as a `p x 1` matrix.
  /// @param[out] sigma Standard deviation of errors, as a positive scalar.
  /// @param[in] lambda Mean vector, as a `p x 1` matrix.
  /// @param[in] Omega Precision matrix of size `p x p`.
  /// @param[in] nu Shape parameter (positive scalar).
  /// @param[in] tau Scale parameter (positive scalar).
  template <class Type>
  inline void mNIX<Type>::simulate(RefMatrix_t beta, Type& sigma,
  				   cRefMatrix_t& lambda,
  				   cRefMatrix_t& Omega,
  				   const Type& nu, const Type& tau) {
    lltx_.compute(Omega);
    simulate(beta, sigma, lambda, lltx_, nu, tau);
    return;
  }

  /// @param[out] beta Regression coefficients, as a `p x 1` matrix.
  /// @param[out] sigma Standard deviation of errors, as a positive scalar.
  ///
  /// @warning Must be run after a call to calc_post.
  template <class Type>
  inline void mNIX<Type>::simulate(RefMatrix_t beta, Type& sigma) {
    simulate(beta, sigma, lambda_hat_, llt_, nu_hat_, tau_hat_);
    return;
  }

  /// @param[in] llt Cholesky decomposition of the precision matrix Omega of size `p x p`.  Must be precomputed.
  /// @param[in] nu Shape parameter (positive scalar).
  /// @param[in] tau Scale parameter (positive scalar).
  /// @return The mNIX normalizing constant, defined as
  /// \f[
  /// 2 \log \Gamma(\nu/2) - nu \log(\tau\nu/2) - \log |\Omega|.
  /// \f]
  template <class Type>
  Type mNIX<Type>::zeta(const Eigen::LLT<MatrixXd_t>& llt,
			const Type& nu, const Type& tau) {
    // log-determinant of chol(Omega)
    Type ldC = 0.0;
    for(int ii=0; ii<p_; ii++) {
      ldC += log(llt.matrixL()(ii,ii));
    }
    Type nu2 = .5 * nu;
    return lgamma(nu2) - ldC - nu2 * log(tau * nu2);
  }
  
  /// @param[in] Omega Precision matrix of size `p x p`.
  /// @param[in] nu Shape parameter (positive scalar).
  /// @param[in] tau Scale parameter (positive scalar).
  /// @return The mNIX normalizing constant, defined as
  /// \f[
  /// 2 \log \Gamma(\nu/2) - nu \log(\tau\nu/2) - \log |\Omega|.
  /// \f]
  template <class Type>
  Type mNIX<Type>::zeta(cRefMatrix_t& Omega,
			const Type& nu, const Type& tau) {
    lltx_.compute(Omega);
    return zeta(lltx_, nu, tau);
    // matrix<Type> Omega_ = Omega;
    // return 2.0 * lgamma(.5 * nu) -
    //   nu * log(.5 * tau*nu) - atomic::logdet(Omega_);
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
    return zeta(llt_, nu_hat_, tau_hat_);
    // // log-determinant of chol(Omega_hat)
    // Type ldC = 0.0;
    // for(int ii=0; ii<p_; ii++) {
    //   ldC += log(llt_.matrixL()(ii,ii));
    // }
    // return 2.0 * (lgamma(.5 * nu_hat_) - ldC) -
    //   nu_hat_ * log(.5 * tau_hat_*nu_hat_);
  }

  /// @return The mNIX marginal (log) likelihood, defined as
  /// \f[
  /// \tfrac 1 2 \big[\zeta(\hat \phi) - \zeta(\phi)\big],
  /// \f]
  /// where `phi = (lambda, Omega, nu, tau)` and `phi_hat = (lambda_hat, Omega_hat, nu_hat, tau_hat)` are the prior and posterior mNIX hyperparameters, and `zeta()` is the mNIX normalizing constant.
  template <class Type>
  inline Type mNIX<Type>::log_marg() {
    return zeta() - zeta(Omega_, nu_, tau_);
  }

} // namespace losmix

#endif
