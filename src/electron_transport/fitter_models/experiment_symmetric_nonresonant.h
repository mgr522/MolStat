/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file experiment_symmetric_nonresonant.h
 * \brief Fitting model for nonresonant tunneling (symmetric coupling) combined
 *    with background tunneling and a pure background component.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#ifndef __experiment_symmetric_nonresonant_h__
#define __experiment_symmetric_nonresonant_h__

#include <general/fitter_tools/fit_model_interface.h>
#include <gsl/gsl_integration.h>

namespace molstat {
namespace transport {

/**
 * \brief The fit model for nonresonant tunneling (symmetric coupling) through
 *    a single site, along with direct, electrode-electrode tunneling. This
 *    model was designed with molecules, combined with the background
 *    signal, in mind.
 *
 * \note This model includes the convolution of the symmetric nonresonant model
 *    with the background model as well as a component for background only.
 *
 * The line shape for nonresonant tunneling through a single site with symmetric
 * coupling to the leads is \cite williams-5937
 * \f[
 * \hat{P}_\mathrm{NR}(g) = \frac{N_\mathrm{NR}}{\sqrt{g(1-g)^3}} \exp\left[ - \frac{\left( c_\varepsilon\sqrt{g} - c_\Gamma\sqrt{1-g} \right)^2}{2(1-g)} \right],
 * \f]
 * where \f$g\f$ is the conductance in atomic units. Similarly, the line shape
 * for the background tunneling is
 * \f[
 * \hat{P}_\mathrm{background}(g) = \frac{N_\mathrm{background} \Theta(g-g_-)}{g},
 * \f]
 * where \f$g_-\f$ is a lower bound for the conductance (perhaps experimental
 * resolution). To make things smoother for computation, we take
 * \f[
 * \hat{P}_\mathrm{background}(g) = \frac{N_\mathrm{background}}{g} \left[ 1 + \mathrm{erf}\left( \frac{g-g_-}{kg_-} \right) \right],
 * \f]
 * where \f$0<k\ll 1\f$. (This essentially replaces the step function by erf
 * with a \"standard deviation\" of \f$kg_-\f$.) \f$k\f$ is not regarded as a
 * fitting parameter; it is taken as a constant.
 *
 * As described in Reference \cite reuter-2243, the line shape through both
 * channels simultaneously is
 * \f[
 * \hat{P}_\mathrm{NR\ast background}(g) = N_\mathrm{signal} \int\limits_0^g \mathrm{d}g' \hat{P}_\mathrm{NR}(g-g') \hat{P}_\mathrm{background}(g') + \frac{N_\mathrm{background}}{g} + N_\mathrm{baseline}.
 * \f]
 * Note that this formula assumes \f$g\ll 1\f$. Ultimately, the fitting
 * parameters are \f$c_\varepsilon\f$, the average level alignment of the
 * nonresonant channel relative to the standard deviation in the coupling;
 * \f$c_\Gamma\f$, the average coupling element of the nonresonant channel
 * relative to the standard deviation in the coupling; \f$g_-\f$, the threshold
 * conductance for vacuum tunneling; \f$N_\mathrm{signal}\f$ and
 * \f$N_\mathrm{background}\f$ are normalization/scaling constants; and
 * \f$N_\mathrm{baseline}\f$ is a constant.
 *
 * \if fullref
 * Additional details on the nonresonant tunneling model can be found at
 * \c molstat::transport::SymmetricNonresonantFitModel.
 * \endif
 */
class ExperimentSymmetricNonresonantFitModel
	: public FitModel<1>
{
protected:
	/// Maximum number of quadrature points for the adaptive GSL routines.
	const std::size_t nquad = 2000;

	/// The effective standard deviation of the step function smoothing.
	constexpr static double k = 0.05;

	/// Integration workspace for GSL numerical integration.
	std::unique_ptr<gsl_integration_workspace,
	                decltype(&gsl_integration_workspace_free)> w;

	/**
	 * \brief Converts a map of names to values to an initial guess (ordered
	 *    vector).
	 *
	 * `cepsilon`, `cgamma`, `gminus`, `nsignal`, `nbackground`, and `nbaseline`
	 * parameters are needed for this model. See \c FitModel::create_initial_guess
	 * for more information.
	 *
	 * \throw invalid_argument_exception if `cepsilon`, `cgamma`, `gminus`,
	 *    `nsignal`, `nbackground`, `nbaseline` parameters are not specified.
	 *
	 * \param[in] values The map of names to values.
	 * \return A vector containing the initial guess.
	 */
	virtual std::vector<double> create_initial_guess(
		const std::map<std::string, double> &values) const override;

	/**
	 * \brief GSL integrand function for the fit function integral.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{0}^{g} \mathrm{d}g' \frac{1+\mathrm{erf}\left( \frac{g'-g_-}{kg_-} \right)}{g' \sqrt{(g-g')(1-g+g')^3}} \exp\left[ - \frac{\left( c_\varepsilon\sqrt{g-g'} - c_\Gamma\sqrt{1-g+g'} \right)^2}{2(1-g+g')} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_p(double gp, void *params);

	/**
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$c_\varepsilon\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{0}^{g} \mathrm{d}g' \frac{\left[ 1+\mathrm{erf}\left( \frac{g'-g_-}{kg_-} \right) \right] \left(c_\Gamma\sqrt{1-g+g'}-c_\varepsilon\sqrt{g-g'}\right)}{g' (1-g+g')^{5/2}} \exp\left[ - \frac{\left( c_\varepsilon\sqrt{g-g'} - c_\Gamma\sqrt{1-g+g'} \right)^2}{2(1-g+g')} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_dp_dcepsilon(double gp, void *params);

	/**
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$c_\Gamma\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{0}^{g} \mathrm{d}g' \frac{\left[ 1+\mathrm{erf}\left( \frac{g'-g_-}{kg_-} \right) \right] \left( c_\varepsilon\sqrt{g-g'}-c_\Gamma\sqrt{1-g+g'}\right)}{g' \sqrt{g-g'} (1-g+g')^2} \exp\left[ - \frac{\left( c_\varepsilon\sqrt{g-g'} - c_\Gamma\sqrt{1-g+g'} \right)^2}{2(1-g+g')} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_dp_dcgamma(double gp, void *params);

	/**
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$g_-\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{0}^{g} \mathrm{d}g' \frac{1}{\sqrt{(g-g')(1-g+g')^3}} \exp\left[ -\left( \frac{g'-g_-}{kg_-} \right)^2 - \frac{\left( c_\varepsilon\sqrt{g-g'} - c_\Gamma\sqrt{1-g+g'} \right)^2}{2(1-g+g')} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_dp_dgminus(double gp, void *params);

public:
	/// Index for the \f$c_\varepsilon\f$ fitting parameter.
	const static int CEPSILON = 0;

	/// Index for the \f$c_\Gamma\f$ fitting parameter.
	const static int CGAMMA = 1;

	/// Index for the \f$g_-\f$ fitting parameter.
	const static int GMINUS = 2;

	/// Index for the norm (composite component) fitting parameter.
	const static int NSIGNAL = 3;

	/// Index for the norm (background component) fitting parameter.
	const static int NBACKGROUND = 4;

	/// Index for the additive constant (baseline) fitting parameter.
	const static int NBASELINE = 5;

	ExperimentSymmetricNonresonantFitModel() = delete;
	virtual ~ExperimentSymmetricNonresonantFitModel() = default;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	ExperimentSymmetricNonresonantFitModel(
		const std::list<std::pair<std::array<double, 1>, double>> &data);

	/**
	 * \brief Evaluates the fit function for this model at a given set of
	 *    independent variables and fitting parameters.
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return The function evaluated at these independent variables and
	 *    fitting parameters.
	 */
	virtual double resid(const std::vector<double> &fitparam,
		const std::array<double, 1> &x, const double f) const override;

	/**
	 * \brief Calculates the Jacobian of the fit function for a given set of
	 *    fitting parameters and a specific set of independent variables.
	 *
	 * Each element of the Jacobian depends on a particular integral, as
	 * follows.
	 * \f{eqnarray}
	 * \frac{\partial \hat{P}(g)}{\partial c_\varepsilon} & = & N_\mathrm{signal} \mathrm{int\_dp\_dcepsilon}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial c_\Gamma} & = & N_\mathrm{signal} \mathrm{int\_dp\_dcgamma}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial g_-} & = & \frac{-N_\mathrm{signal}}{g_- \sqrt{(g-g_-)(1-g+g_-)^3}} \exp\left[ -\frac{\left( c_\varepsilon\sqrt{g-g_-} - c_\Gamma\sqrt{1-g+g_-} \right)^2}{2(1-g+g_-)} \right], \\
	 * \frac{\partial \hat{P}(g)}{\partial N_\mathrm{signal}} & = & \mathrm{int\_p}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial N_\mathrm{background}} & = & \frac{1}{g}, \\
	 * \frac{\partial \hat{P}(g)}{\partial N_\mathrm{baseline}} & = & 1.
	 * \f}
	 * where int_p is evaluated with \c SymmetricNonresonantPlusBackgroundFitModel::int_p,
	 * and likewise for \c SymmetricNonresonantPlusBackgroundFitModel::int_dp_dcepsilon
	 * and \c SymmetricNonresonantPlusBackgroundFitModel::int_dp_dcgamma.
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return The Jacobian evaluated at these independent variables and
	 *    fitting parameters.
	 */
	virtual std::vector<double> jacobian(const std::vector<double> &fitparam,
		const std::array<double, 1> &x, const double f) const override;

	/**
	 * \brief Calculates both the residual and Jacobian of the fit for a given
	 *    set of fitting parameters.
	 *
	 * The residual and the derivative with respect to \f$N_\mathrm{signal}\f$
	 * both require the integral in \c SymmetricNonresonantPlusBackgroundFitModel::int_p.
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return A pair of the residual and Jacobian evaluated at these
	 *    independent variables and fitting parameters.
	 */
	virtual std::pair<double, std::vector<double>> resid_j(
		const std::vector<double> &fitparam, const std::array<double, 1> &x,
		const double f) const override;

	virtual void append_default_guesses(std::list<std::vector<double>> &guess)
		const override;

	virtual void print_fit(std::ostream &out,
		const std::vector<double> &fitparam) const override;

	/**
	 * \brief Perform post-processing on a set of fit parameters.
	 *
	 * \param[in,out] fitparams The fitting parameters.
	 */
	virtual void process_fit_parameters(std::vector<double> &fitparams) const
		override;

	virtual bool is_good_fit(const std::vector<double> &fitparams) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
