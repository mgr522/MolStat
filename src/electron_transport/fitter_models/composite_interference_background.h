/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file composite_interference_background.h
 * \brief Fitting model for interference transport combined with background
 *    tunneling.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#ifndef __composite_interference_vacuum_h__
#define __composite_interference_vacuum_h__

#include <general/fitter_tools/fit_model_interface.h>
#include <gsl/gsl_integration.h>

namespace molstat {
namespace transport {

/**
 * \brief The fit model for an interference feature, along with direct,
 *    electrode-electrode tunneling.
 *
 * \note This model only considers the interference tunneling line shape
 *    convolved with the background line shape. A bare background component is
 *    not included.
 *
 * The line shape for transport around an interference feature is
 * \f[
 * \hat{P}_\mathrm{interference}(g) = \frac{N}{\sqrt{g}} \exp\left[ - \frac{c_\omega^2g}{2} \right],
 * \f]
 * where \f$g\f$ is the conductance in atomic units and \f$c_\omega\f$ is the
 * fitting parameter (the mean \"steepness\" of the interference feature
 * relative to the standard deviation in level alignment). Similarly, the line
 * shape for the background tunneling is
 * \f[
 * \hat{P}_\mathrm{background}(g) = \frac{N \Theta(g-g_-)}{g},
 * \f]
 * where \f$g_-\f$ is a lower bound for the conductance (perhaps experimental
 * resolution). To make things smoother for computation, we take
 * \f[
 * \hat{P}_\mathrm{vacuum}(g) = \frac{N}{g} \left[ 1 + \mathrm{erf}\left( \frac{g-g_-}{kg_-} \right) \right],
 * \f]
 * where \f$0<k\ll 1\f$. (This essentially replaces the step function by erf
 * with a \"standard deviation\" of \f$kg_-\f$.) \f$k\f$ is not regarded as a
 * fitting parameter; it is taken as a constant.
 *
 * As described in Reference \cite reuter-2243, the line shape through both
 * channels simultaneously is
 * \f[
 * \hat{P}_\mathrm{interference\ast background}(g) = N \int\limits_0^g \mathrm{d}g' \hat{P}_\mathrm{interference}(g-g') \hat{P}_\mathrm{background}(g').
 * \f]
 * Note that this formula assumes \f$g\ll 1\f$. Ultimately, the fitting
 * parameters are \f$c_\omega\f$, described above; \f$g_-\f$, the threshold
 * conductance for vacuum tunneling; and \f$N\f$, the normalization constant.
 *
 * Since this functional form is singular as \f$g\to0\f$, the line shape spans
 * several orders of magnitude, possibly causing problems during the fit.
 * Specifically, the fit essentially overweights the points close to \f$g=0\f$;
 * their higher magnitudes leads to larger absolute residuals. To compensate,
 * we scale the residuals by the magnitude of the observed value; this, in
 * effect, provides a more equal weighting for all points and provides more
 * accurate fits.
 */
class CompositeInterferenceBackgroundFitModel : public FitModel<1>
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
	 * `comega` and `gminus` parameters are needed for this model. See
	 * FitModel::create_initial_guess for more information.
	 *
	 * \throw invalid_argument_exception if `comega` and `gminus` parameters are
	 *    not specified.
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
	 * \int\limits_{0}^{g} \mathrm{d}g' \frac{1+\mathrm{erf}\left( \frac{g'-g_-}{kg_-} \right)}{g' \sqrt{g-g'}} \exp\left[ - \frac{c_\omega^2(g-g')}{2} \right].
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
	 *    integral with respect to \f$c_\omega\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{0}^{g} \mathrm{d}g' \frac{\left[ 1 + \mathrm{erf}\left( \frac{g'-g_-}{kg_-} \right) \right]\sqrt{g-g'}}{g'} \exp\left[ - \frac{c_\omega^2(g-g')}{2} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_dp_dcomega(double gp, void *params);

	/**
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$g_-\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{0}^{g} \mathrm{d}g' \frac{1}{\sqrt{g-g'}} \exp\left[ -\frac{c_\omega^2 (g-g')}{2} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_dp_dgminus(double gp, void *params);

public:
	/// Index for the \f$c_\omega\f$ fitting parameter.
	const static int COMEGA = 0;

	/// Index for the \f$g_-\f$ fitting parameter.
	const static int GMINUS = 1;

	/// Index for the norm fitting parameter.
	const static int NORM = 2;

	CompositeInterferenceBackgroundFitModel() = delete;
	virtual ~CompositeInterferenceBackgroundFitModel() = default;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	CompositeInterferenceBackgroundFitModel(
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
	 * \frac{\partial \hat{P}(g)}{\partial c_\omega} & = & -c_\omega N \mathrm{int\_dp\_dcomega}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial g_-} & = & \frac{-N}{g_- \sqrt{g-g_-}} \exp\left[ -\frac{c_\omega^2(g-g_-)}{2} \right], \\
	 * \frac{\partial \hat{P}(g)}{\partial N} & = & \mathrm{int\_p}(g).
	 * \f}
	 * where int_p is evaluated with \c CompositeInterferenceVacuumFitModel::int_p,
	 * and likewise for \c CompositeInterferenceVacuumFitModel::int_dp_dcomega.
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
	 * The residual and the derivative with respect to \f$N\f$ both require
	 * the integral in \c CompositeInterferenceVacuumFitModel::int_p.
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
};

} // namespace molstat::transport
} // namespace molstat

#endif
