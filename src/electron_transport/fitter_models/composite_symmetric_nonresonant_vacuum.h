/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file composite_symmetric_nonresonant_vacuum.h
 * \brief Fitting model for nonresonant tunneling (symmetric coupling) combined
 *    with background \"vacuum\" tunneling.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#ifndef __composite_symmetric_nonresonant_vacuum_h__
#define __composite_symmetric_nonresonant_vacuum_h__

#include <general/fitter_tools/fit_model_interface.h>
#include <gsl/gsl_integration.h>

namespace molstat {
namespace transport {

/**
 * \brief The fit model for nonresonant tunneling (symmetric coupling) through
 *    a single site, along with direct, electrode-electrode tunneling. This
 *    model was designed with molecules, combined with the background
 *    \"vacuum\" signal, in mind.
 *
 * \note This model only considers the nonresonant tunneling line shape
 *    convolved with the vacuum line shape. A bare vacuum component is not
 *    included.
 *
 * The line shape for nonresonant tunneling through a single site with symmetric
 * coupling to the leads is \cite williams-5937
 * \f[
 * \hat{P}_\mathrm{NR}(g) = \frac{N}{\sqrt{g(1-g)^3}} \exp\left[ - \frac{\left( c\sqrt{g} - d\sqrt{1-g} \right)^2}{2(1-g)} \right],
 * \f]
 * where \f$g\f$ is the conductance in atomic units. Similarly, the line shape
 * for the background tunneling is
 * \f[
 * \hat{P}_\mathrm{vacuum}(g) = \frac{N \Theta(g-g_-)}{g},
 * \f]
 * where \f$g_-\f$ is a lower bound for the conductance (perhaps experimental
 * resolution).
 *
 * As described in Reference \cite reuter-2243, the line shape through both
 * channels simultaneously is
 * \f[
 * \hat{P}_\mathrm{NR+vacuum}(g) = N \int\limits_0^g \mathrm{d}g' \hat{P}_\mathrm{NR}(g-g') \hat{P}_\mathrm{vacuum}(g').
 * \f]
 * Note that this formula assumes \f$g\ll 1\f$. Ultimately, the fitting
 * parameters are \f$c\f$, the average level alignment of the nonresonant
 * channel relative to the standard deviation in the coupling;
 * \f$d\f$, the average coupling element of the nonresonant channel relative
 * to the standard deviation in the coupling; \f$g_-\f$, the threshold
 * conductance for vacuum tunneling; and \f$N\f$, the normalization constant.
 *
 * \if fullref
 * Additional details on the nonresonant tunneling model can be found at
 * molstat::transport::SymmetricNonresonantFitModel.
 * \endif
 */
class CompositeSymmetricNonresonantVacuumFitModel : public FitModel<1>
{
protected:
	/// Maximum number of quadrature points for the adaptive GSL routines.
	const std::size_t nquad = 2000;

	/// Integration workspace for GSL numerical integration.
	std::unique_ptr<gsl_integration_workspace,
	                decltype(&gsl_integration_workspace_free)> w;

	/**
	 * \brief Converts a map of names to values to an initial guess (ordered
	 *    vector).
	 *
	 * `c`, `d`, and `gminus` parameters are needed for this model. See
	 * FitModel::create_initial_guess for more information.
	 *
	 * \throw invalid_argument_exception if `c`, `d`, and `gminus`
	 *    parameters are not specified.
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
	 * \int\limits_{g_-}^{g} \mathrm{d}g' \frac{1}{g' \sqrt{(g-g')(1-g+g')^3}} \exp\left[ - \frac{\left( c\sqrt{g-g'} - d\sqrt{1-g+g'} \right)^2}{2(1-g+g')} \right].
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
	 *    integral with respect to \f$c\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{g_-}^{g} \mathrm{d}g' \frac{d\sqrt{1-g+g'}-c\sqrt{g-g'}}{g' (1-g+g')^{5/2}} \exp\left[ - \frac{\left( c\sqrt{g-g'} - d\sqrt{1-g+g'} \right)^2}{2(1-g+g')} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_dp_dc(double gp, void *params);

	/**
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$d\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{g_-}^{g} \mathrm{d}g' \frac{c\sqrt{g-g'}-d\sqrt{1-g+g'}}{g' \sqrt{g-g'} (1-g+g')^2} \exp\left[ - \frac{\left( c\sqrt{g-g'} - d\sqrt{1-g+g'} \right)^2}{2(1-g+g')} \right].
	 * \f]
	 *
	 * \param[in] gp The current value of \f$g'\f$.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of \f$g'\f$.
	 */
	static double int_dp_dd(double gp, void *params);

public:
	/// Index for the \f$c\f$ fitting parameter.
	const static int C = 0;

	/// Index for the \f$d\f$ fitting parameter.
	const static int D = 1;

	/// Index for the \f$g_-\f$ fitting parameter.
	const static int GMINUS = 2;

	/// Index for the norm fitting parameter.
	const static int NORM = 3;

	CompositeSymmetricNonresonantVacuumFitModel() = delete;
	virtual ~CompositeSymmetricNonresonantVacuumFitModel() = default;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	CompositeSymmetricNonresonantVacuumFitModel(
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
	 * \frac{\partial \hat{P}(g)}{\partial c} & = & N \mathrm{int\_dp\_dc}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial d} & = & N \mathrm{int\_dp\_dd}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial g_-} & = & \frac{-N}{g_- \sqrt{(g-g_-)(1-g+g_-)^3}} \exp\left[ -\frac{\left( c\sqrt{g-g_-} - d\sqrt{1-g+g_-} \right)^2}{2(1-g+g_-)} \right], \\
	 * \frac{\partial \hat{P}(g)}{\partial N} & = & \mathrm{int\_p}(g).
	 * \f}
	 * where int_p is evaluated with CompositeSymmetricNonresonantVacuumFitModel::int_p,
	 * and likewise for CompositeSymmetricNonresonantVacuumFitModel::int_dp_dc and
	 * CompositeSymmetricNonresonantVacuumFitModel::int_dp_dd.
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
	 * the integral in CompositeSymmetricNonresonantVacuumFitModel::int_p.
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
