/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file symmetric_nonresonant.h
 * \brief Fitting model for nonresonant tunneling (symmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __symmetric_nonresonant_h__
#define __symmetric_nonresonant_h__

#include <general/fitter_tools/fit_model_interface.h>

namespace molstat {
namespace transport {

/**
 * \brief The fit model for nonresonant tunneling through a single site with
 *    symmetric electrode/site couplings.
 *
 * The line shape for nonresonant tunneling through a single site with
 * symmetric coupling to the leads is \f[
 * \hat{P}(g) = \frac{N}{\sqrt{g(1-g)^3}} \exp\left[ - \frac{(c_\varepsilon \sqrt{g} - c_\Gamma \sqrt{1-g})^2}{2(1-g)} \right],
 * \f] where \f$g\f$ is the conductance in atomic units. The fitting parameters
 * are \f$c_\varepsilon\f$, the average level alignment relative to the
 * standard deviation in the channel coupling; \f$c_\Gamma\f$, the average
 * channel coupling relative to its standard deviation; and \f$N\f$, the
 * normalization constant.
 *
 * This form is based on a symmetric-coupling single-site model
 * (\c SymmetricOneSiteModel) with normal distributions (\c NormalDistribution)
 * on \f$\varepsilon\f$ and \f$\Gamma\f$.
 *
 * This model is detailed in Reference \cite williams-5937.
 */
class SymmetricNonresonantFitModel : public FitModel<1>
{
protected:
	/**
	 * \brief Converts a map of names to values to an initial guess (ordered
	 *    vector).
	 *
	 * `cepsilon` and `cgamma` parameters are needed for this model. See
	 * FitModel::create_initial_guess for more information.
	 *
	 * \throw invalid_argument_exception if `cepsilon` and `cgamma` parameters 
	 *    are not specified.
	 *
	 * \param[in] values The map of names to values.
	 * \return A vector containing the initial guess.
	 */
	virtual std::vector<double> create_initial_guess(
		const std::map<std::string, double> &values) const override;

public:
	/// Index for the \f$c_\varepsilon\f$ fitting parameter.
	const static int CEPSILON = 0;

	/// Index for the \f$c_\Gamma\f$ fitting parameter.
	const static int CGAMMA = 1;

	/// Index for the norm fitting parameter.
	const static int NORM = 2;

	SymmetricNonresonantFitModel() = delete;
	virtual ~SymmetricNonresonantFitModel() = default;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	SymmetricNonresonantFitModel(
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
	 * For reference, the elements of the Jacobian are given by
	 * \f{eqnarray}
	 * \frac{\partial \hat{P}(g)}{\partial c_\varepsilon} & = & -\frac{N (c_\varepsilon\sqrt{g} - c_\Gamma\sqrt{1-g})}{(1-g)^{5/2}} \exp\left[ - \frac{(c_\varepsilon\sqrt{g} - c_\Gamma\sqrt{1-g})^2}{2(1-g)} \right]; \\
	 * \frac{\partial \hat{P}(g)}{\partial c_\Gamma} & = & \frac{N (c_\varepsilon\sqrt{g} - c_\Gamma\sqrt{1-g})}{\sqrt{g} (1-g)^2} \exp\left[ - \frac{(c_\varepsilon\sqrt{g} - c_\Gamma\sqrt{1-g})^2}{2(1-g)} \right]; \\
	 * \frac{\partial \hat{P}(g)}{\partial N} & = & \frac{1}{\sqrt{g(1-g)^3}} \exp\left[ - \frac{(c_\varepsilon\sqrt{g} - c_\Gamma\sqrt{1-g})^2}{2(1-g)} \right]; \\
	 * \f}
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return The Jacobian evaluated at these independent variables and
	 *    fitting parameters.
	 */
	virtual std::vector<double> jacobian(const std::vector<double> &fitparam,
		const std::array<double, 1> &x, const double f) const override;

	/* There is no redundency in calculating fit_function and jacobian, so we
	 * can use the simple default resid_j function in FitModel<1>. */

	virtual void append_default_guesses(std::list<std::vector<double>> &guess)
		const override;

	virtual void print_fit(std::ostream &out,
		const std::vector<double> &fitparam) const override;

	virtual void process_fit_parameters(std::vector<double> &fitparams) const
		override;

	virtual bool is_good_fit(const std::vector<double> &fitparams) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
