/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file symmetric_resonant.h
 * \brief Fitting model for resonant tunneling (symmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __symmetric_resonant_h__
#define __symmetric_resonant_h__

#include <general/fitter_tools/fit_model_interface.h>

namespace molstat {
namespace transport {

/**
 * \brief The fit model for resonant tunneling through a single site with
 *    symmetric electrode/site couplings.
 *
 * The line shape for resonant tunneling through a single site with symmetric
 * coupling to the leads is
 * \f[
 * \hat{P}(g) = \frac{N}{\sqrt{g^3(1-g)}} \exp\left[ -\gamma^2 \frac{1-g}{2g} \right],
 * \f]
 * where \f$g\f$ is the conductance in atomic units. The fitting parameters
 * are \f$\gamma\f$, the average channel coupling relative to the standard
 * deviation in the level alignment, and \f$N\f$, the normalization constant.
 *
 * This form is based on a symmetric-coupling single-site model
 * (SymmetricOneSiteModel) with normal distributions (NormalDistribution) on
 * \f$\varepsilon\f$ and \f$\Gamma\f$.
 *
 * Since this functional form is singular as \f$g\to1\f$, the line shape spans
 * several orders of magnitude, possibly causing problems during the fit.
 * Specifically, the fit essentially overweights the points close to \f$g=1\f$;
 * their higher magnitudes leads to larger absolute residuals. To compensate,
 * we scale the residuals by the magnitude of the observed value; this, in
 * effect, provides a more equal weighting for all points and provides more
 * accurate fits.
 *
 * This model is detailed in Reference \cite williams-5937.
 */
class SymmetricResonantFitModel : public FitModel<1>
{
protected:
	/**
	 * \brief Converts a map of names to values to an initial guess (ordered
	 *    vector).
	 *
	 * The `Gamma` parameter is needed for this model. See
	 * FitModel::create_initial_guess for more information.
	 *
	 * \throw invalid_argument_exception if the `Gamma` parameter is not
	 *    specified.
	 *
	 * \param[in] values The map of names to values.
	 * \return A vector containing the initial guess.
	 */
	virtual std::vector<double> create_initial_guess(
		const std::map<std::string, double> &values) const override;

public:
	/// Index for the gamma fitting parameter.
	const static int GAMMA = 0;

	/// Index for the norm fitting parameter.
	const static int NORM = 1;

	SymmetricResonantFitModel() = delete;
	virtual ~SymmetricResonantFitModel() = default;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	SymmetricResonantFitModel(
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
	 * \frac{\partial \hat{P}(g)}{\partial \gamma} & = & -\frac{N \gamma}{g^2} \sqrt{\frac{1-g}{g}} \exp\left[ -\gamma^2 \frac{1-g}{2g} \right]; \\
	 * \frac{\partial \hat{P}(g)}{\partial N} & = & \frac{1}{\sqrt{g^3(1-g)}} \exp\left[ -\gamma^2 \frac{1-g}{2g} \right].
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

	virtual void print_fit(std::ostream &out, const std::vector<double> &fitparam)
		const override;

	/**
	 * \brief Perform post-processing on a set of fit parameters.
	 *
	 * Makes sure `gamma` is positive.
	 *
	 * \param[in,out] fitparams The fitting parameters.
	 */
	virtual void process_fit_parameters(std::vector<double> &fitparams) const
		override;
};

} // namespace transport
} // namespace molstat

#endif
