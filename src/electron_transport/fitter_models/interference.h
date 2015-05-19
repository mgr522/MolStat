/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file interference.h
 * \brief Fitting model for transport around an interference feature.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#ifndef __interference_h__
#define __interference_h__

#include <general/fitter_tools/fit_model_interface.h>

namespace molstat {
namespace transport {

/**
 * \brief The fit model for transport around an interference feature.
 *
 * The line shape for transport near an interference feature is
 * \f[
 * \hat{P}(g) = \frac{N}{\sqrt{g}} \exp\left[ - \frac{c_\omega^2g}{2} \right],
 * \f]
 * where \f$g\f$ is the conductance in atomic units. The fitting parameters
 * are \f$c_\omega\f$, the steepness of the interference well relative to the
 * standard deviation in the level alignment, and \f$N\f$, the normalization
 * constant.
 *
 * Since this functional form is singular as \f$g\to0\f$, the line shape spans
 * several orders of magnitude, possibly causing problems during the fit.
 * Specifically, the fit essentially overweights the points close to \f$g=0\f$;
 * their higher magnitudes leads to larger absolute residuals. To compensate,
 * we scale the residuals by the magnitude of the observed value; this, in
 * effect, provides a more equal weighting for all points and provides more
 * accurate fits.
 */
class InterferenceFitModel : public FitModel<1>
{
protected:
	/**
	 * \brief Converts a map of names to values to an initial guess (ordered
	 *    vector).
	 *
	 * The `comega` parameter is needed for this model. See
	 * \c FitModel::create_initial_guess for more information.
	 *
	 * \throw invalid_argument_exception if the `comega` parameter is not
	 *    specified.
	 *
	 * \param[in] values The map of names to values.
	 * \return A vector containing the initial guess.
	 */
	virtual std::vector<double> create_initial_guess(
		const std::map<std::string, double> &values) const override;

public:
	/// Index for the \f$c_\omega\f$ fitting parameter.
	const static int COMEGA = 0;

	/// Index for the norm fitting parameter.
	const static int NORM = 1;

	InterferenceFitModel() = delete;
	virtual ~InterferenceFitModel() = default;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	InterferenceFitModel(
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
	 * \frac{\partial \hat{P}(g)}{\partial c_\omega} & = & -N c_\omega \sqrt{g} \exp\left[ -\frac{c_\omega^2 g}{2} \right]; \\
	 * \frac{\partial \hat{P}(g)}{\partial N} & = & \frac{1}{\sqrt{g}} \exp\left[ -\frac{c_\omega^2 g}{2} \right].
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
	 * Makes sure \f$c_\omega\f$ is positive.
	 *
	 * \param[in,out] fitparams The fitting parameters.
	 */
	virtual void process_fit_parameters(std::vector<double> &fitparams) const
		override;
};

} // namespace transport
} // namespace molstat

#endif
