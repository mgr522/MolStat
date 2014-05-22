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

/**
 * \brief The fit model for nonresonant tunneling through a single site with
 *    symmetric electrode/site couplings.
 *
 * The line shape for nonresonant tunneling through a single site with
 * symmetric coupling to the leads is \f[
 * \hat{P}(g) = \frac{N}{\sqrt{g(1-g)^3}} \exp\left[ - \frac{(c \sqrt{g} - d \sqrt{1-g})^2}{2(1-g)} \right],
 * \f] where \f$g\f$ is the conductance in atomic units. The fitting parameters
 * are \f$c\f$, the average level alignment relative to the standard deviation
 * in the channel coupling; \f$d\f$, the average channel coupling relative to
 * its standard deviation; and \f$N\f$, the normalization constant.
 *
 * This form is based on a symmetric-coupling single-site model
 * (SymmetricOneSiteModel) with normal distributions (NormalDistribution) on
 * \f$\varepsilon\f$ and \f$\Gamma\f$.
 *
 * This model is detailed in Reference \cite williams-5937.
 */
class SymmetricNonresonantFitModel : public FitModel<1> {
public:
	/**
	 * \brief Index for the \f$c\f$ fitting parameter.
	 */
	const static int C = 0;

	/**
	 * \brief Index for the \f$d\f$ fitting parameter.
	 */
	const static int D = 1;

	/**
	 * \brief Index for the norm fitting parameter.
	 */
	const static int NORM = 2;

	SymmetricNonresonantFitModel() = delete;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	SymmetricNonresonantFitModel(
		const std::list<std::pair<std::array<double, 1>, double>> &data);

	/**
	 * \brief Destructor.
	 */
	virtual ~SymmetricNonresonantFitModel() = default;

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
		const std::array<double, 1> &x, const double f) const;

	/**
	 * \brief Calculates the Jacobian of the fit function for a given set of
	 *    fitting parameters and a specific set of independent variables.
	 *
	 * For reference, the elements of the Jacobian are given by
	 * \f{eqnarray}
	 * \frac{\partial \hat{P}(g)}{\partial c} & = & -\frac{N (c\sqrt{g} - d\sqrt{1-g})}{(1-g)^{5/2}} \exp\left[ - \frac{(c\sqrt{g} - d\sqrt{1-g})^2}{2(1-g)} \right]; \\
	 * \frac{\partial \hat{P}(g)}{\partial d} & = & \frac{N (c\sqrt{g} - d\sqrt{1-g})}{\sqrt{g} (1-g)^2} \exp\left[ - \frac{(c\sqrt{g} - d\sqrt{1-g})^2}{2(1-g)} \right]; \\
	 * \frac{\partial \hat{P}(g)}{\partial N} & = & \frac{1}{\sqrt{g(1-g)3}} \exp\left[ - \frac{(c\sqrt{g} - d\sqrt{1-g})^2}{2(1-g)} \right]; \\
	 * \f}
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return The Jacobian evaluated at these independent variables and
	 *    fitting parameters.
	 */
	virtual std::vector<double> jacobian(const std::vector<double> &fitparam,
		const std::array<double, 1> &x, const double f) const;

	/* There is no redundency in calculating fit_function and jacobian, so we
	 * can use the simple default resid_j function in FitModel<1>. */

	/**
	 * \brief Returns a list of initial guesses to use when fitting the data.
	 *
	 * \return A list of initial guesses.
	 */
	virtual std::list<std::vector<double>> initial_guesses() const;

	/**
	 * \brief Prints the fit variables from a gsl_vector.
	 *
	 * \param[in] f The output stream.
	 * \param[in] fitparam The fitting parameters.
	 */
	virtual void print_fit(FILE *f, const std::vector<double> &fitparam)
		const;

	/* No post-processing issues encountered. We can just use the default
	 * function that does nothing. */
};

#endif
