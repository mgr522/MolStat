/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asymmetric_resonant.h
 * \brief Fitting model for resonant tunneling (asymmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __asymmetric_resonant_h__
#define __asymmetric_resonant_h__

#include <general/fitter_tools/fit_model_interface.h>
#include <gsl/gsl_integration.h>

/**
 * \brief The fit model for resonant tunneling through a single site with
 *    asymmetric electrode/site couplings.
 *
 * The line shape for resonant tunneling through a single site with asymmetric
 * coupling to the leads is \f[
 * \hat{P}(g) = \frac{N}{g^{3/2}} \int\limits_{(2-g-2\sqrt{1-g})/g}^{(2-g+2\sqrt{1-g})/g} \mathrm{d}x \frac{x(1+x^2)^{-3/2}}{\sqrt{4x-g(1+x)^2}} \left( 1 + \frac{(\gamma_\mathrm{L} + x\gamma_\mathrm{R})^2}{1+x^2} \right) \exp\left[ - \frac{(x\gamma_\mathrm{L} - \gamma_\mathrm{R})^2}{2(1+x^2)} - \frac{r^2(\gamma_\mathrm{L}^2 + \gamma_\mathrm{R}^2) [4x-g(1+x)^2]}{8g(1+x^2)} \right],
 * \f] where \f$g\f$ is the conductance in atomic units. The fitting parameters
 * are \f$\gamma_\mathrm{L/R}\f$, the average channel couplings (to the
 * left/right electrode) relative to the standard deviation in the coupling;
 * \f$r\f$, the standard deviation in the coupling relative to the standard
 * deviation in the level alignment; and \f$N\f$, the normalization constant.
 * Note that the fits are somewhat insensitive to \f$r\f$.
 *
 * This form is based on a asymmetric-coupling single-site model
 * (AsymmetricOneSiteModel) with normal distributions (NormalDistribution) on
 * \f$\varepsilon\f$, \f$\Gamma_\mathrm{L}\f$, and \f$\Gamma_\mathrm{R}\f$.
 *
 * This model is detailed in Reference \cite williams-5937.
 */
class AsymmetricResonantFitModel : public FitModel<1> {
protected:
	/**
	 * \internal
	 * \brief Maximum number of quadrature points for the adaptive quadrature
	 *    GSL routines.
	 * \endinternal
	 */
	const std::size_t nquad = 2000;

	/**
	 * \internal
	 * \brief Integration workspace for GSL numerical integration.
	 * \endinternal
	 */
	std::shared_ptr<gsl_integration_workspace> w;

	/**
	 * \brief Converts a map of names to values to an initial guess (ordered
	 *    vector).
	 *
	 * `GammaL`, `GammaR`, and `r` parameters are needed for this model. See
	 * FitModel::create_initial_guess for more information.
	 *
	 * \throw invalid_argument_exception if `GammaL`, `GammaR`, and `r`
	 *    parameters are not specified.
	 *
	 * \param[in] values The map of names to values.
	 * \return A vector containing the initial guess.
	 */
	virtual std::vector<double> create_initial_guess(
		const std::map<std::string, double> &values) const;

	/**
	 * \internal
	 * \brief GSL integrand function for the fit function integral.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{(2-g-2\sqrt{1-g})/g}^{(2-g+2\sqrt{1-g})/g} \mathrm{d}x \frac{x(1+x^2)^{-3/2}}{\sqrt{4x-g(1+x)^2}} \left( 1 + \frac{(\gamma_\mathrm{L} + x\gamma_\mathrm{R})^2}{1+x^2} \right) \exp\left[ -\frac{(x\gamma_\mathrm{L} - x\gamma_\mathrm{R})^2}{2(1+x^2)} - \frac{r^2(\gamma_\mathrm{L}^2 + \gamma_\mathrm{R}^2) [4x-g(1+x)^2]}{8g(1+x^2)} \right].
	 * \f]
	 *
	 * \param[in] x The current value of x.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of x.
	 * \endinternal
	 */
	static double int_p(double x, void *params);

	/**
	 * \internal
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$\gamma_\mathrm{L}\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f{eqnarray}
	 * \int\limits_{(2-g-2\sqrt{1-g})/g}^{(2-g+2\sqrt{1-g})/g} \mathrm{d}x && \frac{x(1+x^2)^{-5/2}}{\sqrt{4x-g(1+x)^2}} \exp\left[ -\frac{(x\gamma_\mathrm{L} - x\gamma_\mathrm{R})^2}{2(1+x^2)} - \frac{r^2(\gamma_\mathrm{L}^2 + \gamma_\mathrm{R}^2) [4x-g(1+x)^2]}{8g(1+x^2)} \right] \nonumber \\
	 * && \times \left[ (2-x^2)\gamma_\mathrm{L} + 3x\gamma_\mathrm{L} - \frac{x(\gamma_\mathrm{L}+x\gamma_\mathrm{R})^2(x\gamma_\mathrm{L}-\gamma_\mathrm{R})}{1+x^2} - \left( 1 + \frac{(\gamma_\mathrm{L} + x\gamma_\mathrm{R})^2}{1+x^2} \right) \frac{r^2\gamma_\mathrm{L}[4x-g(1+x)^2]}{4g} \right] \nonumber
	 * \f}
	 *
	 * \param[in] x The current value of x.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of x.
	 * \endinternal
	 */
	static double int_dp_dgammaL(double x, void *params);

	/**
	 * \internal
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$\gamma_\mathrm{R}\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f{eqnarray}
	 * \int\limits_{(2-g-2\sqrt{1-g})/g}^{(2-g+2\sqrt{1-g})/g} \mathrm{d}x && \frac{x(1+x^2)^{-5/2}}{\sqrt{4x-g(1+x)^2}} \exp\left[ -\frac{(x\gamma_\mathrm{L} - x\gamma_\mathrm{R})^2}{2(1+x^2)} - \frac{r^2(\gamma_\mathrm{L}^2 + \gamma_\mathrm{R}^2) [4x-g(1+x)^2]}{8g(1+x^2)} \right] \nonumber \\
	 * && \times \left[ 3x\gamma_\mathrm{L} + (2x^2-1)\gamma_\mathrm{L} + \frac{(\gamma_\mathrm{L}+x\gamma_\mathrm{R})^2(x\gamma_\mathrm{L}-\gamma_\mathrm{R})}{1+x^2} - \left( 1 + \frac{(\gamma_\mathrm{L} + x\gamma_\mathrm{R})^2}{1+x^2} \right) \frac{r^2\gamma_\mathrm{R}[4x-g(1+x)^2]}{4g} \right] \nonumber
	 * \f}
	 *
	 * \param[in] x The current value of x.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of x.
	 * \endinternal
	 */
	static double int_dp_dgammaR(double x, void *params);

	/**
	 * \internal
	 * \brief GSL integrand function for the derivative of the fit function
	 *    integral with respect to \f$r\f$.
	 *
	 * Function to be used in the GSL QAGS routine to evaluate \f[
	 * \int\limits_{(2-g-2\sqrt{1-g})/g}^{(2-g+2\sqrt{1-g})/g} \mathrm{d}x \frac{x\sqrt{4x-g(1+x)^2}}{(1+x^2)^{5/2}} \left( 1 + \frac{(\gamma_\mathrm{L} + x\gamma_\mathrm{R})^2}{1+x^2} \right) \exp\left[ -\frac{(x\gamma_\mathrm{L} - x\gamma_\mathrm{R})^2}{2(1+x^2)} - \frac{r^2(\gamma_\mathrm{L}^2 + \gamma_\mathrm{R}^2) [4x-g(1+x)^2]}{8g(1+x^2)} \right].
	 * \f]
	 *
	 * \param[in] x The current value of x.
	 * \param[in] params The fit parameters, assumed to be in vector<double>
	 *    form (although passed as void* to satisfy GSL requirements).
	 * \return The integrand evaluated at this value of x.
	 * \endinternal
	 */
	static double int_dp_dr(double x, void *params);

public:
	/**
	 * \brief Index for the \f$\gamma_\mathrm{L}\f$ fitting parameter.
	 */
	const static int GAMMAL = 0;

	/**
	 * \brief Index for the \f$\gamma_\mathrm{R}\f$ fitting parameter.
	 */
	const static int GAMMAR = 1;

	/**
	 * \brief Index for the \f$r\f$ fitting parameter.
	 */
	const static int R = 2;

	/**
	 * \brief Index for the norm fitting parameter.
	 */
	const static int NORM = 3;

	AsymmetricResonantFitModel() = delete;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	AsymmetricResonantFitModel(
		const std::list<std::pair<std::array<double, 1>, double>> &data);

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
	virtual ~AsymmetricResonantFitModel() = default;

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
	 * Each element of the Jacobian depends on a particular integral, as
	 * follows.
	 * \f{eqnarray}
	 * \frac{\partial \hat{P}(g)}{\partial \gamma_\mathrm{L}} & = & \frac{N}{g^{3/2}} \mathrm{int\_dp\_dgammaL}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial \gamma_\mathrm{R}} & = & \frac{N}{g^{3/2}} \mathrm{int\_dp\_dgammaR}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial r} & = & \frac{-Nr(\gamma_\mathrm{L}^2+\gamma_\mathrm{R}^2)}{4g^{5/2}} \mathrm{int\_dp\_dr}(g), \\
	 * \frac{\partial \hat{P}(g)}{\partial N} & = & \frac{1}{g^{3/2}} \mathrm{int\_p}(g).
	 * \f}
	 * where int_p is evaluated with AsymmetricResonantFitModel::int_p, and
	 * likewise for AsymmetricResonantFitModel::int_dp_dgammaL,
	 * AsymmetricResonantFitModel::int_dp_dgammaR, and
	 * AsymmetricResonantFitModel::int_dp_dr.
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return The Jacobian evaluated at these independent variables and
	 *    fitting parameters.
	 */
	virtual std::vector<double> jacobian(const std::vector<double> &fitparam,
		const std::array<double, 1> &x, const double f) const;

	/**
	 * \brief Calculates both the residual and Jacobian of the fit for a given
	 *    set of fitting parameters.
	 *
	 * The residual and the derivative with respect to \f$N\f$ both require
	 * the integral in AsymmetricResonantFitModel::int_p.
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return A pair of the residual and Jacobian evaluated at these
	 *    independent variables and fitting parameters.
	 */
	virtual std::pair<double, std::vector<double>> resid_j(
		const std::vector<double> &fitparam, const std::array<double, 1> &x,
		const double f) const;

	/**
	 * \brief Appends default initial guesses to a list.
	 *
	 * \param[in,out] guess A list of initial guesses.
	 */
	virtual void append_default_guesses(std::list<std::vector<double>> &guess)
		const;

	/**
	 * \brief Prints the fit variables from a gsl_vector.
	 *
	 * \param[in] f The output stream.
	 * \param[in] fitparam The fitting parameters.
	 */
	virtual void print_fit(FILE *f, const std::vector<double> &fitparam)
		const;

	/**
	 * \brief Perform post-processing on a set of fit parameters.
	 *
	 * This basic implementation does nothing.
	 *
	 * \param[in,out] fitparams The fitting parameters.
	 */
	virtual void process_fit_parameters(std::vector<double> &fitparams) const;
};

#endif
