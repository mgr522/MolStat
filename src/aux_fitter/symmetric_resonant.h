/**
 * \file symmetric_resonant.h
 * \brief Fitting model for resonant tunneling (symmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __symmetric_resonant_h__
#define __symmetric_resonant_h__

#include "fit_model.h"

/**
 * \brief Abstract class encapsulating models that data be fit to.
 *
 */
class SymmetricResonantFitModel : public FitModel<1> {
public:
	/**
	 * \brief Index for the gamma fitting parameter.
	 */
	const static int GAMMA = 0;

	/**
	 * \brief Index for the norm fitting parameter.
	 */
	const static int NORM = 1;

	SymmetricResonantFitModel() = delete;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] data The data, organized in pairs of array<double, 1>,
	 *    double objects.
	 */
	SymmetricResonantFitModel(
		const std::list<std::pair<std::array<double, 1>, double>> &data);

	/**
	 * \brief Destructor.
	 */
	virtual ~SymmetricResonantFitModel() = default;

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
	virtual void print_fit(FILE *f, const std::shared_ptr<gsl_vector> fitparam)
		const;
};

#endif
