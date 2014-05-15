/**
 * \file fit_model.h
 * \brief Defines an abstract class encapsulating a model for nonlinear
 *    fitting.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __fit_model_h__
#define __fit_model_h__

#include <memory>
#include <vector>
#include <array>
#include <utility>
#include <list>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>

/**
 * \brief Abstract class encapsulating models that data be fit to.
 *
 * Throughout, the function we fit to will be denoted \f$f\f$. For
 * extensibility, this class is templated over the number of independent
 * variables of \f$f\f$ (note that this is not the number of fitting
 * parameters). In other words, \f$f:\mathbb{R}^n\to\mathbb{R}\f$. The fit
 * function has nfit fitting parameters.
 *
 * Derived classes do not have to deal with the GSL particulars. Rather, they
 * only need to implement functions to evaluate the fit function and/or its
 * Jacobian for specific values of the independent variables and fitting
 * parameters.
 *
 * \tparam N The number of independent variables in the fit function.
 */
template<std::size_t N>
class FitModel {
private:
	/**
	 * \brief The data we fit against.
	 */
	const std::vector<std::pair<std::array<double, N>, double>> &data;

protected:
	/**
	 * \brief The number of fitting parameters in the model.
	 */
	const std::size_t nfit;

public:
	FitModel() = delete;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] nfit_ The number of fitting parameters in the model.
	 * \param[in] data_ The data, organized in pairs of array<double, N>,
	 *    double objects.
	 */
	FitModel(const std::size_t nfit_,
		const std::vector<std::pair<std::array<double, N>, double>> &data_);

	/**
	 * \brief Destructor.
	 */
	virtual ~FitModel() = default;

	/**
	 * \brief Calculates the residuals of the fit for each set of model
	 *    parameters and for a given set of fitting parameters.
	 *
	 * Conforms to the GSL functional form for gsl_multifit_function_f.
	 *
	 * \param[in] x Vector of fit variables.
	 * \param[in] params Model-dependent parameters.
	 * \param[out] f Vector of residuals for each data point.
	 * \return GSL_SUCCESS for success, otherwise a GSL error code.
	 */
	int f(const gsl_vector *x, void *params, gsl_vector *f) const;

	/**
	 * \brief Calculates the Jacobian of the fit function for a given set of
	 *    fitting parameters.
	 *
	 * Conforms to the GSL functional form for gsl_multifit_function_df.
	 *
	 * \param[in] x Vector of fit variables.
	 * \param[in] params Model-dependent parameters.
	 * \param[out] J Jacobian matrix.
	 * \return 0 for success, otherwise a GSL error code.
	 */
	int df(const gsl_vector *x, void *params, gsl_matrix *J) const;

	/**
	 * \brief Calculates both the residuals and Jacobian of the fit for a given
	 *    set of fitting parameters.
	 *
	 * Conforms to the GSL functional form for gsl_multifit_function_fdf.
	 *
	 * \param[in] x Vector of fit variables.
	 * \param[in] params Model-dependent parameters.
	 * \param[out] f Vector of residuals for each data point.
	 * \param[out] J Jacobian matrix.
	 * \return GSL_SUCCESS for success, otherwise a GSL error code.
	 */
	int fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
		const;

	/**
	 * \brief Evaluates the residual for this data point at a given set of
	 *    independent variables and fitting parameters.
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return The function evaluated at these independent variables and
	 *    fitting parameters.
	 */
	virtual double resid(const std::vector<double> &fitparam,
		const std::array<double, N> &x, const double f) const = 0;

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
		const std::array<double, N> &x, const double f) const = 0;

	/**
	 * \brief Calculates both the residuals and Jacobian of the fit for a given
	 *    set of fitting parameters.
	 *
	 * This basic implementation simply calls the two independent routines and
	 * is provided for convenience. Depending on the fit function, it may be
	 * more efficient to override this function and calculate them both
	 * together.
	 *
	 * \param[in] fitparam The fitting parameters.
	 * \param[in] x The independent variables for the function.
	 * \param[in] f The observed value of the fit function at x.
	 * \return A pair of the residual and Jacobian evaluated at these
	 *    independent variables and fitting parameters.
	 */
	virtual std::pair<double, std::vector<double>> resid_j(
		const std::vector<double> &fitparam, const std::array<double, N> &x,
		const double f) const;

	/**
	 * \brief Gets a GSL handle for fitting data to this model.
	 *
	 * \return The GSL nonlinear fitting handle.
	 */
	gsl_multifit_function_fdf gsl_handle() const;

	/**
	 * \brief Returns a list of initial guesses to use when fitting the data.
	 *
	 * \return A list of initial guesses.
	 */
	virtual std::list<std::vector<double>> initial_guesses() const = 0;
};

// Other function prototypes
/**
 * \brief Converts a gsl_vector to a std::vector<double>.
 *
 * \param[in] gslv The gsl_vector.
 * \return The std::vector<double>.
 */
std::vector<double> gsl_to_std(const gsl_vector *gslv);

// Implementation of templated class functions
template<std::size_t N>
FitModel<N>::FitModel(const std::size_t nfit_,
	const std::vector<std::pair<std::array<double, N>, double>> &data_)
	: data(data_), nfit(nfit_) {
}

template<std::size_t N>
int FitModel<N>::f(const gsl_vector *x, void *params, gsl_vector *f) const {
	std::size_t i;
	const std::size_t n = data.size();
	const std::vector<double> fitparam(gsl_to_std(x));

	for(i = 0; i < n; ++i) {
		// get the data point
		const std::pair<std::array<double, N>, double> &datai = data[i];

		gsl_vector_set(f, i, resid(fitparam, datai.first, datai.second));
	}

	return GSL_SUCCESS;
}

template<std::size_t N>
int FitModel<N>::df(const gsl_vector *x, void *params, gsl_matrix *J) const {
	std::size_t i, j;
	const std::size_t n = data.size();
	const std::vector<double> fitparam(gsl_to_std(x));

	for(i = 0; i < n; ++i) {
		// get the data point
		const std::pair<std::array<double, N>, double> &datai = data[i];

		const std::vector<double> jac(jacobian(fitparam, datai.first,
			datai.second));

		// set the matrix elements
		for(j = 0; j < nfit; ++j)
			gsl_matrix_set(J, i, j, jac[j]);
	}

	return GSL_SUCCESS;
}

template<std::size_t N>
int FitModel<N>::fdf(const gsl_vector *x, void *params, gsl_vector *f,
	gsl_matrix *J) const {

	std::size_t i, j;
	const std::size_t n = data.size();
	const std::vector<double> fitparam(gsl_to_std(x));

	for(i = 0; i < n; ++i) {
		// get the data point
		const std::pair<std::array<double, N>, double> &datai = data[i];

		const std::pair<double, std::vector<double>> 
			vals(resid_j(fitparam, datai.first, datai.second));

		// set the residual
		gsl_vector_set(f, i, vals.first);

		// set the matrix elements
		for(j = 0; j < nfit; ++j)
			gsl_matrix_set(J, i, j, vals.second[j]);
	}

	return GSL_SUCCESS;
}

template<std::size_t N>
std::pair<double, std::vector<double>> FitModel<N>::resid_j(
	const std::vector<double> &fitparam, const std::array<double, N> &x,
	const double f) const {

	std::pair<double, std::vector<double>> ret;

	ret.first = resid(fitparam, x, f);
	ret.second = jacobian(fitparam, x, f);

	return ret;
}

template<std::size_t N>
gsl_multifit_function_fdf FitModel<N>::gsl_handle() const {
	gsl_multifit_function_fdf fit;

	fit.f = &this->f;
	fit.df = &this->df;
	fit.fdf = &this->fdf;
	fit.n = data.size();
	fit.p = nfit;
	fit.params = nullptr;

	return fit;
}

#endif
