/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file fit_model_interface.h
 * \brief Defines an abstract class encapsulating a model for nonlinear
 *    fitting.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __fit_model_interface_h__
#define __fit_model_interface_h__

#include <memory>
#include <vector>
#include <array>
#include <utility>
#include <list>
#include <string>
#include <map>
#include <functional>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include "../string_tools.h"

namespace molstat {

/**
 * \brief Abstract class encapsulating models can fit data.
 *
 * Throughout, the function we fit to will be denoted \f$f\f$. For
 * extensibility, this class is templated over the number of independent
 * variables of \f$f\f$ (note that this is not the number of fitting
 * parameters). In other words, \f$f:\mathbb{R}^N\to\mathbb{R}\f$. The fit
 * function has \f$N=\f$`nfit` fitting parameters.
 *
 * Derived classes do not have to deal with the GSL particulars. Rather, they
 * only need to implement functions to evaluate the fit function and its
 * Jacobian for specific values of the independent variables and fitting
 * parameters.
 *
 * \tparam N The number of independent variables in the fit function.
 */
template<std::size_t N>
class FitModel {
private:
	/**
	 * \internal
	 * \brief The data we fit against.
	 *
	 * This is a list of pairs; each pair contains the independent variables
	 * (the array) and the observed value of the fit function.
	 * \endinternal
	 *
	 * \todo Rework this class and main-fitter.cc so that the data is moved
	 *    into this class.
	 */
	const std::list<std::pair<std::array<double, N>, double>> &data;

protected:
	/**
	 * \brief Produces an initial guess from a map of names to values.
	 *
	 * When specifying initial guesses, the user puts something like
	 * "guess c c-value d d-value" in the input file.
	 * FitModel::append_initial_guess converts the list of names and values
	 * into a map<string, double> object. This function, which depends on the
	 * model, makes sure all appropriate variables are set and orders them
	 * correctly.
	 *
	 * \param[in] values The map of names to values.
	 * \return A vector containing the initial guess.
	 */
	virtual std::vector<double> create_initial_guess(
		const std::map<std::string, double> &values) const = 0;

public:
	/**
	 * \brief The number of fitting parameters in the model.
	 */
	const std::size_t nfit;

	FitModel() = delete;
	virtual ~FitModel() = default;

	/**
	 * \brief Constructor requiring the data that we will be fitting against.
	 *
	 * \param[in] nfit_ The number of fitting parameters in the model.
	 * \param[in] data_ The data, organized in pairs of array<double, N>,
	 *    double objects.
	 */
	FitModel(const std::size_t nfit_,
		const std::list<std::pair<std::array<double, N>, double>> &data_);

	/**
	 * \internal
	 * \brief Calculates the residuals of the fit for each set of model
	 *    parameters and for a given set of fitting parameters.
	 *
	 * Conforms to the GSL functional form for gsl_multifit_function_f.
	 *
	 * \param[in] x Vector of fit variables.
	 * \param[in] model The FitModel<N> model object.
	 * \param[out] f Vector of residuals for each data point.
	 * \return GSL_SUCCESS for success, otherwise a GSL error code.
	 * \endinternal
	 */
	static int f(const gsl_vector *x, void *model, gsl_vector *f);

	/**
	 * \internal
	 * \brief Calculates the Jacobian of the fit function for a given set of
	 *    fitting parameters.
	 *
	 * Conforms to the GSL functional form for gsl_multifit_function_df.
	 *
	 * \param[in] x Vector of fit variables.
	 * \param[in] model The FitModel<N> model object.
	 * \param[out] J Jacobian matrix.
	 * \return 0 for success, otherwise a GSL error code.
	 * \endinternal
	 */
	static int df(const gsl_vector *x, void *model, gsl_matrix *J);

	/**
	 * \internal
	 * \brief Calculates both the residuals and Jacobian of the fit for a given
	 *    set of fitting parameters.
	 *
	 * Conforms to the GSL functional form for gsl_multifit_function_fdf.
	 *
	 * \param[in] x Vector of fit variables.
	 * \param[in] model The FitModel<N> model object.
	 * \param[out] f Vector of residuals for each data point.
	 * \param[out] J Jacobian matrix.
	 * \return GSL_SUCCESS for success, otherwise a GSL error code.
	 * \endinternal
	 */
	static int fdf(const gsl_vector *x, void *model, gsl_vector *f,
		gsl_matrix *J);

	/**
	 * \brief Evaluates the residual for the specified data point at a given set
	 *    of independent variables and fitting parameters.
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
	 * \brief Calculates both the residual and Jacobian of the fit for a given
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
	 * \internal
	 * \brief Gets a GSL handle for fitting data to this model.
	 *
	 * \return The GSL nonlinear fitting handle.
	 * \endinternal
	 */
	gsl_multifit_function_fdf gsl_handle() const;

	/**
	 * \brief Appends default initial guesses to a list.
	 *
	 * \param[in,out] guess A list of initial guesses.
	 */
	virtual void append_default_guesses(std::list<std::vector<double>> &guess)
		const = 0;

	/**
	 * \brief Appends a user-specified initial guess to a list.
	 *
	 * This general routine first converts the user-specified tokens to a
	 * name-to-value list and then calls the model-dependent
	 * FitModel::create_initial_guess. More details can be found in the
	 * documentation for FitModel::create_initial_guess.
	 *
	 * \param[in] tokens The list of user-specified tokens. Even elements
	 *    (0, 2, etc.) are the names of fit variables; odd elements are the
	 *    specified values for the preceding names.
	 * \param[in,out] guess A list of initial guesses.
	 */
	void append_initial_guess(const std::vector<std::string> &tokens,
		std::list<std::vector<double>> &guess) const;

	/**
	 * \brief Prints the fit variables from a gsl_vector.
	 *
	 * \param[in] f The output stream.
	 * \param[in] fitparam The fitting parameters.
	 */
	virtual void print_fit(FILE *f, const std::vector<double> &fitparam)
		const = 0;

	/**
	 * \brief Perform post-processing on a set of fit parameters.
	 *
	 * Sometimes the fit produces unphysical parameters. For example, the
	 * product of \f$\Gamma_\mathrm{L}\f$ and \f$\Gamma_\mathrm{R}\f$ may be
	 * the dominant quantity in the fit (conductance from non-resonant
	 * tunneling), and the fit produces \f$\Gamma_\mathrm{L} < 0\f$ and
	 * \f$ \Gamma_\mathrm{R} < 0\f$. We can fix that to make them both positive,
	 * as required physically.
	 *
	 * This basic implementation does nothing.
	 *
	 * \param[in,out] fitparams The fitting parameters.
	 */
	virtual void process_fit_parameters(std::vector<double> &fitparams) const;
};

// Other function prototypes
/**
 * \brief Shortcut for a factory that makes a FitModel object.
 *
 * FitModel objects are created by passing in the list of data points. This
 * factory is a function type that produces the FitModel from a list of data
 * points.
 *
 * \tparam N The number of independent variables in the fit function.
 */
template<std::size_t N>
using FitModelFactory =
	std::function< std::unique_ptr< FitModel<N> >
		(const std::list<std::pair<std::array<double, N>, double>> &) >;

/**
 * \brief Creates a FitModelFactory for a particular FitModel.
 *
 * \tparam T The type of FitModel we wish to instantiate.
 * \tparam N The number of independent variables in the fit function.
 * \return A function for instantiating the class from a list of data points.
 */
template<typename T, std::size_t N>
inline FitModelFactory<N> GetFitModelFactory() {
	return []
		(const std::list<std::pair<std::array<double, N>, double>> &data)
		-> std::unique_ptr<FitModel<N>> {

		std::unique_ptr<FitModel<N>> ret(new T(data));

		return ret;
	};
}

/**
 * \internal
 * \brief Converts a gsl_vector to a std::vector<double>.
 *
 * \param[in] gslv The gsl_vector.
 * \return The std::vector<double>.
 * \endinternal
 */
std::vector<double> gsl_to_std(const gsl_vector *gslv);

// Implementation of templated class functions
template<std::size_t N>
FitModel<N>::FitModel(const std::size_t nfit_,
	const std::list<std::pair<std::array<double, N>, double>> &data_)
	: data(data_), nfit(nfit_) {
}

template<std::size_t N>
int FitModel<N>::f(const gsl_vector *x, void *model, gsl_vector *f) {
	const FitModel<N> *fitmodel = (FitModel<N>*)model;
	std::size_t i;
	const std::vector<double> fitparam(gsl_to_std(x));
	typename std::list<std::pair<std::array<double, N>, double>>::const_iterator
		iter;

	i = 0;
	for(iter = fitmodel->data.cbegin(); iter != fitmodel->data.cend(); ++iter) {
		// get the data point
		const std::pair<std::array<double, N>, double> &datai = *iter;

		gsl_vector_set(f, i,
			fitmodel->resid(fitparam, datai.first, datai.second));
		++i;
	}

	return GSL_SUCCESS;
}

template<std::size_t N>
int FitModel<N>::df(const gsl_vector *x, void *model, gsl_matrix *J) {
	const FitModel<N> *fitmodel = (FitModel<N>*)model;
	std::size_t i, j;
	const std::vector<double> fitparam(gsl_to_std(x));
	typename std::list<std::pair<std::array<double, N>, double>>::const_iterator
		iter;

	i = 0;
	for(iter = fitmodel->data.cbegin(); iter != fitmodel->data.cend(); ++iter) {
		// get the data point
		const std::pair<std::array<double, N>, double> &datai = *iter;

		const std::vector<double> jac(
			fitmodel->jacobian(fitparam, datai.first, datai.second));

		// set the matrix elements
		for(j = 0; j < fitmodel->nfit; ++j)
			gsl_matrix_set(J, i, j, jac[j]);

		++i;
	}

	return GSL_SUCCESS;
}

template<std::size_t N>
int FitModel<N>::fdf(const gsl_vector *x, void *model, gsl_vector *f,
	gsl_matrix *J) {

	const FitModel<N> *fitmodel = (FitModel<N>*)model;
	std::size_t i, j;
	const std::vector<double> fitparam(gsl_to_std(x));
	typename std::list<std::pair<std::array<double, N>, double>>::const_iterator
		iter;

	i = 0;
	for(iter = fitmodel->data.cbegin(); iter != fitmodel->data.cend(); ++iter) {
		// get the data point
		const std::pair<std::array<double, N>, double> &datai = *iter;

		const std::pair<double, std::vector<double>> vals(
			fitmodel->resid_j(fitparam, datai.first, datai.second));

		// set the residual
		gsl_vector_set(f, i, vals.first);

		// set the matrix elements
		for(j = 0; j < fitmodel->nfit; ++j)
			gsl_matrix_set(J, i, j, vals.second[j]);

		++i;
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

	fit.f = &f;
	fit.df = &df;
	fit.fdf = &fdf;
	fit.n = data.size();
	fit.p = nfit;
	fit.params = (void*)this;

	return fit;
}

template<std::size_t N>
void FitModel<N>::append_initial_guess(const std::vector<std::string> &tokens,
	std::list<std::vector<double>> &guess) const {

	std::map<std::string, double> values;
	size_t i, size = tokens.size();

	if(size % 2 == 1)
		--size; // there's a name but no value to follow it; ignore

	for(i = 0; i < size; i += 2) {
		std::string name = tokens[i];
		name = to_lower(name);
		double value;

		try {
			value = stod(tokens[i+1]);
		}
		catch(const std::invalid_argument &e) {
			// error converting to a double; ignore this pair
			continue;
		}

		values[name] = value;
	}

	guess.emplace_back(create_initial_guess(values));
}

template<std::size_t N>
void FitModel<N>::process_fit_parameters(std::vector<double> &fitparams)
	const {
}

} // namespace molstat

#endif
