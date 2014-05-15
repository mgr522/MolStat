/**
 * \file symmetric_resonant.cc
 * \brief Implementation of the fitting model for resonant tunneling
 *        (symmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "symmetric_resonant.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

using namespace std;

SymmetricResonantFitModel::SymmetricResonantFitModel(
	const std::vector<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(2, data) {
}

double SymmetricResonantFitModel::fit_function(
	const std::vector<double> &fitparam,
	const std::array<double, 1> &x) const {

	// get the current fit parameters and independent variable
	const double g = x[0];
	const double gamma = fitparam[GAMMA];
	const double norm = fitparam[NORM];
	
	return norm / sqrt(g*g*g*(1.0 - g))
		* exp(-0.5*gamma*gamma*(1.0 - g) / g);

		// owing to the singularity in the form -- the data can span several
		// orders of magnitude with most points much smaller than a few --
		// scale by the size of the point to give more weight to the smaller
		// points
		//gsl_vector_set(resid, i, (pdfi - pdf[i]) / pdf[i]);
}

std::vector<double> SymmetricResonantFitModel::jacobian(
	const std::vector<double> &fitparam,
	const std::array<double, 1> &x) const {

	// get the current fit parameters and independent variable
	const double g = x[0];
	const double gamma = fitparam[GAMMA];
	const double norm = fitparam[NORM];

	vector<double> ret(2);

	ret[GAMMA] = -gamma * norm * sqrt((1.0 - g) / g)
		* exp(-0.5*gamma*gamma*(1.0-g)/g) / (g*g);
	//ret[GAMMA] /= pdf[i]; // scaling as described above

	ret[NORM] = 1.0 / sqrt(g*g*g*(1.0 - g))
		* exp(-0.5*gamma*gamma*(1.0 - g) / g);
	//ret[NORM] /= pdf[i]; // scaling, again

	return ret;
}

std::list<std::vector<double>> SymmetricResonantFitModel::initial_guesses()
	const {

	list<vector<double>> ret;

	return ret;
}
