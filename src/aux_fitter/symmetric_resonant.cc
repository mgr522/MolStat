/**
 * \file symmetric_resonant.cc
 * \brief Implementation of the fitting model for resonant tunneling
 *        (symmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "symmetric_resonant.h"

using namespace std;

SymmetricResonantFitModel::SymmetricResonantFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(2, data) {
}

double SymmetricResonantFitModel::resid(const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const {

	// get the current fit parameters and independent variable
	const double g = x[0];
	const double gamma = fitparam[GAMMA];
	const double norm = fitparam[NORM];
	
	const double model = norm / sqrt(g*g*g*(1.0 - g))
		* exp(-0.5*gamma*gamma*(1.0 - g) / g);

	// owing to the singularity in the form -- the data can span several
	// orders of magnitude with most points much smaller than a few --
	// scale by the size of the point to give more weight to the smaller
	// points
	return (model - f) / f;
}

std::vector<double> SymmetricResonantFitModel::jacobian(
	const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const {

	// get the current fit parameters and independent variable
	const double g = x[0];
	const double gamma = fitparam[GAMMA];
	const double norm = fitparam[NORM];

	vector<double> ret(nfit);

	ret[GAMMA] = -gamma * norm * sqrt((1.0 - g) / g)
		* exp(-0.5*gamma*gamma*(1.0-g)/g) / (g*g);
	ret[GAMMA] /= f; // scaling as described above

	ret[NORM] = 1.0 / sqrt(g*g*g*(1.0 - g))
		* exp(-0.5*gamma*gamma*(1.0 - g) / g);
	ret[NORM] /= f; // scaling, again

	return ret;
}

std::list<std::vector<double>> SymmetricResonantFitModel::initial_guesses()
	const {

	const list<double> list_gamma{5., 10., 20., 35., 50.};
	list<vector<double>> ret;

	for(list<double>::const_iterator gamma = list_gamma.cbegin();
		gamma != list_gamma.cend(); ++gamma) {

		vector<double> init(nfit);
		init[GAMMA] = *gamma;
		init[NORM] = 1.;

		ret.emplace_back(init);
	}

	return ret;
}

void SymmetricResonantFitModel::print_fit(FILE *f,
	const std::vector<double> &fitparam) const {

	fprintf(f, "gamma=% .3e, norm=% .3e", fitparam[GAMMA], fitparam[NORM]);
}

void SymmetricResonantFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const {

	// sometimes gamma goes negative.
	if(fitparams[GAMMA] < 0.)
		fitparams[GAMMA] = -fitparams[GAMMA];
}
