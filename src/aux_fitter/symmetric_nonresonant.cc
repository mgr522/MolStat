/**
 * \file symmetric_nonresonant.cc
 * \brief Implementation of the fitting model for nonresonant tunneling
 *        (symmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "symmetric_nonresonant.h"

using namespace std;

SymmetricNonresonantFitModel::SymmetricNonresonantFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(3, data) {
}

double SymmetricNonresonantFitModel::resid(const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const {

	// get the current parameters and independent variable
	const double g = x[0];
	const double c = fitparam[C];
	const double d = fitparam[D];
	const double norm = fitparam[NORM];

	const double cd = c*sqrt(g) - d*sqrt(1. - g);
	const double expcd = exp(-0.5*cd*cd / (1. - g));

	const double model = norm / sqrt(g*(1.-g)*(1.-g)*(1.-g)) * expcd;

	return model - f;
}

std::vector<double> SymmetricNonresonantFitModel::jacobian(
	const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const {

	// get the current parameters and independent variable
	const double g = x[0];
	const double c = fitparam[C];
	const double d = fitparam[D];
	const double norm = fitparam[NORM];

	const double cd = c*sqrt(g) - d*sqrt(1. - g);
	const double expcd = exp(-0.5*cd*cd / (1. - g));

	vector<double> ret(nfit);

	ret[C] = -norm * cd * expcd / ((1.-g)*(1.-g)*sqrt(1.-g));

	ret[D] = norm * cd * expcd / ((1.-g)*(1.-g)*sqrt(g));

	ret[NORM] = expcd / ((1.-g)*sqrt(g*(1.-g)));

	return ret;
}

std::list<std::vector<double>> SymmetricNonresonantFitModel::initial_guesses()
	const {

	const list<double> list_c{50., 100., 200., 300., 400., 500.},
		list_d{5., 10., 20., 30., 40., 50.};
	list<vector<double>> ret;

	for(list<double>::const_iterator c = list_c.cbegin();
		c != list_c.cend(); ++c) {
	for(list<double>::const_iterator d = list_d.cbegin();
		d != list_d.cend(); ++d) {

		vector<double> init(nfit);
		init[C] = *c;
		init[D] = *d;
		init[NORM] = 1.;

		ret.emplace_back(init);
	}}

	return ret;
}

void SymmetricNonresonantFitModel::print_fit(FILE *f,
	const std::vector<double> &fitparam) const {

	fprintf(f, "c=% .3e, d=% .3e, norm=% .3e", fitparam[C], fitparam[D],
		fitparam[NORM]);
}