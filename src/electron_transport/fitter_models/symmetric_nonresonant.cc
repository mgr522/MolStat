/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file symmetric_nonresonant.cc
 * \brief Implementation of the fitting model for nonresonant tunneling
 *        (symmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "symmetric_nonresonant.h"
#include <iomanip>

using namespace std;

namespace molstat {
namespace transport {

std::vector<double> SymmetricNonresonantFitModel::create_initial_guess(
	const std::map<std::string, double> &values) const
{
	vector<double> ret(3);

	try
	{
		ret[C] = values.at("c");
		ret[D] = values.at("d");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the Symmetric" \
			"NonresonantFitModel must specify \"c\" and \"d\" parameters.");
	}

	try
	{
		ret[NORM] = values.at("norm");
	}
	catch(const out_of_range &e)
	{
		// norm not specified
		ret[NORM] = 1.;
	}

	return ret;
}

SymmetricNonresonantFitModel::SymmetricNonresonantFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(3, data)
{
}

double SymmetricNonresonantFitModel::resid(const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const
{
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
	const std::array<double, 1> &x, const double f) const
{
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

void SymmetricNonresonantFitModel::append_default_guesses(
	std::list<std::vector<double>> &guess) const
{
	const list<double> list_c{50., 100., 200., 300., 400., 500.},
		list_d{5., 10., 20., 30., 40., 50.};

	for(list<double>::const_iterator c = list_c.cbegin();
		c != list_c.cend(); ++c)
	{
	for(list<double>::const_iterator d = list_d.cbegin();
		d != list_d.cend(); ++d)
	{

		vector<double> init(nfit);
		init[C] = *c;
		init[D] = *d;
		init[NORM] = 1.;

		guess.emplace_back(init);
	}}
}

void SymmetricNonresonantFitModel::print_fit(std::ostream &out,
	const std::vector<double> &fitparam) const
{
	out << "c=" << scientific << setprecision(4) << fitparam[C] << ", d=" <<
		scientific << setprecision(4) << fitparam[D] << ", norm=" << scientific
		<< setprecision(4) << fitparam[NORM];
}

void SymmetricNonresonantFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const
{
	// sometimes both c and d go negative
	if(fitparams[C] < 0. && fitparams[D] < 0.)
	{
		fitparams[C] = -fitparams[C];
		fitparams[D] = -fitparams[D];
	}
}

} // namespace molstat::transport
} // namespace molstat
