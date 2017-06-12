/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

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
		ret[CEPSILON] = values.at("cepsilon");
		ret[CGAMMA] = values.at("cgamma");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the Symmetric" \
			"NonresonantFitModel must specify \"cepsilon\" and \"cgamma\" " \
			"parameters.");
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
	const double ceps = fitparam[CEPSILON];
	const double cgamma = fitparam[CGAMMA];
	const double norm = fitparam[NORM];

	const double cd = ceps*sqrt(g) - cgamma*sqrt(1. - g);
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
	const double ceps = fitparam[CEPSILON];
	const double cgamma = fitparam[CGAMMA];
	const double norm = fitparam[NORM];

	const double cd = ceps*sqrt(g) - cgamma*sqrt(1. - g);
	const double expcd = exp(-0.5*cd*cd / (1. - g));

	vector<double> ret(nfit);

	ret[CEPSILON] = -norm * cd * expcd / ((1.-g)*(1.-g)*sqrt(1.-g));

	ret[CGAMMA] = norm * cd * expcd / ((1.-g)*(1.-g)*sqrt(g));

	ret[NORM] = expcd / ((1.-g)*sqrt(g*(1.-g)));

	return ret;
}

void SymmetricNonresonantFitModel::append_default_guesses(
	std::list<std::vector<double>> &guess) const
{
	const list<double> list_ceps{50., 100., 200., 300., 400., 500.},
		list_cgamma{5., 10., 20., 30., 40., 50.};

	for(list<double>::const_iterator ceps = list_ceps.cbegin();
		ceps != list_ceps.cend(); ++ceps)
	{
	for(list<double>::const_iterator cgamma = list_cgamma.cbegin();
		cgamma != list_cgamma.cend(); ++cgamma)
	{

		vector<double> init(nfit);
		init[CEPSILON] = *ceps;
		init[CGAMMA] = *cgamma;
		init[NORM] = 1.;

		guess.emplace_back(init);
	}}
}

void SymmetricNonresonantFitModel::print_fit(std::ostream &out,
	const std::vector<double> &fitparam) const
{
	out << "cepsilon=" << scientific << setprecision(4) << fitparam[CEPSILON] <<
		", cgamma=" << scientific << setprecision(4) << fitparam[CGAMMA] <<
		", norm=" << scientific << setprecision(4) << fitparam[NORM];
}

void SymmetricNonresonantFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const
{
	// sometimes both cepsilon and cgamma go negative
	if(fitparams[CEPSILON] < 0. && fitparams[CGAMMA] < 0.)
	{
		fitparams[CEPSILON] = -fitparams[CEPSILON];
		fitparams[CGAMMA] = -fitparams[CGAMMA];
	}
}

bool SymmetricNonresonantFitModel::is_good_fit(
	const std::vector<double> &fitparams) const
{
	return (fitparams[CEPSILON] > 0 && fitparams[CGAMMA] > 0.);
}

} // namespace molstat::transport
} // namespace molstat
