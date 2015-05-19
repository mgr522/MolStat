/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file interference.cc
 * \brief Implementation of the fitting model for transport around an
 *    interference feature.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#include "interference.h"
#include <iomanip>

using namespace std;

namespace molstat {
namespace transport {

std::vector<double> InterferenceFitModel::create_initial_guess(
	const std::map<std::string, double> &values) const
{
	vector<double> ret(2);

	try
	{
		ret[F] = values.at("f");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the Interference" \
			"FitModel must specify the \"f\" parameter.");
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

InterferenceFitModel::InterferenceFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(2, data)
{
}

double InterferenceFitModel::resid(const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const
{
	// get the current fit parameters and independent variable
	const double g = x[0];
	const double fparam = fitparam[F];
	const double norm = fitparam[NORM];
	
	const double model = norm / sqrt(g) * exp(-0.5*fparam*fparam * g);

	// owing to the singularity in the form -- the data can span several
	// orders of magnitude with most points much smaller than a few --
	// scale by the size of the point to give more weight to the smaller
	// points
	return (model - f) / f;
}

std::vector<double> InterferenceFitModel::jacobian(
	const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const
{
	// get the current fit parameters and independent variable
	const double g = x[0];
	const double fparam = fitparam[F];
	const double norm = fitparam[NORM];

	vector<double> ret(nfit);

	ret[F] = -norm * fparam * sqrt(g) * exp(-0.5*fparam*fparam*g);
	ret[F] /= f; // scaling as described above

	ret[NORM] = exp(-0.5*fparam*fparam*g) / sqrt(g);
	ret[NORM] /= f; // scaling, again

	return ret;
}

void InterferenceFitModel::append_default_guesses(
	std::list<std::vector<double>> &guess) const
{
	const list<double> list_f{1., 10., 100.};

	for(list<double>::const_iterator f = list_f.cbegin();
		f != list_f.cend(); ++f)
	{

		vector<double> init(nfit);
		init[F] = *f;
		init[NORM] = 1.;

		guess.emplace_back(init);
	}
}

void InterferenceFitModel::print_fit(std::ostream &out,
	const std::vector<double> &fitparam) const
{
	out << "f=" << scientific << setprecision(4) << fitparam[F] <<
		", norm=" << scientific << setprecision(4) << fitparam[NORM];
}

void InterferenceFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const
{
	// sometimes f goes negative.
	if(fitparams[F] < 0.)
		fitparams[F] = -fitparams[F];
}

} // namespace molstat::transport
} // namespace molstat
