/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file composite_interference_vacuum.cc
 * \brief Implementation of the fitting model for transport around an
 *        interference feature combined with background \"vacuum\" tunneling.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#include "composite_interference_vacuum.h"
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

using namespace std;

namespace molstat {
namespace transport {

std::vector<double> CompositeInterferenceVacuumFitModel
	::create_initial_guess(const std::map<std::string, double> &values) const
{
	vector<double> ret(3);

	try
	{
		ret[F] = values.at("f");
		ret[GMINUS] = values.at("gminus");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the CompositeInterference" \
			"VacuumFitModel must specify \"f\" and \"gminus\" parameters.");
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

CompositeInterferenceVacuumFitModel
	::CompositeInterferenceVacuumFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(3, data), w(gsl_integration_workspace_alloc(nquad),
		&gsl_integration_workspace_free) 
{}

double CompositeInterferenceVacuumFitModel::resid(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1);
	double error, integral;
	gsl_function func;
	func.function = &CompositeInterferenceVacuumFitModel::int_p;
	func.params = &params;

	const double g = x[0];

	params[F] = fitparam[F];
	params[GMINUS] = fitparam[GMINUS];
	params[NORM] = fitparam[NORM];
	params[nfit] = g; // need to pass in the conductance value

	// calculate the integral
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad,
		w.get(), &integral, &error);

	const double model = integral * fitparam[NORM];

	// owing to the singularity in the form -- the data can span several
	// orders of magnitude with most points much smaller than a few --
	// scale by the size of the point to give more weight to the smaller
	// points
	return (model - f) / f;
}

std::vector<double> CompositeInterferenceVacuumFitModel::jacobian(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1), ret(nfit);
	double error, integral, intf, intgminus;
	gsl_function func;
	func.params = &params;

	const double g = x[0];
	const double fparam = fitparam[F];
	const double gminus = fitparam[GMINUS];
	const double norm = fitparam[NORM];

	params[F] = fparam;
	params[GMINUS] = gminus;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	// evaluate the integrals
	func.function = &CompositeInterferenceVacuumFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &CompositeInterferenceVacuumFitModel::int_dp_df;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intf, &error);

	func.function = &CompositeInterferenceVacuumFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the derivatives
	ret[F] = -fparam * norm * intf;
	ret[F] /= f; // scaling as described above

	ret[GMINUS] = -2.*norm / (k * gminus * gminus * M_SQRTPI) * intgminus;
	ret[GMINUS] /= f; // scaling, again

	ret[NORM] = integral;
	ret[NORM] /= f; // scaling, again

	return ret;
}

std::pair<double, std::vector<double>>
	CompositeInterferenceVacuumFitModel::resid_j(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	pair<double, vector<double>> ret;
	vector<double> params(nfit+1);
	double error, integral, intf, intgminus;
	gsl_function func;
	func.params = &params;

	ret.second.resize(nfit);

	const double g = x[0];
	const double fparam = fitparam[F];
	const double gminus = fitparam[GMINUS];
	const double norm = fitparam[NORM];

	params[F] = fparam;
	params[GMINUS] = gminus;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	// evaluate the integrals
	func.function = &CompositeInterferenceVacuumFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &CompositeInterferenceVacuumFitModel::int_dp_df;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intf, &error);

	func.function = &CompositeInterferenceVacuumFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the residual and derivatives
	ret.first = norm * integral - f;
	ret.first /= f; // scaling as described above

	ret.second[F] = -fparam * norm * intf;
	ret.second[F] /= f; // scaling, again

	ret.second[GMINUS] = -2.*norm / (k * gminus * gminus * M_SQRTPI)
		* intgminus;
	ret.second[GMINUS] /= f; // scaling, again

	ret.second[NORM] = integral;
	ret.second[NORM] /= f; // scaling, again

	return ret;
}

void CompositeInterferenceVacuumFitModel::append_default_guesses(
	std::list<std::vector<double>> &guess) const
{
	const list<double> list_f{1., 10., 100.},
		list_gminus{1.e-7, 1.e-6, 1.e-5};

	for(list<double>::const_iterator f = list_f.cbegin();
		f != list_f.cend(); ++f)
	{
	for(list<double>::const_iterator gminus = list_gminus.cbegin();
		gminus != list_gminus.cend(); ++gminus)
	{

		vector<double> init(nfit);
		init[F] = *f;
		init[GMINUS] = *gminus;
		init[NORM] = 1.;

		guess.emplace_back(init);
	}}
}

void CompositeInterferenceVacuumFitModel::print_fit(std::ostream &out,
	const std::vector<double> &fitparam) const
{
	out << "f=" << scientific << setprecision(4) << fitparam[F] <<
		", gminus=" << scientific << setprecision(4) << fitparam[GMINUS] <<
		", norm=" << scientific << setprecision(4) << fitparam[NORM];
}

void CompositeInterferenceVacuumFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const
{
	// make sure `f` is positive
	if(fitparams[F] < 0.)
		fitparams[F] = -fitparams[F];
}

double CompositeInterferenceVacuumFitModel::int_p(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[3];
	const double f = fitparams[F];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus)))
		* exp(-0.5 * f * f * temp1) / (gp * sqrt(temp1));
}

double CompositeInterferenceVacuumFitModel::int_dp_df(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[3];
	const double f = fitparams[F];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus))) * sqrt(temp1)
		* exp(-0.5 * f * f * temp1) / gp;
}

double CompositeInterferenceVacuumFitModel::int_dp_dgminus(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[3];
	const double f = fitparams[F];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = (gp - gminus) / (k * gminus);

	return exp(-temp2*temp2 - 0.5 * f * f * temp1) / sqrt(temp1);
}

} // namespace molstat::transport
} // namespace molstat
