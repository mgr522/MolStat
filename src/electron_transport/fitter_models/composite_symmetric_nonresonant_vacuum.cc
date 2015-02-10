/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file composite_symmetric_nonresonant_vacuum.cc
 * \brief Implementation of the fitting model for nonresonant tunneling
 *        (symmetric coupling) combined with background \"vacuum\" tunneling.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include "composite_symmetric_nonresonant_vacuum.h"
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

using namespace std;

namespace molstat {
namespace transport {

std::vector<double> CompositeSymmetricNonresonantVacuumFitModel
	::create_initial_guess(const std::map<std::string, double> &values) const
{
	vector<double> ret(4);

	try
	{
		ret[C] = values.at("c");
		ret[D] = values.at("d");
		ret[GMINUS] = values.at("gminus");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the CompositeSymmetric" \
			"NonresonantVacuumFitModel must specify \"c\", \"d\", and " \
			"\"gminus\" parameters.");
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

CompositeSymmetricNonresonantVacuumFitModel
	::CompositeSymmetricNonresonantVacuumFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(4, data), w(gsl_integration_workspace_alloc(nquad),
		&gsl_integration_workspace_free) 
{}

double CompositeSymmetricNonresonantVacuumFitModel::resid(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1);
	double error, integral;
	gsl_function func;
	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_p;
	func.params = &params;

	const double g = x[0];

	params[C] = fitparam[C];
	params[D] = fitparam[D];
	params[GMINUS] = fitparam[GMINUS];
	params[NORM] = fitparam[NORM];
	params[nfit] = g; // need to pass in the conductance value

	// calculate the integral
	gsl_integration_qags(&func, 0., g, 0., 1.0e-7, nquad, w.get(), &integral,
		&error);

	return integral * fitparam[NORM] - f;
}

std::vector<double> CompositeSymmetricNonresonantVacuumFitModel::jacobian(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1), ret(nfit);
	double error, integral, intc, intd, intgminus;
	gsl_function func;
	func.params = &params;

	const double g = x[0];
	const double c = fitparam[C];
	const double d = fitparam[D];
	const double gminus = fitparam[GMINUS];
	const double norm = fitparam[NORM];

	params[C] = c;
	params[D] = d;
	params[GMINUS] = gminus;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	// evaluate the integrals
	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_dp_dc;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intc, &error);

	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_dp_dd;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intd, &error);

	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the derivatives
	ret[C] = -norm * intc;

	ret[D] = norm * intd;

	ret[GMINUS] = -2.*norm / (k * gminus * gminus * M_SQRTPI) * intgminus;

	ret[NORM] = integral;

	return ret;
}

std::pair<double, std::vector<double>>
	CompositeSymmetricNonresonantVacuumFitModel::resid_j(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	pair<double, vector<double>> ret;
	vector<double> params(nfit+1);
	double error, integral, intc, intd, intgminus;
	gsl_function func;
	func.params = &params;

	ret.second.resize(nfit);

	const double g = x[0];
	const double c = fitparam[C];
	const double d = fitparam[D];
	const double gminus = fitparam[GMINUS];
	const double norm = fitparam[NORM];

	params[C] = c;
	params[D] = d;
	params[GMINUS] = gminus;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	// evaluate the integrals
	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_dp_dc;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intc, &error);

	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_dp_dd;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intd, &error);

	func.function = &CompositeSymmetricNonresonantVacuumFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the residual and derivatives
	ret.first = norm * integral - f;

	ret.second[C] = -norm * intc;

	ret.second[D] = norm * intd;

	ret.second[GMINUS] = -2.*norm / (k * gminus * gminus * M_SQRTPI) * intgminus;

	ret.second[NORM] = integral;

	return ret;
}

void CompositeSymmetricNonresonantVacuumFitModel::append_default_guesses(
	std::list<std::vector<double>> &guess) const
{
	const list<double> list_c{40., 50., 75., 100., 150.},
		list_d{4., 5., 6.},
		list_gminus{1.e-7, 1.e-6, 1.e-5};

	for(list<double>::const_iterator c = list_c.cbegin();
		c != list_c.cend(); ++c)
	{
	for(list<double>::const_iterator d = list_d.cbegin();
		d != list_d.cend(); ++d)
	{
	for(list<double>::const_iterator gminus = list_gminus.cbegin();
		gminus != list_gminus.cend(); ++gminus)
	{

		vector<double> init(nfit);
		init[C] = *c;
		init[D] = *d;
		init[GMINUS] = *gminus;
		init[NORM] = 1.;

		guess.emplace_back(init);
	}}}
}

void CompositeSymmetricNonresonantVacuumFitModel::print_fit(std::ostream &out,
	const std::vector<double> &fitparam) const
{
	out << "c=" << scientific << setprecision(4) << fitparam[C] <<
		", d=" << scientific << setprecision(4) << fitparam[D] <<
		", gminus=" << scientific << setprecision(4) << fitparam[GMINUS] <<
		", norm=" << scientific << setprecision(4) << fitparam[NORM];
}

void CompositeSymmetricNonresonantVacuumFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const
{
	// sometimes both c and d go negative
	if(fitparams[C] < 0. && fitparams[D] < 0.)
	{
		fitparams[C] = -fitparams[C];
		fitparams[D] = -fitparams[D];
	}
}

double CompositeSymmetricNonresonantVacuumFitModel::int_p(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double c = fitparams[C];
	const double d = fitparams[D];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = c*sqrt(temp1) - d*sqrt(temp2);

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus)))
		* exp(-0.5 * temp3 * temp3 / temp2)
		/ (gp * sqrt(temp1 * temp2 * temp2 * temp2));
}

double CompositeSymmetricNonresonantVacuumFitModel::int_dp_dc(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double c = fitparams[C];
	const double d = fitparams[D];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = c*sqrt(temp1) - d*sqrt(temp2);

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus))) * temp3
		* exp(-0.5 * temp3 * temp3 / temp2)
		/ (gp * sqrt(temp2) * temp2 * temp2);
}

double CompositeSymmetricNonresonantVacuumFitModel::int_dp_dd(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double c = fitparams[C];
	const double d = fitparams[D];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = c*sqrt(temp1) - d*sqrt(temp2);

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus))) * temp3
		* exp(-0.5 * temp3 * temp3 / temp2)
		/ (gp * sqrt(temp1) * temp2 * temp2);
}

double CompositeSymmetricNonresonantVacuumFitModel::int_dp_dgminus(
	double gp, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[6];
	const double c = fitparams[C];
	const double d = fitparams[D];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = c*sqrt(temp1) - d*sqrt(temp2);
	const double temp4 = (gp - gminus) / (k * gminus);

	return exp(-temp4*temp4 -0.5 * temp3 * temp3 / temp2)
		/ sqrt(temp1 * temp2 * temp2 * temp2);
}

} // namespace molstat::transport
} // namespace molstat
