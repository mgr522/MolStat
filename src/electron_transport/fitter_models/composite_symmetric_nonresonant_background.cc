/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file composite_symmetric_nonresonant_background.cc
 * \brief Implementation of the fitting model for nonresonant tunneling
 *        (symmetric coupling) combined with background tunneling.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include "composite_symmetric_nonresonant_background.h"
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

using namespace std;

namespace molstat {
namespace transport {

std::vector<double> CompositeSymmetricNonresonantBackgroundFitModel
	::create_initial_guess(const std::map<std::string, double> &values) const
{
	vector<double> ret(4);

	try
	{
		ret[CEPSILON] = values.at("cepsilon");
		ret[CGAMMA] = values.at("cgamma");
		ret[GMINUS] = values.at("gminus");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the CompositeSymmetric" \
			"NonresonantBackgroundFitModel must specify \"cepsilon\", \"cgamma\"," \
			" and \"gminus\" parameters.");
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

CompositeSymmetricNonresonantBackgroundFitModel
	::CompositeSymmetricNonresonantBackgroundFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(4, data), w(gsl_integration_workspace_alloc(nquad),
		&gsl_integration_workspace_free) 
{}

double CompositeSymmetricNonresonantBackgroundFitModel::resid(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1);
	double error, integral;
	gsl_function func;
	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_p;
	func.params = &params;

	const double g = x[0];

	params[CEPSILON] = fitparam[CEPSILON];
	params[CGAMMA] = fitparam[CGAMMA];
	params[GMINUS] = fitparam[GMINUS];
	params[NORM] = fitparam[NORM];
	params[nfit] = g; // need to pass in the conductance value

	// calculate the integral
	gsl_integration_qags(&func, 0., g, 0., 1.0e-7, nquad, w.get(), &integral,
		&error);

	return integral * fitparam[NORM] - f;
}

std::vector<double> CompositeSymmetricNonresonantBackgroundFitModel::jacobian(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1), ret(nfit);
	double error, integral, intceps, intcgamma, intgminus;
	gsl_function func;
	func.params = &params;

	const double g = x[0];
	const double ceps = fitparam[CEPSILON];
	const double cgamma = fitparam[CGAMMA];
	const double gminus = fitparam[GMINUS];
	const double norm = fitparam[NORM];

	params[CEPSILON] = ceps;
	params[CGAMMA] = cgamma;
	params[GMINUS] = gminus;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	// evaluate the integrals
	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dcepsilon;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intceps, &error);

	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dcgamma;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intcgamma, &error);

	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the derivatives
	ret[CEPSILON] = -norm * intceps;

	ret[CGAMMA] = norm * intcgamma;

	ret[GMINUS] = -2.*norm / (k * gminus * gminus * M_SQRTPI) * intgminus;

	ret[NORM] = integral;

	return ret;
}

std::pair<double, std::vector<double>>
	CompositeSymmetricNonresonantBackgroundFitModel::resid_j(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	pair<double, vector<double>> ret;
	vector<double> params(nfit+1);
	double error, integral, intceps, intcgamma, intgminus;
	gsl_function func;
	func.params = &params;

	ret.second.resize(nfit);

	const double g = x[0];
	const double ceps = fitparam[CEPSILON];
	const double cgamma = fitparam[CGAMMA];
	const double gminus = fitparam[GMINUS];
	const double norm = fitparam[NORM];

	params[CEPSILON] = ceps;
	params[CGAMMA] = cgamma;
	params[GMINUS] = gminus;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	// evaluate the integrals
	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dcepsilon;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intceps, &error);

	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dcgamma;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intcgamma, &error);

	func.function = &CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the residual and derivatives
	ret.first = norm * integral - f;

	ret.second[CEPSILON] = -norm * intceps;

	ret.second[CGAMMA] = norm * intcgamma;

	ret.second[GMINUS] = -2.*norm / (k * gminus * gminus * M_SQRTPI) * intgminus;

	ret.second[NORM] = integral;

	return ret;
}

void CompositeSymmetricNonresonantBackgroundFitModel::append_default_guesses(
	std::list<std::vector<double>> &guess) const
{
	const list<double> list_ceps{40., 50., 75., 100., 150.},
		list_cgamma{4., 5., 6.},
		list_gminus{1.e-7, 1.e-6, 1.e-5};

	for(list<double>::const_iterator ceps = list_ceps.cbegin();
		ceps != list_ceps.cend(); ++ceps)
	{
	for(list<double>::const_iterator cgamma = list_cgamma.cbegin();
		cgamma != list_cgamma.cend(); ++cgamma)
	{
	for(list<double>::const_iterator gminus = list_gminus.cbegin();
		gminus != list_gminus.cend(); ++gminus)
	{

		vector<double> init(nfit);
		init[CEPSILON] = *ceps;
		init[CGAMMA] = *cgamma;
		init[GMINUS] = *gminus;
		init[NORM] = 1.;

		guess.emplace_back(init);
	}}}
}

void CompositeSymmetricNonresonantBackgroundFitModel::print_fit(std::ostream &out,
	const std::vector<double> &fitparam) const
{
	out << "cepsilon=" << scientific << setprecision(4) << fitparam[CEPSILON] <<
		", cgamma=" << scientific << setprecision(4) << fitparam[CGAMMA] <<
		", gminus=" << scientific << setprecision(4) << fitparam[GMINUS] <<
		", norm=" << scientific << setprecision(4) << fitparam[NORM];
}

void CompositeSymmetricNonresonantBackgroundFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const
{
	// sometimes both cepsilon and cgamma go negative
	if(fitparams[CEPSILON] < 0. && fitparams[CGAMMA] < 0.)
	{
		fitparams[CEPSILON] = -fitparams[CEPSILON];
		fitparams[CGAMMA] = -fitparams[CGAMMA];
	}
}

bool CompositeSymmetricNonresonantBackgroundFitModel::is_good_fit(
	const std::vector<double> &fitparams) const
{
	return (fitparams[CEPSILON] > 0. && fitparams[CGAMMA] > 0. &&
		fitparams[GMINUS] > 0.);
}

double CompositeSymmetricNonresonantBackgroundFitModel::int_p(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double ceps = fitparams[CEPSILON];
	const double cgamma = fitparams[CGAMMA];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = ceps*sqrt(temp1) - cgamma*sqrt(temp2);

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus)))
		* exp(-0.5 * temp3 * temp3 / temp2)
		/ (gp * sqrt(temp1 * temp2 * temp2 * temp2));
}

double CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dcepsilon(
	double gp, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double ceps = fitparams[CEPSILON];
	const double cgamma = fitparams[CGAMMA];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = ceps*sqrt(temp1) - cgamma*sqrt(temp2);

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus))) * temp3
		* exp(-0.5 * temp3 * temp3 / temp2)
		/ (gp * sqrt(temp2) * temp2 * temp2);
}

double CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dcgamma(
	double gp, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double ceps = fitparams[CEPSILON];
	const double cgamma = fitparams[CGAMMA];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = ceps*sqrt(temp1) - cgamma*sqrt(temp2);

	return (1. + gsl_sf_erf((gp - gminus) / (k * gminus))) * temp3
		* exp(-0.5 * temp3 * temp3 / temp2)
		/ (gp * sqrt(temp1) * temp2 * temp2);
}

double CompositeSymmetricNonresonantBackgroundFitModel::int_dp_dgminus(
	double gp, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double ceps = fitparams[CEPSILON];
	const double cgamma = fitparams[CGAMMA];
	const double gminus = fitparams[GMINUS];

	const double temp1 = g-gp;
	const double temp2 = 1.-g+gp;
	const double temp3 = ceps*sqrt(temp1) - cgamma*sqrt(temp2);
	const double temp4 = (gp - gminus) / (k * gminus);

	return exp(-temp4*temp4 -0.5 * temp3 * temp3 / temp2)
		/ sqrt(temp1 * temp2 * temp2 * temp2);
}

} // namespace molstat::transport
} // namespace molstat
