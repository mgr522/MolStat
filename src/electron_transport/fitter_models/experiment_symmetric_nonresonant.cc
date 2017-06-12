/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file experiment_symmetric_nonresonant.cc
 * \brief Implementation of the fitting model for nonresonant tunneling
 *    (symmetric coupling) combined with background tunneling (for describing
 *    experiment).
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#include "experiment_symmetric_nonresonant.h"
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

using namespace std;

namespace molstat {
namespace transport {

std::vector<double> ExperimentSymmetricNonresonantFitModel
	::create_initial_guess(const std::map<std::string, double> &values) const
{
	vector<double> ret(6);

	try
	{
		ret[CEPSILON] = values.at("cepsilon");
		ret[CGAMMA] = values.at("cgamma");
		ret[GMINUS] = values.at("gminus");
		ret[NSIGNAL] = values.at("nsignal");
		ret[NBACKGROUND] = values.at("nbackground");
		ret[NBASELINE] = values.at("nbaseline");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the ExperimentSymmetric" \
			"NonresonantFitModel must specify \"cepsilon\", \"cgamma\", \"gminus" \
			"\", \"nsignal\", \"nbackground\", and \"nbaseline\" parameters.");
	}

	return ret;
}

ExperimentSymmetricNonresonantFitModel
	::ExperimentSymmetricNonresonantFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(6, data), w(gsl_integration_workspace_alloc(nquad),
		&gsl_integration_workspace_free) 
{}

double ExperimentSymmetricNonresonantFitModel::resid(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1);
	double error, integral;
	gsl_function func;
	func.function = &ExperimentSymmetricNonresonantFitModel::int_p;
	func.params = &params;

	const double g = x[0];

	params[CEPSILON] = fitparam[CEPSILON];
	params[CGAMMA] = fitparam[CGAMMA];
	params[GMINUS] = fitparam[GMINUS];
	params[NSIGNAL] = fitparam[NSIGNAL];
	params[NBACKGROUND] = fitparam[NBACKGROUND];
	params[NBASELINE] = fitparam[NBASELINE];
	params[nfit] = g; // need to pass in the conductance value

	// calculate the integral
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad,
		w.get(), &integral, &error);

	return integral * fitparam[NSIGNAL] + fitparam[NBACKGROUND] / g
		+ fitparam[NBASELINE] - f;
}

std::vector<double> ExperimentSymmetricNonresonantFitModel
	::jacobian(const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const
{
	vector<double> params(nfit+1), ret(nfit);
	double error, integral, intceps, intcgamma, intgminus;
	gsl_function func;
	func.params = &params;

	const double g = x[0];
	const double ceps = fitparam[CEPSILON];
	const double cgamma = fitparam[CGAMMA];
	const double gminus = fitparam[GMINUS];
	const double nsignal = fitparam[NSIGNAL];
	const double nbackground = fitparam[NBACKGROUND];
	const double nbaseline = fitparam[NBASELINE];

	params[CEPSILON] = ceps;
	params[CGAMMA] = cgamma;
	params[GMINUS] = gminus;
	params[NSIGNAL] = nsignal;
	params[NBACKGROUND] = nbackground;
	params[NBASELINE] = nbaseline;
	params[nfit] = g; // need to pass in the conductance value

	func.function = &ExperimentSymmetricNonresonantFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &ExperimentSymmetricNonresonantFitModel::int_dp_dcepsilon;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intceps, &error);

	func.function = &ExperimentSymmetricNonresonantFitModel::int_dp_dcgamma;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intcgamma, &error);

	func.function = &ExperimentSymmetricNonresonantFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the derivatives
	ret[CEPSILON] = -nsignal * intceps;

	ret[CGAMMA] = nsignal * intcgamma;

	ret[GMINUS] = -2.*nsignal / (k * gminus * gminus * M_SQRTPI) * intgminus;

	ret[NSIGNAL] = integral;

	ret[NBACKGROUND] = 1. / g;

	ret[NBASELINE] = 1.;

	return ret;
}

std::pair<double, std::vector<double>>
	ExperimentSymmetricNonresonantFitModel::resid_j(
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
	const double nsignal = fitparam[NSIGNAL];
	const double nbackground = fitparam[NBACKGROUND];
	const double nbaseline = fitparam[NBASELINE];

	params[CEPSILON] = ceps;
	params[CGAMMA] = cgamma;
	params[GMINUS] = gminus;
	params[NSIGNAL] = nsignal;
	params[NBACKGROUND] = nbackground;
	params[NBASELINE] = nbaseline;
	params[nfit] = g; // need to pass in the conductance value

	// evaluate the integrals
	func.function = &ExperimentSymmetricNonresonantFitModel::int_p;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &ExperimentSymmetricNonresonantFitModel::int_dp_dcepsilon;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intceps, &error);

	func.function = &ExperimentSymmetricNonresonantFitModel::int_dp_dcgamma;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intcgamma, &error);

	func.function = &ExperimentSymmetricNonresonantFitModel::int_dp_dgminus;
	gsl_integration_qags(&func, 0., g, 0.0, 1.0e-7, nquad, w.get(),
		&intgminus, &error);

	// set the residual and derivatives
	ret.first = nsignal * integral + nbackground / g + nbaseline - f;

	ret.second[CEPSILON] = -nsignal * intceps;

	ret.second[CGAMMA] = nsignal * intcgamma;

	ret.second[GMINUS] = -2.*nsignal / (k * gminus * gminus * M_SQRTPI) * intgminus;

	ret.second[NSIGNAL] = integral;

	ret.second[NBACKGROUND] = 1. / g;

	ret.second[NBASELINE] = 1.;

	return ret;
}

void ExperimentSymmetricNonresonantFitModel
	::append_default_guesses(std::list<std::vector<double>> &guess) const
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
		init[NSIGNAL] = init[NBACKGROUND] = 1.;
		init[NBASELINE] = 0.;

		guess.emplace_back(init);
	}}}
}

void ExperimentSymmetricNonresonantFitModel
	::print_fit(std::ostream &out, const std::vector<double> &fitparam) const
{
	out << "cepsilon=" << scientific << setprecision(4) << fitparam[CEPSILON] <<
		", cgamma=" << scientific << setprecision(4) << fitparam[CGAMMA] <<
		", gminus=" << scientific << setprecision(4) << fitparam[GMINUS] <<
		", nsignal=" << scientific << setprecision(4) << fitparam[NSIGNAL] <<
		", nbackground=" << scientific << setprecision(4) << fitparam[NBACKGROUND] <<
		", nbaseline=" << scientific << setprecision(4) << fitparam[NBASELINE];
}

void ExperimentSymmetricNonresonantFitModel
	::process_fit_parameters(std::vector<double> &fitparams) const
{
	// sometimes both cepsilon and cgamma go negative
	if(fitparams[CEPSILON] < 0. && fitparams[CGAMMA] < 0.)
	{
		fitparams[CEPSILON] = -fitparams[CEPSILON];
		fitparams[CGAMMA] = -fitparams[CGAMMA];
	}
}

bool ExperimentSymmetricNonresonantFitModel::is_good_fit(
	const std::vector<double> &fitparams) const
{
	return (fitparams[CEPSILON] > 0. && fitparams[CGAMMA] > 0.
		&& fitparams[GMINUS] > 0.);
}

double ExperimentSymmetricNonresonantFitModel::int_p(double gp,
	void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[6];
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

double ExperimentSymmetricNonresonantFitModel::int_dp_dcepsilon(
	double gp, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[6];
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

double ExperimentSymmetricNonresonantFitModel::int_dp_dcgamma(
	double gp, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[6];
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

double ExperimentSymmetricNonresonantFitModel::int_dp_dgminus(
	double gp, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[6];
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
