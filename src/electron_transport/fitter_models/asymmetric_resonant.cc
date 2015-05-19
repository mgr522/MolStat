/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file asymmetric_resonant.cc
 * \brief Implementation of the fitting model for resonant tunneling
 *        (asymmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "asymmetric_resonant.h"
#include <iomanip>

using namespace std;

namespace molstat {
namespace transport {

std::vector<double> AsymmetricResonantFitModel::create_initial_guess(
	const std::map<std::string, double> &values) const
{
	vector<double> ret(4);

	try
	{
		ret[GAMMAL] = values.at("gammal");
		ret[GAMMAR] = values.at("gammar");
		ret[R] = values.at("r");
	}
	catch(const out_of_range &e)
	{
		throw invalid_argument("Initial guesses for the Asymmetric" \
			"ResonantFitModel must specify \"gammal\", \"gammar\", and \"r\"" \
			" parameters.");
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

AsymmetricResonantFitModel::AsymmetricResonantFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(4, data), w(gsl_integration_workspace_alloc(nquad),
		&gsl_integration_workspace_free) 
{}

double AsymmetricResonantFitModel::resid(const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const
{
	vector<double> params(nfit+1);
	double error, intmin, intmax, integral;
	gsl_function func;
	func.function = &AsymmetricResonantFitModel::int_p;
	func.params = &params;

	const double g = x[0];

	params[GAMMAL] = fitparam[GAMMAL];
	params[GAMMAR] = fitparam[GAMMAR];
	params[R] = fitparam[R];
	params[NORM] = fitparam[NORM];
	params[nfit] = g; // need to pass in the conductance value

	intmin = (2.0 - g - 2.0*sqrt(1.0-g)) / g;
	intmax = (2.0 - g + 2.0*sqrt(1.0-g)) / g;

	// calculate the integral
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);
	integral *= fitparam[NORM] / (g*sqrt(g));
	return integral - f;
}

std::vector<double> AsymmetricResonantFitModel::jacobian(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	vector<double> params(nfit+1), ret(nfit);
	double error, intmin, intmax, integral, intgl, intgr, intr;
	gsl_function func;
	func.params = &params;

	const double g = x[0];
	const double gammaL = fitparam[GAMMAL];
	const double gammaR = fitparam[GAMMAR];
	const double r = fitparam[R];
	const double norm = fitparam[NORM];

	params[GAMMAL] = gammaL;
	params[GAMMAR] = gammaR;
	params[R] = r;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	intmin = (2.0 - g - 2.0*sqrt(1.0-g)) / g;
	intmax = (2.0 - g + 2.0*sqrt(1.0-g)) / g;

	// evaluate the four integrals
	func.function = &AsymmetricResonantFitModel::int_p;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &AsymmetricResonantFitModel::int_dp_dgammaL;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&intgl, &error);

	func.function = &AsymmetricResonantFitModel::int_dp_dgammaR;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&intgr, &error);

	func.function = &AsymmetricResonantFitModel::int_dp_dr;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&intr, &error);

	// set the derivatives
	ret[GAMMAL] = norm / (g*sqrt(g)) * intgl;

	ret[GAMMAR] = norm / (g*sqrt(g)) * intgr;

	ret[R] = -0.25 * norm * r * intr * (gammaL*gammaL + gammaR*gammaR)
		/ (g*g*sqrt(g));

	ret[NORM] = integral / (g*sqrt(g));

	return ret;
}

std::pair<double, std::vector<double>> AsymmetricResonantFitModel::resid_j(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const
{
	pair<double, vector<double>> ret;
	vector<double> params(nfit+1);
	double error, intmin, intmax, integral, intgl, intgr, intr;
	gsl_function func;
	func.params = &params;

	ret.second.resize(nfit);

	const double g = x[0];
	const double gammaL = fitparam[GAMMAL];
	const double gammaR = fitparam[GAMMAR];
	const double r = fitparam[R];
	const double norm = fitparam[NORM];

	params[GAMMAL] = gammaL;
	params[GAMMAR] = gammaR;
	params[R] = r;
	params[NORM] = norm;
	params[nfit] = g; // need to pass in the conductance value

	intmin = (2.0 - g - 2.0*sqrt(1.0-g)) / g;
	intmax = (2.0 - g + 2.0*sqrt(1.0-g)) / g;

	// evaluate the four integrals
	func.function = &AsymmetricResonantFitModel::int_p;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&integral, &error);

	func.function = &AsymmetricResonantFitModel::int_dp_dgammaL;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&intgl, &error);

	func.function = &AsymmetricResonantFitModel::int_dp_dgammaR;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&intgr, &error);

	func.function = &AsymmetricResonantFitModel::int_dp_dr;
	gsl_integration_qags(&func, intmin, intmax, 0.0, 1.0e-7, nquad, w.get(),
		&intr, &error);

	// set the residual and derivatives
	ret.first = norm * integral / (g*sqrt(g)) - f;

	ret.second[GAMMAL] = norm / (g*sqrt(g)) * intgl;

	ret.second[GAMMAR] = norm / (g*sqrt(g)) * intgr;

	ret.second[R] = -0.25 * norm * r * intr * (gammaL*gammaL + gammaR*gammaR)
		/ (g*g*sqrt(g));

	ret.second[NORM] = integral / (g*sqrt(g));

	return ret;
}

void AsymmetricResonantFitModel::append_default_guesses(
	std::list<std::vector<double>> &guess) const
{
	const list<double> list_gamma{5., 10., 20., 30., 40.},
		list_r{0.1, 0.5, 1., 2., 10.};

	for(list<double>::const_iterator gammaL = list_gamma.cbegin();
		gammaL != list_gamma.cend(); ++gammaL)
	{
	for(list<double>::const_iterator gammaR = list_gamma.cbegin();
		gammaR != list_gamma.cend(); ++gammaR)
	{
	for(list<double>::const_iterator r = list_r.cbegin();
		r != list_r.cend(); ++r)
	{

		vector<double> init(nfit);
		init[GAMMAL] = *gammaL;
		init[GAMMAR] = *gammaR;
		init[R] = *r;
		init[NORM] = 1.;

		guess.emplace_back(init);
	}}}
}

void AsymmetricResonantFitModel::print_fit(std::ostream &out,
	const std::vector<double> &fitparam) const
{
	out << "gammaL=" << scientific << setprecision(4) << fitparam[GAMMAL] <<
		", gammaR=" << scientific << setprecision(4) << fitparam[GAMMAR] <<
		", r=" << scientific << setprecision(4) << fitparam[R] << ", norm=" <<
		scientific << setprecision(4) << fitparam[NORM];
}

void AsymmetricResonantFitModel::process_fit_parameters(
	std::vector<double> &fitparams) const
{
	// sometimes both gammas go negative...
	if(fitparams[GAMMAL] < 0. && fitparams[GAMMAR] < 0.)
	{
		fitparams[GAMMAL] = -fitparams[GAMMAL];
		fitparams[GAMMAR] = -fitparams[GAMMAR];
	}

	// for convenience, make sure gammaL is smaller than gammaR
	if(fitparams[GAMMAL] > fitparams[GAMMAR])
		swap(fitparams[GAMMAL], fitparams[GAMMAR]);

	// sometimes r is negative
	if(fitparams[R] < 0.)
		fitparams[R] = -fitparams[R];
}

double AsymmetricResonantFitModel::int_p(double x, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];

	const double temp1 = 4.0*x - g*(1.0+x)*(1.0+x);
	const double temp2 = 1.0 + x*x;

	return x / (temp2 * sqrt(temp1 * temp2)) *
		(1.0 + (gammaL + x*gammaR)*(gammaL + x*gammaR) / temp2) *
		exp(-0.5*(x*gammaL - gammaR)*(x*gammaL - gammaR) / temp2
			-0.125*r*r*(gammaL*gammaL + gammaR*gammaR)*temp1 / (temp2 * g));
}

double AsymmetricResonantFitModel::int_dp_dr(double x, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];

	const double temp1 = 4.0*x - g*(1.0+x)*(1.0+x);
	const double temp2 = 1.0 + x*x;

	return x*sqrt(temp1 / temp2) / (temp2 * temp2) *
		(1.0 + (gammaL + x*gammaR)*(gammaL + x*gammaR) / temp2) *
		exp(-0.5*(x*gammaL - gammaR)*(x*gammaL - gammaR) / temp2
			-0.125*r*r*(gammaL*gammaL + gammaR*gammaR)*temp1 / (temp2 * g));
}

double AsymmetricResonantFitModel::int_dp_dgammaL(double x, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];

	const double temp1 = 4.0*x - g*(1.0+x)*(1.0+x);
	const double temp2 = 1.0 + x*x;
	double temp3 = gammaL + x*gammaR;
	temp3 *= temp3 / temp2;

	return x / (temp2 * temp2 * sqrt(temp2)) *
		(((2.0-x*x)*gammaL + 3.0*x*gammaR - temp3*x*(x*gammaL-gammaR))
			/ sqrt(temp1)
			-0.25*(1.0 + temp3)*r*r*gammaL*sqrt(temp1)/g) *
		exp(-0.5*(x*gammaL - gammaR)*(x*gammaL - gammaR) / temp2
			-0.125*r*r*(gammaL*gammaL + gammaR*gammaR)*temp1 / (temp2 * g));
}

double AsymmetricResonantFitModel::int_dp_dgammaR(double x, void *params)
{
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];

	const double temp1 = 4.0*x - g*(1.0+x)*(1.0+x);
	const double temp2 = 1.0 + x*x;
	double temp3 = gammaL + x*gammaR;
	temp3 *= temp3 / temp2;

	return x / (temp2 * temp2 * sqrt(temp2)) *
		((3.0*x*gammaL + (2.0*x*x-1.0)*gammaR + temp3*(x*gammaL-gammaR))
			/ sqrt(temp1)
			-0.25*(1.0 + temp3)*r*r*gammaR*sqrt(temp1)/g) *
		exp(-0.5*(x*gammaL - gammaR)*(x*gammaL - gammaR) / temp2
			-0.125*r*r*(gammaL*gammaL + gammaR*gammaR)*temp1 / (temp2 * g));
}

} // namespace molstat::transport
} // namespace molstat
