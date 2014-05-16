/**
 * \file asymmetric_resonant.cc
 * \brief Implementation of the fitting model for resonant tunneling
 *        (asymmetric coupling).
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "asymmetric_resonant.h"

using namespace std;

AsymmetricResonantFitModel::AsymmetricResonantFitModel(
	const std::list<std::pair<std::array<double, 1>, double>> &data)
	: FitModel<1>(4, data), w(gsl_integration_workspace_alloc(nquad),
		&gsl_integration_workspace_free) {
}

double AsymmetricResonantFitModel::resid(const std::vector<double> &fitparam,
	const std::array<double, 1> &x, const double f) const {

	vector<double> params(nfit+1);
	double error, intmin, intmax, integral;
	gsl_function f;
	f.function = &AsymmetricResonantFitModel::int_p;
	f.params = &params;

	const double g = x[0];

	params[GAMMAL] = fitparam[GAMMAL];
	params[GAMMAR] = fitparam[GAMMAR];
	params[R] = fitparam[R];
	params[NORM] = fitparam[NORM];
	params[nfit] = g; // need to pass in the conductance value

	intmin = (2.0 - g - 2.0*sqrt(1.0-g)) / g;
	intmax = (2.0 - g + 2.0*sqrt(1.0-g)) / g;

	// calculate the integral
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, nquad, w, &integral,
		&error);
	integral *= norm / (g*sqrt(g));
	return integral - f;
}

std::vector<double> AsymmetricResonantFitModel::jacobian(
	const std::vector<double> &fitparam, const std::array<double, 1> &x,
	const double f) const {

	vector<double> params(nfit+1), ret(nfit);
	double error, intmin, intmax, integral, intgl, intgr, intr;
	gsl_function f;
	f.params = &params;

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
	f.function = &AsymmetricResonantFitModel::int_p;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
		&integral, &error);

	f.function = &AsymmetricResonantFitModel::int_dp_dgammal;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
		&intgl, &error);

	f.function = &AsymmetricResonantFitModel::int_dp_dgammar;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
		&intgr, &error);

	f.function = &AsymmetricResonantFitModel::int_dp_dr;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
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
	const double f) const {

	pair<double, vector<double>> ret;
	vector<double> params(nfit+1);
	double error, intmin, intmax, integral, intgl, intgr, intr;
	gsl_function f;
	f.params = &params;

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
	f.function = &AsymmetricResonantFitModel::int_p;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
		&integral, &error);

	f.function = &AsymmetricResonantFitModel::int_dp_dgammal;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
		&intgl, &error);

	f.function = &AsymmetricResonantFitModel::int_dp_dgammar;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
		&intgr, &error);

	f.function = &AsymmetricResonantFitModel::int_dp_dr;
	gsl_integration_qags(&f, intmin, intmax, 0.0, 1.0e-7, 2000, w,
		&intr, &error);

	// set the residual and derivatives
	ret.first = integral - f;

	ret.second[GAMMAL] = norm / (g*sqrt(g)) * intgl;

	ret.second[GAMMAR] = norm / (g*sqrt(g)) * intgr;

	ret.second[R] = -0.25 * norm * r * intr * (gammaL*gammaL + gammaR*gammaR)
		/ (g*g*sqrt(g));

	ret.second[NORM] = integral / (g*sqrt(g));

	return ret;
}

double AsymmetricResonantFitModel::int_p(double x, void *params) {
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];
	const double norm = fitparams[NORM];

	const double temp1 = 4.0*x - g*(1.0+x)*(1.0+x);
	const double temp2 = 1.0 + x*x;

	return x / (temp2 * sqrt(temp1 * temp2)) *
		(1.0 + (gammaL + x*gammaR)*(gammaL + x*gammaR) / temp2) *
		exp(-0.5*(x*gammaL - gammaR)*(x*gammaL - gammaR) / temp2
			-0.125*r*r*(gammaL*gammaL + gammaR*gammaR)*temp1 / (temp2 * g));
}

double AsymmetricResonantFitModel::int_dp_dr(double x, void *params) {
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];
	const double norm = fitparams[NORM];

	const double temp1 = 4.0*x - g*(1.0+x)*(1.0+x);
	const double temp2 = 1.0 + x*x;

	return x*sqrt(temp1 / temp2) / (temp2 * temp2) *
		(1.0 + (gammaL + x*gammaR)*(gammaL + x*gammaR) / temp2) *
		exp(-0.5*(x*gammaL - gammaR)*(x*gammaL - gammaR) / temp2
			-0.125*r*r*(gammaL*gammaL + gammaR*gammaR)*temp1 / (temp2 * g));
}

double AsymmetricResonantFitModel::int_dp_dgammaL(double x, void *params) {
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];
	const double norm = fitparams[NORM];

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

double AsymmetricResonantFitModel::int_dp_dgammaR(double x, void *params) {
	const vector<double> &fitparams = *(const vector<double>*)params;

	const double g = fitparams[4];
	const double gammaL = fitparams[GAMMAL];
	const double gammaR = fitparams[GAMMAR];
	const double r = fitparams[R];
	const double norm = fitparams[NORM];

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