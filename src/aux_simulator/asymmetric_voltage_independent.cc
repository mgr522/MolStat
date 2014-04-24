/**
 * \file asymmetric_voltage_independent.cc
 * \brief Implementation for the asymmetric-coupling, voltage-independent
 * tight-binding model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include "asymmetric_voltage_independent.h"
#include <cmath>

using namespace std;

AsymmetricVoltageIndependentModel::AsymmetricVoltageIndependentModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gammaL,
	const shared_ptr<const RandomDistribution> &gammaR)
	: dist_eps(eps), dist_gammaL(gammaL), dist_gammaR(gammaR) {}

double AsymmetricVoltageIndependentModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta, const double V)
	const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);

	return AsymmetricVoltageIndependentModel::static_conductance(EF, V, eta,
		eps, gammaL, gammaR);
}

double AsymmetricVoltageIndependentModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta, const double V)
	const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);

	return AsymmetricVoltageIndependentModel::diff_conductance(EF, V, eta, eps,
		gammaL, gammaR);
}

double AsymmetricVoltageIndependentModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);

	return AsymmetricVoltageIndependentModel::transmission(EF, eps, gammaL,
		gammaR);
}

double AsymmetricVoltageIndependentModel::transmission(const double E,
	const double eps, const double gammaL, const double gammaR) {

	return 4.*gammaL*gammaR /
		(4.*(E - eps)*(E - eps) + (gammaL + gammaR)*(gammaL + gammaR));
}

double AsymmetricVoltageIndependentModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR) {

	double gsum = gammaL + gammaR;

	return 2.*gammaL*gammaR / (V * gsum) *
		(atan(2.*(EF-eps+eta*V) / gsum) - atan(2.*(EF-eps+(eta-1.)*V) / gsum));
}

double AsymmetricVoltageIndependentModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR) {

	return eta * AsymmetricVoltageIndependentModel::transmission(EF + eta*V,
		eps, gammaL, gammaR) + (1.-eta) * AsymmetricVoltageIndependentModel::
		transmission(EF + (eta-1.)*V, eps, gammaL, gammaR);
}
