/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file simulator_models/asymmetric_voltage_one_site.cc
 * \brief Implementation for the voltage-dependent, asymmetric-coupling,
 * tight-binding model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 * \endinternal
 */

#include "asymmetric_voltage_one_site.h"
#include "asymmetric_one_site.h"
#include <cmath>

using namespace std;

AsymmetricVoltageOneSiteModel::AsymmetricVoltageOneSiteModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gammaL,
	const shared_ptr<const RandomDistribution> &gammaR,
	const shared_ptr<const RandomDistribution> &a)
	: dist_eps(eps), dist_gammaL(gammaL), dist_gammaR(gammaR), dist_a(a) {}

double AsymmetricVoltageOneSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta, const double V)
	const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);
	double a = dist_a->sample(r);

	return AsymmetricVoltageOneSiteModel::static_conductance(EF, V, eta,
		eps, gammaL, gammaR, a);
}

double AsymmetricVoltageOneSiteModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta, const double V)
	const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);
	double a = dist_a->sample(r);

	return AsymmetricVoltageOneSiteModel::diff_conductance(EF, V, eta, eps,
		gammaL, gammaR, a);
}

double AsymmetricVoltageOneSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);
	double a = dist_a->sample(r);

	return AsymmetricVoltageOneSiteModel::transmission(EF, 0., eps, gammaL,
		gammaR, a);
}

double AsymmetricVoltageOneSiteModel::transmission(const double E,
	const double V, const double eps, const double gammaL, const double gammaR,
	const double a) {

	return AsymmetricOneSiteModel::transmission(E, eps+a*V, gammaL, gammaR);
}

double AsymmetricVoltageOneSiteModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR, const double a) {

	return AsymmetricOneSiteModel::static_conductance(EF, V, eta, eps+a*V,
		gammaL, gammaR);
}

double AsymmetricVoltageOneSiteModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR, const double a) {

	return (eta-a) * AsymmetricVoltageOneSiteModel::transmission(EF + eta*V, V,
		eps, gammaL, gammaR, a) + (a+1.-eta) *
		AsymmetricVoltageOneSiteModel::transmission(EF + (eta-1.)*V, V, eps,
		gammaL, gammaR, a);
}
