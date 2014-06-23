/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file simulator_models/asymmetric_one_site.cc
 * \brief Implementation for the asymmetric-coupling, tight-binding model for
 * calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 * \endinternal
 */

#include "asymmetric_one_site.h"
#include <cmath>

using namespace std;

AsymmetricOneSiteModel::AsymmetricOneSiteModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gammaL,
	const shared_ptr<const RandomDistribution> &gammaR)
	: dist_eps(eps), dist_gammaL(gammaL), dist_gammaR(gammaR) {}

double AsymmetricOneSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta, const double V)
	const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);

	return AsymmetricOneSiteModel::static_conductance(EF, V, eta,
		eps, gammaL, gammaR);
}

double AsymmetricOneSiteModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta, const double V)
	const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);

	return AsymmetricOneSiteModel::diff_conductance(EF, V, eta, eps, gammaL,
		gammaR);
}

double AsymmetricOneSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);

	return AsymmetricOneSiteModel::transmission(EF, eps, gammaL, gammaR);
}

double AsymmetricOneSiteModel::transmission(const double E,
	const double eps, const double gammaL, const double gammaR) {

	return 4.*gammaL*gammaR /
		(4.*(E - eps)*(E - eps) + (gammaL + gammaR)*(gammaL + gammaR));
}

double AsymmetricOneSiteModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR) {

	double gsum = gammaL + gammaR;

	return 2.*gammaL*gammaR / (V * gsum) *
		(atan(2.*(EF-eps+eta*V) / gsum) - atan(2.*(EF-eps+(eta-1.)*V) / gsum));
}

double AsymmetricOneSiteModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR) {

	return eta * AsymmetricOneSiteModel::transmission(EF + eta*V, eps, gammaL,
		gammaR) + (1.-eta) * AsymmetricOneSiteModel::transmission(
		EF + (eta-1.)*V, eps, gammaL, gammaR);
}
