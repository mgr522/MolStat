/**
 * \file simulator_models/asymmetric_two_site.cc
 * \brief Implementation of the asymmetric-coupling, two-site tight-binding
 *    model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "asymmetric_two_site.h"
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

AsymmetricTwoSiteModel::AsymmetricTwoSiteModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gammaL,
	const shared_ptr<const RandomDistribution> &gammaR,
	const shared_ptr<const RandomDistribution> &beta)
	: dist_eps(eps), dist_gammaL(gammaL), dist_gammaR(gammaR), dist_beta(beta) {
}

double AsymmetricTwoSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);
	double beta = dist_beta->sample(r);

	return AsymmetricTwoSiteModel::static_conductance(EF, V, eta, eps,
		gammaL, gammaR, beta);
}

double AsymmetricTwoSiteModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);
	double beta = dist_beta->sample(r);

	return AsymmetricTwoSiteModel::diff_conductance(EF, V, eta, eps,
		gammaL, gammaR, beta);
}

double AsymmetricTwoSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gammaL = dist_gammaL->sample(r);
	double gammaR = dist_gammaR->sample(r);
	double beta = dist_beta->sample(r);

	return AsymmetricTwoSiteModel::transmission(EF, eps, gammaL, gammaR,
		beta);
}

double AsymmetricTwoSiteModel::transmission(const double E, const double eps,
	const double gammaL, const double gammaR, const double beta) {

	double temp = 4.*(E-eps)*(E-eps) - 4.*beta*beta - gammaL*gammaR;

	return 16.*gammaL*gammaR*beta*beta /
		(temp*temp + 4.*(gammaL+gammaR)*(gammaL+gammaR)*(E-eps)*(E-eps));
}

double AsymmetricTwoSiteModel::static_c_integral(const double z,
	const double eps, const double gammaL, const double gammaR,
	const double beta) {

	complex<double> bgg = sqrt(complex<double>((gammaL-gammaR)*(gammaL-gammaR)
		- 16.*beta*beta, 0.));

	complex<double> denom1 = sqrt(-8.*beta*beta + gammaL*gammaL + gammaR*gammaR
		- (gammaL + gammaR) * bgg);

	complex<double> denom2 = sqrt(-8.*beta*beta + gammaL*gammaL + gammaR*gammaR
		+ (gammaL + gammaR) * bgg);

	return sqrt(128.)*gammaL*gammaR*beta*beta / (gammaL + gammaR) * real(
		(atan(sqrt(8.)*(z-eps) / denom1) / denom1 -
		 atan(sqrt(8.)*(z-eps) / denom2) / denom2) / bgg
		);
}

double AsymmetricTwoSiteModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR, const double beta) {

	return (AsymmetricTwoSiteModel::static_c_integral(EF+eta*V, eps, gammaL,
		gammaR, beta) - AsymmetricTwoSiteModel::static_c_integral(
		EF+(eta-1.)*V, eps, gammaL, gammaR, beta)) / V;
}

double AsymmetricTwoSiteModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gammaL,
	const double gammaR, const double beta) {

	return eta * AsymmetricTwoSiteModel::transmission(EF + eta*V, eps, gammaL,
		gammaR, beta)

		+ (1.-eta) * AsymmetricTwoSiteModel::transmission(EF + (eta-1.)*V,
		eps, gammaL, gammaR, beta);
}
