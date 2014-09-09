/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file simulator_models/symmetric_voltage_one_site.cc
 * \brief Implementation of the symmetric-coupling, voltage-dependent one-site
 *    tight-binding model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date July 2014
 * \endinternal
 */

#include "symmetric_voltage_one_site.h"
#include "symmetric_one_site.h"
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

SymmetricVoltageOneSiteModel::SymmetricVoltageOneSiteModel(
	const shared_ptr<const RandomDistribution> &eps,
	const shared_ptr<const RandomDistribution> &gamma,
	const shared_ptr<const RandomDistribution> &a)
	: dist_eps(eps), dist_gamma(gamma), dist_a(a) {}

double SymmetricVoltageOneSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double a = dist_a->sample(r);

	return SymmetricVoltageOneSiteModel::static_conductance(EF, V, eta, eps,
		gamma, a);
}

double SymmetricVoltageOneSiteModel::diff_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double a = dist_a->sample(r);

	return SymmetricVoltageOneSiteModel::diff_conductance(EF, V, eta, eps,
		gamma, a);
}

double SymmetricVoltageOneSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	// get model parameters from the random distributions
	double eps = dist_eps->sample(r);
	double gamma = dist_gamma->sample(r);
	double a = dist_a->sample(r);

	return SymmetricVoltageOneSiteModel::transmission(EF, 0., eps, gamma, a);
}

double SymmetricVoltageOneSiteModel::transmission(const double E,
	const double V, const double eps, const double gamma, const double a) {

	return SymmetricOneSiteModel::transmission(E, eps+a*V, gamma);
}

double SymmetricVoltageOneSiteModel::static_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma,
	const double a) {

	return SymmetricOneSiteModel::static_conductance(EF, V, eta,
		eps+a*V, gamma);
}

double SymmetricVoltageOneSiteModel::diff_conductance(const double EF,
	const double V, const double eta, const double eps, const double gamma,
	const double a) {

	return (eta-a) * SymmetricVoltageOneSiteModel::transmission(EF + eta*V, V,
		eps, gamma, a) + (a+1.-eta) * SymmetricVoltageOneSiteModel::transmission(
		EF + (eta-1.)*V, V, eps, gamma, a);
}