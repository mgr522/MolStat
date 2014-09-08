/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file independent_voltage_one_site.cc
 * \brief A sum of two independent, symmetric-coupling, voltage-dependent
 *    one-site tight-binding models for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date July 2014
 * \endinternal
 */

#include "independent_voltage_two_site.h"

using namespace std;

IndependentVoltageTwoSiteModel::IndependentVoltageTwoSiteModel(
	const shared_ptr<const RandomDistribution> &eps1,
	const shared_ptr<const RandomDistribution> &gamma1,
	const shared_ptr<const RandomDistribution> &a1,
	const shared_ptr<const RandomDistribution> &eps2,
	const shared_ptr<const RandomDistribution> &gamma2,
	const shared_ptr<const RandomDistribution> &a2)
	: channel1(eps1, gamma1, a1), channel2(eps2, gamma2, a2) {}

double IndependentVoltageTwoSiteModel::static_conductance(
	shared_ptr<gsl_rng> r, const double EF, const double eta,
	const double V) const {

	return channel1.static_conductance(r, EF, eta, V) +
		channel2.static_conductance(r, EF, eta, V);
}

double IndependentVoltageTwoSiteModel::diff_conductance(shared_ptr<gsl_rng> r,
	const double EF, const double eta, const double V) const {

	return channel1.diff_conductance(r, EF, eta, V) +
		channel2.diff_conductance(r, EF, eta, V);
}

double IndependentVoltageTwoSiteModel::zero_bias_conductance(
	shared_ptr<gsl_rng> r, const double EF) const {

	return channel1.zero_bias_conductance(r, EF) +
		channel2.zero_bias_conductance(r, EF);
}
