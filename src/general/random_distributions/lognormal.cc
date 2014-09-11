/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file lognormal.cc
 * \brief Implementation of the lognormal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 * \endinternal
 */

#include <gsl/gsl_randist.h>
#include "lognormal.h"

LognormalDistribution::LognormalDistribution(const double zeta_,
	const double sigma_)
	: RandomDistribution(), zeta(zeta_), sigma(sigma_) {}

double LognormalDistribution::sample(shared_ptr<gsl_rng> r) const {
	return gsl_ran_lognormal(r.get(), zeta, sigma);
}

std::string LognormalDistribution::info() const {
	return "Lognormal distribution with mean " + std::to_string(zeta) +
		" and standard deviation " + std::to_string(sigma) + " in log space.";
}
