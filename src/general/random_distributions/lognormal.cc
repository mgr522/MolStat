/**
 * \file lognormal.cc
 * \brief Implementation of the lognormal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <gsl/gsl_randist.h>
#include "lognormal.h"

LognormalDistribution::LognormalDistribution(const double zeta_,
	const double sigma_)
	: RandomDistribution(), zeta(zeta_), sigma(sigma_) {}

double LognormalDistribution::sample(shared_ptr<gsl_rng> r) const {
	return gsl_ran_lognormal(r.get(), zeta, sigma);
}
