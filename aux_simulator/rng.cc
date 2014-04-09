/**
 * \file rng.cc
 * \brief Implementation of the interface for random number generation.
 *
 * Provides a common interface for random number generation with GSL.
 * The following distributions are implemented:
 *    - Constant distribution (the value is fixed).
 *    - Uniform distribution.
 *    - Normal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include "rng.h"
#include <gsl/gsl_randist.h>

ConstantDistribution::ConstantDistribution(const double val)
	: RandomDistribution(), value(val) {}

double ConstantDistribution::sample(shared_ptr<gsl_rng> r) const {
	return value;
}

UniformDistribution::UniformDistribution(const double low, const double up)
	: RandomDistribution(), lower(low), upper(up) {}

double UniformDistribution::sample(shared_ptr<gsl_rng> r) const {
	return lower + (upper - lower) * gsl_rng_uniform(r.get());
}

NormalDistribution::NormalDistribution(const double mean_, const double stdev_)
	: RandomDistribution(), mean(mean_), stdev(stdev_) {}

double NormalDistribution::sample(shared_ptr<gsl_rng> r) const {
	return gsl_ran_gaussian(r.get(), stdev) + mean;
}
