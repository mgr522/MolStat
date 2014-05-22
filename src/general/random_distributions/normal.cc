/**
 * \file normal.cc
 * \brief Implementation of the normal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <gsl/gsl_randist.h>
#include "normal.h"

NormalDistribution::NormalDistribution(const double mean_, const double stdev_)
	: RandomDistribution(), mean(mean_), stdev(stdev_) {}

double NormalDistribution::sample(shared_ptr<gsl_rng> r) const {
	return gsl_ran_gaussian(r.get(), stdev) + mean;
}
