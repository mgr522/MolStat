/**
 * \file gamma.cc
 * \brief Implementation of the gamma distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include <gsl/gsl_randist.h>
#include "gamma.h"

GammaDistribution::GammaDistribution(const double shape_, const double scale_)
	: RandomDistribution(), shape(shape_), scale(scale_) {}

double GammaDistribution::sample(shared_ptr<gsl_rng> r) const {
	return gsl_ran_gamma(r.get(), shape, scale);
}
