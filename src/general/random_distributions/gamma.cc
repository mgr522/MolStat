/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file gamma.cc
 * \brief Implementation of the gamma distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 * \endinternal
 */

#include <gsl/gsl_randist.h>
#include "gamma.h"

GammaDistribution::GammaDistribution(const double shape_, const double scale_)
	: RandomDistribution(), shape(shape_), scale(scale_) {}

double GammaDistribution::sample(gsl_rng_ptr &r) const {
	return gsl_ran_gamma(r.get(), shape, scale);
}

std::string GammaDistribution::info() const {
	return "Gamma distribution with shape and scale factors of " +
		std::to_string(shape) + " and " + std::to_string(scale) + ".";
}
