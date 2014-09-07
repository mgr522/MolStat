/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file uniform.cc
 * \brief Implementation of the uniform distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 * \endinternal
 */

#include <gsl/gsl_randist.h>
#include "uniform.h"

UniformDistribution::UniformDistribution(const double low, const double up)
	: RandomDistribution(), lower(low), upper(up) {}

double UniformDistribution::sample(shared_ptr<gsl_rng> r) const {
	return lower + (upper - lower) * gsl_rng_uniform(r.get());
}

std::string UniformDistribution::info() const {
	return "Uniform distribution between " + std::to_string(lower) + " and " +
		std::to_string(upper) + ".";
}
