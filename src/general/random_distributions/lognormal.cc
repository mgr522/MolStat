/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file lognormal.cc
 * \brief Implementation of the lognormal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 */

#include "lognormal.h"

namespace molstat {

LognormalDistribution::LognormalDistribution(const double zeta,
	const double sigma)
	: RandomDistribution(), dist(zeta, sigma)
{
	if(sigma <= 0.)
		throw std::invalid_argument("Lognormal Distribution: The standard " \
				"deviation (sigma) must be positive.");
}

double LognormalDistribution::sample(Engine &engine) const
{
	return dist(engine);
}

std::string LognormalDistribution::info() const
{
	return "Lognormal: mean = " + std::to_string(dist.m()) +
		" and stdev = " + std::to_string(dist.s()) + " (log space).";
}

} // namespace molstat
