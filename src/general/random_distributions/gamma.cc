/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file gamma.cc
 * \brief Implementation of the gamma distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 */

#include "gamma.h"

namespace molstat {

GammaDistribution::GammaDistribution(const double shape, const double scale)
	: RandomDistribution(), dist(shape, scale)
{
	if(shape <= 0. || scale <= 0.)
		throw std::invalid_argument("Gamma Distribution: The shape and scale " \
			"factors must be positive.");
}

double GammaDistribution::sample(Engine &engine) const
{
	return dist(engine);
}

std::string GammaDistribution::info() const
{
	return "Gamma: shape = " + std::to_string(dist.alpha()) + " and scale = " +
		std::to_string(dist.beta()) + ".";
}

} // namespace molstat
