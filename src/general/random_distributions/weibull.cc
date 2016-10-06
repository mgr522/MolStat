/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2016 Stony Brook University. */

/**
 * \file weibull.cc
 * \brief Implementation of the Weibull distribution.
 *
 * \author Matthew G.\ Reuter
 * \date October 2016
 */

#include "weibull.h"

namespace molstat {

WeibullDistribution::WeibullDistribution(const double shape, const double scale)
	: RandomDistribution(), dist(shape, scale)
{
	if(scale <= 0.)
		throw std::invalid_argument("Weibull Distribution: The scale factor" \
			" must be positive.");
	if(shape <= 0.)
		throw std::invalid_argument("Weibull Distribution: The shape factor" \
			" must be positive.");
}

double WeibullDistribution::sample(Engine &engine) const
{
	return dist(engine);
}

std::string WeibullDistribution::info() const
{
	return "Weibull: shape = " + std::to_string(dist.a()) + " and scale = " +
		std::to_string(dist.b()) + ".";
}

} // namespace molstat
