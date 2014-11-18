/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file normal.cc
 * \brief Implementation of the normal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 * \endinternal
 */

#include "normal.h"

namespace molstat {

NormalDistribution::NormalDistribution(const double mean, const double stdev)
	: RandomDistribution(), dist(mean, stdev)
{
	if(stdev <= 0.)
		throw std::invalid_argument("Normal Distribution: The standard " \
			"deviation must be positive.");
}

double NormalDistribution::sample(Engine &engine) const
{
	return dist(engine);
}

std::string NormalDistribution::info() const
{
	return "Normal: mean = " + std::to_string(dist.mean()) + " and stdev = " +
		std::to_string(dist.stddev()) + ".";
}

} // namespace molstat
