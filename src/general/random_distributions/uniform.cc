/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file uniform.cc
 * \brief Implementation of the uniform distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 * \endinternal
 */

#include "uniform.h"

namespace molstat {

UniformDistribution::UniformDistribution(const double low, const double up)
	: RandomDistribution(), dist(low, up)
{
	if(low >= up)
		throw std::invalid_argument("Uniform distribution: The lower bound " \
			"must be lower than the upper bound.");
}

double UniformDistribution::sample(Engine &engine) const
{
	return dist(engine);
}

std::string UniformDistribution::info() const
{
	return "Uniform between " + std::to_string(dist.a()) + " and " +
		std::to_string(dist.b()) + ".";
}

} // namespace molstat
