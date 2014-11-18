/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file constant.cc
 * \brief Implementation of the constant distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 * \endinternal
 */

#include "constant.h"

namespace molstat {

ConstantDistribution::ConstantDistribution(const double val)
	: RandomDistribution(), value(val)
{
}

double ConstantDistribution::sample(Engine &engine) const
{
	return value;
}

std::string ConstantDistribution::info() const
{
	return "Constant = " + std::to_string(value) + ".";
}

} // namespace molstat
