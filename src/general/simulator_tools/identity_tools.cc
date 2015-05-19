/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file identity_tools.cc
 * \brief Implementations of an observable and model used to test the
 *    simulator.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "identity_tools.h"

namespace molstat {

std::vector<std::string> IdentityModel::get_names() const
{
	return { "parameter" };
}

double IdentityModel::identity(const std::valarray<double> &params) const
{
	return params[0];
}

} // namespace molstat
