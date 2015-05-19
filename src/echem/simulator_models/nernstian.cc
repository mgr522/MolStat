/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Northwestern University. */

/**
 * \file nernstian.cc
 * \brief Simulator model for electron transfer with a Nernstian reaction.
 *
 * \author Bo Fu, Matthew G.\ Reuter
 * \date January 2015
 */

#include "nernstian.h"
#include <cmath>

namespace molstat {
namespace echem {

const std::size_t NernstianReaction::Index_Eref = 0;
const std::size_t NernstianReaction::Index_Af = 1;
const std::size_t NernstianReaction::Index_Ab = 2;

std::vector<std::string> NernstianReaction::get_names() const
{
	std::vector<std::string> ret(3);

	ret[Index_Eref] = "eref";
	ret[Index_Af] = "af";
	ret[Index_Ab] = "ab";

	return ret;
}

double NernstianReaction::ForwardETP(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &eref = params[Index_Eref];
	const double &af = params[Index_Af];
	const double &ab = params[Index_Ab];
	
	return eref - std::log(ab / af); 
}

double NernstianReaction::BackwardETP(const std::valarray<double> &params) const
{
	return ForwardETP(params);
}

} // namespace molstat::echem
} // namespace molstat
