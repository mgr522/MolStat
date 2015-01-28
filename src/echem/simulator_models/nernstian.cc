/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file nernstian.cc
 * \brief Simulator model for electron transfer with a Nernstian reaction.
 *
 * \author Bo Fu
 * \date December 2014
 */

#include "nernstian.h"
#include <cmath>

#define kB 1.38066e-23 // Unit J/K
#define e_charge 1.602189e-19 
#define T 300.0 //Temperature
#define n 1.0 // The number of electrons involved in the Nernstian reaction.

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


double NernstianReaction::RedoxETP(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &eref = params[Index_Eref];
	const double &af = params[Index_Af];
	const double &ab = params[Index_Ab];
	
	return eref - (kB * T / (n * e_charge )) * std::log(ab / af); 
}


} // namespace molstat::echem
} // namespace molstat
