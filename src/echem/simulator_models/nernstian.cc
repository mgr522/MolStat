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

#define kB 1.38066e-23
#define e_charge 1.602189e-19 

namespace molstat {
namespace echem {

const std::size_t SingleMoleculeEchemNernstian::Index_E0 = 0;
const std::size_t SingleMoleculeEchemNernstian::Index_Af = 1;
const std::size_t SingleMoleculeEchemNernstian::Index_Ab = 2;
const std::size_t SingleMoleculeEchemNernstian::Index_T = 3;
const std::size_t SingleMoleculeEchemNernstian::Index_n = 4;

std::vector<std::string> SingleMoleculeEchemNernstian::get_names() const
{
	std::vector<std::string> ret(5);

	ret[Index_E0] = "e0";
	ret[Index_Af] = "af";
	ret[Index_Ab] = "ab";
	ret[Index_T] = "temp";
	ret[Index_n] = "n";

	return ret;
}


double SingleMoleculeEchemNernstian::PeakV(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &e0 = params[Index_E0];
	const double &af = params[Index_Af];
	const double &ab = params[Index_Ab];
	const double &T = params[Index_T];
	const double &n = params[Index_n];
	
	return e0 - (kB * T / (n * e_charge )) * std::log(ab / af); 
}


} // namespace molstat::echem
} // namespace molstat
