/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file rectangular_barrier.cc
 * \brief Channel with transmission given by a rectangular barrier.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include "rectangular_barrier.h"
#include <cmath>

namespace molstat {
namespace transport {

const std::size_t RectangularBarrier::Index_EF = TransportJunction::Index_EF;
const std::size_t RectangularBarrier::Index_V = TransportJunction::Index_V;
const std::size_t RectangularBarrier::Index_h = 2;
const std::size_t RectangularBarrier::Index_w = 3;

std::vector<std::string> RectangularBarrier::get_names() const
{
	std::vector<std::string> ret(2);

	// subtract two because we expect 2 parameters from the TransportJunction
	ret[Index_h - 2] = "height";
	ret[Index_w - 2] = "width";

	return ret;
}

double RectangularBarrier::transmission(const double e, const double h,
	const double w)
{
	// sqrt(2m (eV)) / hbar = 5.12317 / nm
	const double sinhval = std::sinh(5.12317 * w) * h;

	return 1. / (1. + 0.25 * sinhval * sinhval / (e * (h - e)));
}

double RectangularBarrier::ZeroBiasG(const std::valarray<double> &params) const
{
	// unpack the parameters
	const double &ef = params[Index_EF];
	const double &h = params[Index_h];
	const double &w = params[Index_w];
	
	return transmission(ef, h, w);
}

} // namespace molstat::transport
} // namespace molstat