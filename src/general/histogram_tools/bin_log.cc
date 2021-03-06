/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file bin_log.cc
 * \brief Implements logarithmic binning.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 */

#include "bin_log.h"
#include <cmath>

namespace molstat
{

BinLog::BinLog(const std::size_t nbin_, const double b_)
	: BinStyle(nbin_), b(b_)
{
}

double BinLog::mask(const double x) const
{
	return log(x) / log(b);
}

double BinLog::invmask(const double u) const
{
	return pow(b, u);
}

double BinLog::dmaskdx(const double x) const
{
	return 1. / (x * log(b));
}

std::string BinLog::info() const
{
	return std::to_string(nbins) + " logarithmic bins, base " +
		std::to_string(b);
}

} // namespace molstat
