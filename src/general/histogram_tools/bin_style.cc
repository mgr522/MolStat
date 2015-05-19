/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file bin_style.cc
 * \brief Implements a function to create a BinStyle from a vector of tokens.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 */

#include "bin_style.h"
#include "bin_linear.h"
#include "bin_log.h"
#include <general/string_tools.h>

using namespace std;

namespace molstat
{

BinStyle::BinStyle(const std::size_t nbins_)
	: nbins(nbins_)
{
}

std::unique_ptr<BinStyle> BinStyleFactory(TokenContainer &&tokens)
{
	unique_ptr<BinStyle> ret;

	if(tokens.size() == 0)
		throw invalid_argument("Empty line.");

	// the first token is the number of bins to use
	size_t nbins;
	try
	{
		nbins = cast_string<size_t>(tokens.front());
	}
	catch(const bad_cast &e)
	{
		throw invalid_argument("Unable to determine the number of bins.");
	}
	tokens.pop();

	// the next token is the name of the binning style
	if(tokens.size() == 0)
		throw invalid_argument("No binning style specified.");

	string name = to_lower(tokens.front());
	tokens.pop();

	if(name == "linear")
	{
		ret.reset(new BinLinear(nbins));
	}
	else if(name == "log")
	{
		// need to read the base, if available. If not, use 10.
		double b;

		if(tokens.size() > 0)
		{
			try
			{
				b = cast_string<double>(tokens.front());
			}
			catch(const bad_cast &e)
			{
				throw invalid_argument(
					"Unable to convert the base to a numerical value.");
			}

			if(b <= 0.)
				throw invalid_argument("The logarithm base must be positive.");
		}
		else
			b = 10.;

		ret.reset(new BinLog(nbins, b));
	}
	else
		throw invalid_argument(
			"Unrecognized binning style: \"" + name + "\".\n" \
			"Possible options are:\n" \
			"   Linear - Linear binning.\n" \
			"   Log - Logarithmic binning (base defaults to 10).\n");

	return ret;
}

} // namespace molstat
