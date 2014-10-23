/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
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

using namespace std;

namespace molstat {

std::unique_ptr<BinStyle> BinStyleFactory(
	TokenContainer &&tokens) {

	unique_ptr<BinStyle> ret;

	if(tokens.size() == 0)
		throw invalid_argument("Empty line.");

	// the first token is the name of the binning style
	string name = to_lower(tokens.front());
	tokens.pop();

	if(name == "linear") {
		ret.reset(new BinLinear());
	}
	else if(name == "log") {
		// need to read the base, if available. If not, use 10.
		double b;

		if(tokens.size() > 0) {
			b = cast_string<double>(tokens.front());
			if(b <= 0.)
				throw invalid_argument("The logarithm base must be positive.");
		}
		else
			b = 10.;

		ret.reset(new BinLog(b));
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
