/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file rng.cc
 * \brief Implementation of the function for generating a random number
 *    distribution from a tokenized line of input.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "rng.h"
#include <gsl/gsl_randist.h>
#include "constant.h"
#include "uniform.h"
#include "normal.h"
#include "lognormal.h"
#include "gamma.h"

using namespace std;

namespace molstat {

std::unique_ptr<RandomDistribution> RandomDistributionFactory(
	TokenContainer &&tokens)
{
	unique_ptr<RandomDistribution> ret;

	if(tokens.size() == 0)
		throw invalid_argument("Empty line.");

	// the first token is the type of distribution to form
	string type = to_lower(tokens.front());
	tokens.pop();

	if(type == "constant")
	{
		// tokens[0] is the value to return; cast it to a double
		if(tokens.size() < 1)
			throw invalid_argument("Invalid constant distribution. Use\n" \
				"   constant value\n" \
				"where value is the value to be returned.");

		double val;
		try
		{
			val = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"value\" to a numeric value.");
		}

		ret = unique_ptr<RandomDistribution>(new ConstantDistribution(val));
	}
	else if(type == "uniform")
	{
		// tokens[0] and tokens[1] are the lower and upper bounds of the range,
		// respectively
		if(tokens.size() < 2)
			throw invalid_argument("Invalid uniform distribution. Use\n" \
				"   uniform lower upper\n" \
				"where lower and upper are the bounds, respectively.");

		double lower, upper;
		try
		{
			lower = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"lower\" to a numeric value.");
		}
		tokens.pop();

		try
		{
			upper = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"upper\" to a numeric value.");
		}

		if(lower >= upper)
			throw invalid_argument("Uniform Distribution: The lower bound must " \
				"be lower than the upper bound.");

		ret = unique_ptr<RandomDistribution>(
			new UniformDistribution(lower, upper));
	}
	else if(type == "normal" || type == "gaussian")
	{
		// tokens[0] and tokens[1] are the mean and standard deviation of the
		// distribution, respectively
		if(tokens.size() < 2)
			throw invalid_argument("Invalid normal distribution. Use\n" \
				"   normal mean standard-deviation");

		double mean, stdev;
		try
		{
			mean = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"mean\" to a numeric value.");
		}
		tokens.pop();

		try
		{
			stdev = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"stdev\" to a numeric value.");
		}

		if(stdev <= 0.)
			throw invalid_argument("Normal Distribution: The standard deviation" \
				" must be positive.");

		ret = unique_ptr<RandomDistribution>(
			new NormalDistribution(mean, stdev));
	}
	else if(type == "lognormal")
	{
		// tokens[0] and tokens[1] are the mean and standard deviation (in log-
		// space), respectively
		if(tokens.size() < 2)
			throw invalid_argument("Invalid lognormal distribution. Use\n" \
				"   lognormal zeta sigma");

		double zeta, sigma;
		try
		{
			zeta = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"zeta\" to a numeric value.");
		}
		tokens.pop();

		try
		{
			sigma = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"sigma\" to a numeric value.");
		}

		if(sigma <= 0.)
			throw invalid_argument("Lognormal Distribution: The standard " \
				"deviation (sigma) must be positive.");

		ret = unique_ptr<RandomDistribution>(
			new LognormalDistribution(zeta, sigma));
	}
	else if(type == "gamma")
	{
		// tokens[0] and tokens[1] are the shape and scale factors, respectively
		if(tokens.size() < 2)
			throw invalid_argument("Invalid gamma distribution. Use\n" \
				"   gamma shape scale");

		double shape, scale;
		try
		{
			shape = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"shape\" to a numeric value.");
		}
		tokens.pop();

		try
		{
			scale = cast_string<double>(tokens.front());
		}
		catch(const bad_cast &e)
		{
			throw invalid_argument(
				"Unable to convert \"scale\" to a numeric value.");
		}

		if(shape <= 0. || scale <= 0.)
			throw invalid_argument("Gamma Distribution: The shape and scale " \
				"factors must be positive.");

		ret = unique_ptr<RandomDistribution>(
			new GammaDistribution(shape, scale));
	}
	else
		throw invalid_argument(
			string("Unrecognized probability distribution: \"") + type + "\".\n" \
			"Possible options are:\n" \
			"   Constant - Specify a value.\n" \
			"   Uniform - Uniform distribution.\n" \
			"   Normal - Normal (Gaussian) distribution.\n" \
			"   Gaussian - Normal (Gaussian) distribution.\n" \
			"   Lognormal - Lognormal distribution.\n" \
			"   Gamma - Gamma distribution.\n");

	return ret;
}

} // namespace molstat