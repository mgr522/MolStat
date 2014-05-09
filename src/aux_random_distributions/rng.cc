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
#include "../string_tools.h"
#include "constant.h"
#include "uniform.h"
#include "normal.h"
#include "lognormal.h"

using namespace std;

/**
 * \brief Turns a string into a double.
 *
 * \throw invalid_argument if the string cannot be cast to a double.
 *
 * \param[in] str The string to be cast.
 * \return The double.
 */
static double read_double(const string &str) {
	double ret;

	try {
		ret = stod(str);
	}
	catch(const invalid_argument &e) {
		throw invalid_argument("Unable to convert \"" + str + "\" to a double.");
	}

	return ret;
}

shared_ptr<RandomDistribution> distribution_from_tokens(
	const std::vector<std::string> &tokens) {

	shared_ptr<RandomDistribution> ret;

	if(tokens.size() < 1)
		throw invalid_argument("Empty line.");

	// the first token is the type of distribution to form
	string type = tokens[0];
	make_lower(type);

	if(type == "constant") {
		double val;

		// tokens[1] is the value to return; cast it to a double
		if(tokens.size() < 2)
			throw invalid_argument("Invalid constant distribution. Use\n" \
				"   constant value\n" \
				"where value is the value to be returned.");

		val = read_double(tokens[1]);

		ret = make_shared<ConstantDistribution>(val);
	}
	else if(type == "uniform") {
		double lower, upper;

		// tokens[1] and tokens[2] are the lower and upper bounds of the range,
		// respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid uniform distribution. Use\n" \
				"   uniform lower upper\n" \
				"where lower and upper are the bounds, respectively.");

		lower = read_double(tokens[1]);
		upper = read_double(tokens[2]);

		if(lower >= upper)
			throw invalid_argument("Uniform Distribution: The lower bound must " \
				"be lower than the upper bound.");

		ret = make_shared<UniformDistribution>(lower, upper);
	}
	else if(type == "normal" || type == "gaussian") {
		double mean, stdev;

		// tokens[1] and tokens[2] are the mean and standard deviation of the
		// distribution, respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid normal distribution. Use\n" \
				"   normal mean standard-deviation");

		mean = read_double(tokens[1]);
		stdev = read_double(tokens[2]);

		if(stdev <= 0.)
			throw invalid_argument("Normal Distribution: The standard deviation" \
				" must be positive.");

		ret = make_shared<NormalDistribution>(mean, stdev);
	}
	else if(type == "lognormal") {
		double zeta, sigma;

		// tokens[1] and tokens[2] are the mean and standard deviation (in log-
		// space), respectively
		if(tokens.size() < 3)
			throw invalid_argument("Invalid lognormal distribution. Use\n" \
				"   lognormal zeta sigma");

		zeta = read_double(tokens[1]);
		sigma = read_double(tokens[2]);

		if(sigma <= 0.)
			throw invalid_argument("Lognormal Distribution: The standard " \
				"deviation (sigma) must be positive.");

		ret = make_shared<LognormalDistribution>(zeta, sigma);
	}
	else
		throw invalid_argument("Unrecognized probability distribution.\n" \
			"Possible options are:\n" \
			"   Constant - Specify a value.\n" \
			"   Uniform - Uniform distribution.\n" \
			"   Normal - Normal (Gaussian) distribution.\n" \
			"   Gaussian - Normal (Gaussian) distribution.\n");

	return ret;
}
