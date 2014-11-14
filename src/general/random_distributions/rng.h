/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file rng.h
 * \brief Interface for random number generation.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __rng_h__
#define __rng_h__

#include <memory>
#include <queue>
#include <string>
#include <gsl/gsl_rng.h>
#include <general/string_tools.h>

namespace molstat {

/**
 * \brief Shortcut for a handle for GSL random number generation.
 *
 * Need to specify the deleter for unique pointer; that is, the gsl_rng_free
 * function.
 */
using gsl_rng_ptr =
	std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)>;

/// Interface for random number generation.
class RandomDistribution
{
public:
	RandomDistribution() = default;
	virtual ~RandomDistribution() = default;

	/**
	 * \brief Samples from the random number distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	virtual double sample(gsl_rng_ptr &r) const = 0;

	/**
	 * \brief A description of this random number distribution.
	 *
	 * \return A string containing the description.
	 */
	virtual std::string info() const = 0;
};

/**
 * \brief Factory for random number distributions.
 *
 * Gets parameters from a tokenized string and forms a RandomDistribution.
 *
 * The first token is the name of the distribution. The remaining tokens are
 * the parameters for the distribution.
 *
 * The list of tokens is destroyed by this function.
 *
 * \exception std::invalid_argument if the parameters are invalid for the
 *    specified distribution.
 *
 * \param[in] tokens The tokens.
 * \return The RandomDistribution.
 */
std::unique_ptr<RandomDistribution> RandomDistributionFactory(
	TokenContainer &&tokens);

} // namespace molstat

#endif