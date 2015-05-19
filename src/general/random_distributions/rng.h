/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file rng.h
 * \brief Interface for random number generation.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 */

#ifndef __rng_h__
#define __rng_h__

#include <memory>
#include <queue>
#include <string>
#include <random>
#include <general/string_tools.h>

namespace molstat {

/// The type of C++11 random number engine to use.
using Engine = std::default_random_engine;

/// Interface for random number generation.
class RandomDistribution
{
public:
	RandomDistribution() = default;
	virtual ~RandomDistribution() = default;

	/**
	 * \brief Samples from the random number distribution.
	 *
	 * \param[in] engine The random number engine.
	 * \return The random number.
	 */
	virtual double sample(Engine &engine) const = 0;

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
 * \return The molstat::RandomDistribution.
 */
std::unique_ptr<RandomDistribution> RandomDistributionFactory(
	TokenContainer &&tokens);

} // namespace molstat

#endif
