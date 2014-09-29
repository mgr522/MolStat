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
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>

/**
 * \brief Shortcut for a handle for GSL random number generation.
 *
 * Needed to specify the deleter for unique pointer
 */
using gsl_rng_ptr =
	std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)>;

/**
 * \brief Interface for random number generation.
 */
class RandomDistribution {
public:
	RandomDistribution() = default;

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
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
 * \exception std::invalid_argument if the parameters are invalid for the
 *    specified distribution.
 *
 * \param[in] tokens The tokens.
 * \return The RandomDistribution.
 */
std::unique_ptr<RandomDistribution> RandomDistributionFactory(
	const std::vector<std::string> &tokens);

#endif