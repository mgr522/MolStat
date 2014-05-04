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

using std::shared_ptr;

/**
 * \brief Interface for random number generation.
 */
class RandomDistribution {
public:
	RandomDistribution() = default;

	/**
	 * \brief Destructor.
	 */
	virtual ~RandomDistribution() = default;

	/**
	 * \brief Samples from the random number distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	virtual double sample(shared_ptr<gsl_rng> r) const = 0;
};

/**
 * \brief Gets parameters from a tokenized string and forms a
 *    RandomDistribution.
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
shared_ptr<RandomDistribution> distribution_from_tokens(
	const std::vector<std::string> &tokens);

#endif
