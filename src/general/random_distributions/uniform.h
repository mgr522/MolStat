/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file uniform.h
 * \brief Interface for the uniform distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __uniform_h__
#define __uniform_h__

#include <memory>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include "rng.h"

/**
 * \brief Uniform distribution.
 */
class UniformDistribution : public RandomDistribution {
protected:
	/**
	 * \brief The lower bound of possible values.
	 */
	double lower;

	/**
	 * \brief The upper bound of possible values.
	 */
	double upper;

public:
	UniformDistribution() = delete;
	~UniformDistribution() = default;

	/**
	 * \brief Constructor specifying the range.
	 *
	 * \param[in] low The lower bound.
	 * \param[in] up The upper bound.
	 */
	UniformDistribution(const double low, const double up);

	/**
	 * \brief Samples from the uniform distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	virtual double sample(gsl_rng_ptr &r) const override;

	/**
	 * \brief A description of this random number distribution.
	 *
	 * \return A string containing the description.
	 */
	virtual std::string info() const override;
};

#endif