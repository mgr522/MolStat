/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file constant.h
 * \brief Interface for the constant distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __constant_h__
#define __constant_h__

#include <memory>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include "rng.h"

/**
 * \brief Constant distribution.
 *
 * This distribution is included for the case where the value should be fixed.
 */
class ConstantDistribution : public RandomDistribution {
protected:
	/**
	 * \brief The fixed value to return.
	 */
	double value;

public:
	ConstantDistribution() = delete;
	~ConstantDistribution() = default;

	/**
	 * \brief Constructor specifying the fixed value.
	 *
	 * \param[in] val The fixed value.
	 */
	ConstantDistribution(const double val);

	/**
	 * \brief Samples from the constant distribution (returns value).
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