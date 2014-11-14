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

namespace molstat {

/// Uniform distribution.
class UniformDistribution : public RandomDistribution
{
protected:
	/// The lower bound of possible values.
	double lower;

	/// The upper bound of possible values.
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

	virtual double sample(gsl_rng_ptr &r) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif