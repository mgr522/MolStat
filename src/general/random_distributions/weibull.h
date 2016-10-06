/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2016 Stony Brook University. */

/**
 * \file weibull.h
 * \brief Interface for the Weibull distribution.
 *
 * \author Matthew G.\ Reuter
 * \date October 2016
 */

#ifndef __weibull_h__
#define __weibull_h__

#include <memory>
#include <vector>
#include <string>
#include <random>
#include "rng.h"

namespace molstat {

/// Weibull distribution.
class WeibullDistribution : public RandomDistribution
{
protected:
	/// The C++11 Weibull distribution.
	mutable std::weibull_distribution<double> dist;

public:
	WeibullDistribution() = delete;
	~WeibullDistribution() = default;

	/**
	 * \brief Constructor specifying the shape and scale factors.
	 *
	 * \throw std::invalid_argument if either parameter is non-positive.
	 *
	 * \param[in] shape The shape factor.
	 * \param[in] scale The scale factor.
	 */
	WeibullDistribution(const double shape, const double scale);

	virtual double sample(Engine &engine) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
