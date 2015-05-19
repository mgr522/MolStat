/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file uniform.h
 * \brief Interface for the uniform distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 */

#ifndef __uniform_h__
#define __uniform_h__

#include <memory>
#include <vector>
#include <string>
#include <random>
#include "rng.h"

namespace molstat {

/// Uniform distribution.
class UniformDistribution : public RandomDistribution
{
protected:
	/// The C++11 uniform distribution
	mutable std::uniform_real_distribution<double> dist;
	
public:
	UniformDistribution() = delete;
	~UniformDistribution() = default;

	/**
	 * \brief Constructor specifying the range.
	 *
	 * \throw std::invalid_argument if the lower bound is not less than the
	 *    upper bound.
	 *
	 * \param[in] low The lower bound.
	 * \param[in] up The upper bound.
	 */
	UniformDistribution(const double low, const double up);

	virtual double sample(Engine &engine) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
