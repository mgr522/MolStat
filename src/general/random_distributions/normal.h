/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file normal.h
 * \brief Interface for the normal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 */

#ifndef __normal_h__
#define __normal_h__

#include <memory>
#include <vector>
#include <string>
#include <random>
#include "rng.h"

namespace molstat {

/// Normal distribution.
class NormalDistribution : public RandomDistribution
{
protected:
	/// The C++11 normal distribution
	mutable std::normal_distribution<double> dist;

public:
	NormalDistribution() = delete;
	~NormalDistribution() = default;

	/**
	 * \brief Constructor specifying the mean and standard deviation.
	 *
	 * \throw std::invalid_argument if the standard deviation is non-positive.
	 *
	 * \param[in] mean The mean.
	 * \param[in] stdev The standard deviation.
	 */
	NormalDistribution(const double mean, const double stdev);

	virtual double sample(Engine &engine) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
