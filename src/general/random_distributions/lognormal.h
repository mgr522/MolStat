/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file lognormal.h
 * \brief Interface for the lognormal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 */

#ifndef __lognormal_h__
#define __lognormal_h__

#include <memory>
#include <vector>
#include <string>
#include <random>
#include "rng.h"

namespace molstat {

/// Lognormal distribution.
class LognormalDistribution : public RandomDistribution
{
protected:
	/// The C++11 lognormal distribution.
	mutable std::lognormal_distribution<double> dist;

public:
	LognormalDistribution() = delete;
	~LognormalDistribution() = default;

	/**
	 * \brief Constructor specifying the average and standard deviation.
	 *
	 * \throw std::invalid_argument if the standard deviation is non-positive.
	 *
	 * \param[in] zeta The average (in log-space).
	 * \param[in] sigma The standard deviation (in log-space).
	 */
	LognormalDistribution(const double zeta, const double sigma);

	virtual double sample(Engine &engine) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
