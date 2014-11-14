/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file lognormal.h
 * \brief Interface for the lognormal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __lognormal_h__
#define __lognormal_h__

#include <memory>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include "rng.h"

namespace molstat {

/// Lognormal distribution.
class LognormalDistribution : public RandomDistribution
{
protected:
	/// The average, in log-space.
	double zeta;

	/// The standard deviation, in log-space.
	double sigma;

public:
	LognormalDistribution() = delete;
	~LognormalDistribution() = default;

	/**
	 * \brief Constructor specifying the average and standard deviation.
	 *
	 * \param[in] zeta_ The average (in log-space).
	 * \param[in] sigma_ The standard deviation (in log-space).
	 */
	LognormalDistribution(const double zeta_, const double sigma_);

	virtual double sample(gsl_rng_ptr &r) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
