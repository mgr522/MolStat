/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file gamma.h
 * \brief Interface for the gamma distribution.
 *
 * \author Matthew G.\ Reuter
 * \date November 2014
 */

#ifndef __gamma_h__
#define __gamma_h__

#include <memory>
#include <vector>
#include <string>
#include <random>
#include "rng.h"

namespace molstat {

/// Gamma distribution.
class GammaDistribution : public RandomDistribution
{
protected:
	/// The C++11 gamma distribution
	mutable std::gamma_distribution<double> dist;

public:
	GammaDistribution() = delete;
	~GammaDistribution() = default;

	/**
	 * \brief Constructor specifying the shape and scale factors.
	 *
	 * \throw std::invalid_argument if either the shape or scale factor is
	 *    non-negative.
	 *
	 * \param[in] shape The shape factor.
	 * \param[in] scale The scale factor.
	 */
	GammaDistribution(const double shape, const double scale);

	virtual double sample(Engine &engine) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
