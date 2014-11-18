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
#include <random>
#include "rng.h"

namespace molstat {

/**
 * \brief Constant distribution.
 *
 * This distribution is included for the case where the value should be fixed.
 */
class ConstantDistribution : public RandomDistribution
{
protected:
	/// The fixed value to return.
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

	virtual double sample(Engine &engine) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif