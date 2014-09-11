/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file normal.h
 * \brief Interface for the normal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __normal_h__
#define __normal_h__

#include <memory>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include "rng.h"

using std::shared_ptr;

/**
 * \brief Normal distribution.
 */
class NormalDistribution : public RandomDistribution {
protected:
	/**
	 * \brief The mean.
	 */
	double mean;

	/**
	 * \brief The standard deviation.
	 */
	double stdev;

public:
	NormalDistribution() = delete;

	/**
	 * \brief Constructor specifying the mean and standard deviation.
	 *
	 * \param[in] mean_ The mean.
	 * \param[in] stdev_ The standard deviation.
	 */
	NormalDistribution(const double mean_, const double stdev_);

	/**
	 * \brief Destructor.
	 */
	~NormalDistribution() = default;

	/**
	 * \brief Samples from the uniform distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	virtual double sample(shared_ptr<gsl_rng> r) const override;

	/**
	 * \brief A description of this random number distribution.
	 *
	 * \return A string containing the description.
	 */
	virtual std::string info() const override;
};

#endif
