/**
 * \file gamma.h
 * \brief Interface for the gamma distribution.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __gamma_h__
#define __gamma_h__

#include <memory>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include "rng.h"

using std::shared_ptr;

/**
 * \brief Gamma distribution.
 */
class GammaDistribution : public RandomDistribution {
protected:
	/**
	 * \brief The shape factor.
	 */
	double shape;

	/**
	 * \brief The scale factor.
	 */
	double scale;

public:
	GammaDistribution() = delete;

	/**
	 * \brief Constructor specifying the shape and scale factors.
	 *
	 * \param[in] shape_ The shape factor.
	 * \param[in] scale_ The scale factor.
	 */
	GammaDistribution(const double shape_, const double scale_);

	/**
	 * \brief Destructor.
	 */
	~GammaDistribution() = default;

	/**
	 * \brief Samples from the lognormal distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	virtual double sample(shared_ptr<gsl_rng> r) const;
};

#endif
