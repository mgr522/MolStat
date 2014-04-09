/**
 * \file rng.h
 * \brief Interface for random number generation.
 *
 * Provides a common interface for random number generation with GSL.
 * The following distributions are implemented:
 *    - Constant distribution (the value is fixed).
 *    - Uniform distribution.
 *    - Normal distribution.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include <memory>
#include <gsl/gsl_rng.h>

using std::shared_ptr;

/**
 * \brief Interface for random number generation.
 */
class RandomDistribution {
public:
	RandomDistribution() = default;

	/**
	 * \brief Destructor.
	 */
	virtual ~RandomDistribution() = default;

	/**
	 * \brief Samples from the random number distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	virtual double sample(shared_ptr<gsl_rng> r) const = 0;
};

/**
 * \brief Constant distribution.
 *
 * This distribution is included for the case where the value should be fixed.
 */
class ConstantDistribution : public RandomDistribution {
protected:
	/**
	 * \brief The fixed value to return.
	 */
	double value;

public:
	ConstantDistribution() = delete;

	/**
	 * \brief Constructor specifying the fixed value.
	 *
	 * \param[in] val The fixed value.
	 */
	ConstantDistribution(const double val);

	/**
	 * \brief Samples from the constant distribution (returns value).
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	double sample(shared_ptr<gsl_rng> r) const;
};

/**
 * \brief Uniform distribution.
 */
class UniformDistribution : public RandomDistribution {
protected:
	/**
	 * \brief The lower bound of possible values.
	 */
	double lower;

	/**
	 * \brief The upper bound of possible values.
	 */
	double upper;

public:
	UniformDistribution() = delete;

	/**
	 * \brief Constructor specifying the range.
	 *
	 * \param[in] low The lower bound.
	 * \param[in] up The upper bound.
	 */
	UniformDistribution(const double low, const double up);

	/**
	 * \brief Samples from the uniform distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	double sample(shared_ptr<gsl_rng> r) const;
};

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
	 * \brief Samples from the uniform distribution.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The random number.
	 */
	double sample(shared_ptr<gsl_rng> r) const;
};
