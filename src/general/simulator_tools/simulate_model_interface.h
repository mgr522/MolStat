/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface.h
 * \brief Defines an abstract class encapsulating a model for simulating
 *    histograms.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __simulate_model_interface_h__
#define __simulate_model_interface_h__

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <general/random_distributions/rng.h>

using std::shared_ptr;

/**
 * \brief Base class encapsulating a model for simulating data to construct a
 *    histogram.
 *
 * This base class stores the random distributions for the variables as a map
 * from variable name to distribution.
 */
class SimulateModel {
private:
	/**
	 * \internal
	 * \brief The random distributions, ordered as specified during
	 *    construction.
	 * \endinternal
	 */
	std::vector<shared_ptr<const RandomDistribution>> dists;

public:
	SimulateModel() = delete;

	/**
	 * \brief Constructor requiring a list of available distributions and an
	 *    ordered list of needed distributions.
	 *
	 * \throw runtime_error If there a required distribution is not found among
	 *    the available distributions.
	 *
	 * \param[in] avail The available distributions, keyed by name,
	 * \param[in] names The names of required distributions, in a particular
	 *    order.
	 */
	SimulateModel(
		const std::map<std::string, shared_ptr<RandomDistribution>> &avail,
		const std::vector<std::string> &names);

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
	virtual ~SimulateModel() = default;

	/**
	 * \brief Samples from the random distributions.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[out] vals The random numbers, in the order specified during
	 *    construction.
	 */
	void sample(shared_ptr<gsl_rng> r, std::vector<double> &vals) const;
};

#endif
