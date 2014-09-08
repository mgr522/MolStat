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

/**
 * \brief Shortcut for the function signature of an "instantiator" for a
 *    SimulateModel in the simulator.
 *
 * SimulateModel objects are created by passing in a map of available
 * RandomDistributions; the specification of a model should supply the names
 * of required parameters. This function type produces the SimulateModel
 * from the map of RandomDistributions.
 */
typedef std::function<shared_ptr<SimulateModel>
	(const std::map<std::string, shared_ptr<RandomDistribution>> &)>
	SimulateModelInstantiator;

/**
 * \brief Creates a SimulateModelInstantiator for a particular model.
 *
 * \tparam T The type of SimulateModel we wish to instantiate.
 * \return A function for instantiating the class from a map of available
 *    random number distributions.
 */
template<typename T>
inline SimulateModelInstantiator SimulateModelInstance() {
	return []
		(const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
		-> shared_ptr<SimulateModel> {

		return std::make_shared<T>(avail);
	};
}

#endif
