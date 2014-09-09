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
#include <stdexcept>
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
	virtual ~SimulateModel() = default;

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
inline SimulateModelInstantiator SimulateModelAdd() {
	return []
		(const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
		-> shared_ptr<SimulateModel> {

		return std::make_shared<T>(avail);
	};
}

/**
 * \brief Shortcut for the function signature of a observable.
 *
 * Observables take in a GSL random number generator handle and produce the
 * value of the observable, as determined by the model.
 *
 * \tparam N The number of variables in the observable.
 */
template<size_t N>
using Observable = std::function<std::array<double, N>(shared_ptr<gsl_rng>)>;

/**
 * \brief Returns a function that checks compatibility of a model with an
 *    observable.
 *
 * The returned function first checks that a given model is of a class that
 * implements the desired observable. If not, it throws an exception. If it
 * is, it then creates a wrapper to the observable's function that can be
 * called to simulate data.
 *
 * \throw runtime_error If the model and observable are incompatible; that is,
 *    the model does not implement the observable.
 *
 * \tparam N The dimensionality of the observable.
 * \tparam T The class name of the observable's interface.
 * \param[in] memfunc Pointer to the member function of the observable
 *    interface that simulates the observable.
 * \return The described function.
 */
template<size_t N, typename T>
std::function<Observable<N>(const shared_ptr<SimulateModel>)> ObservableCheck(
	std::array<double, N> (T::*memfunc)(shared_ptr<gsl_rng>) const) {

	return [=] (const shared_ptr<SimulateModel> model) {
		shared_ptr<T> cast = std::dynamic_pointer_cast<T>(model);

		if(cast == nullptr)
			throw std::runtime_error("Incompatible model and observable.");

		return std::bind(memfunc, cast, std::placeholders::_1);
	};
}

/**
 * \brief Retypes a function that returns array<double, 1> to return
 *    array<double, 2>.
 *
 * This simply adds a zero to the 0th element; the 1D element becomes element
 * [1] in the 2-array.
 *
 * \param[in] f The 1D observable function to be wrapped.
 * \return The wrapped function.
 */
std::function<Observable<2>(const shared_ptr<SimulateModel>)> Obs2(
	const std::function<Observable<1>(const shared_ptr<SimulateModel>)> &f);

#endif
