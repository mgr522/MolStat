/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model.h
 * \brief Defines the molstat::SimulateModel class for simulating histograms.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __simulate_model_h__
#define __simulate_model_h__

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <valarray>
#include <typeinfo>
#include <typeindex>
#include <general/random_distributions/rng.h>

namespace molstat {

// forward declaration
class SimulateModel;

/**
 * \brief The signature of a function that calculates an observable.
 */
using ObservableFunction =
	std::function<double(const std::valarray<double> &)>;

/**
 * \brief The signature of a function that produces an ObservableFunction,
 *    given the model.
 *
 * \throw molstat::IncompatibleObservable may be thrown if the observable
 *    and model are incompatible.
 */
using ObservableFactory =
	std::function<ObservableFunction(std::shared_ptr<const SimulateModel>)>;

/**
 * \brief Shortcut for the index type (alias for std::type_index) of an
 *    Observable.
 */
using ObservableIndex = std::type_index;

/**
 * \brief Base class for a model that uses model parameters to calculate
 *    observables.
 *
 * All models for simulating data should derive from this class. Derived
 * classes must implement functions that provide the number of model parameters
 * needed, as well as a name for each parameter. Derived classes should also
 * probably derive from molstat::Observable so that the MolStat simulator knows
 * that this model is compatible with that observable.
 *
 * Construction via the molstat::SimulateModelFactory interface ensures that
 * all required random number distributions are specified before using the
 * model to simulate data.
 */
class SimulateModel : public std::enable_shared_from_this<SimulateModel> {
protected:
	/**
	 * \brief Factories that produce an observable's function, assuming the
	 *    observable and model are compatible.
	 *
	 * The map is keyed by the molstat::ObservableIndex for the observable's
	 * class.
	 */
	std::map<ObservableIndex, ObservableFactory> compatible_observables;

	/**
	 * \brief Ordered vector of random number distributions for the various
	 *    model parameters.
	 */
	std::vector<std::shared_ptr<const RandomDistribution>> dists;

	/**
	 * \brief Gets a map of parameter name to index.
	 *
	 * \return The ordered list of names of distributions.
	 */
	virtual std::vector<std::string> get_names() const = 0;

	SimulateModel() = default;

public:
	virtual ~SimulateModel() = default;

	/**
	 * \brief Gets the number of model parameters for this model.
	 *
	 * \return The number of parameters used by the model for calculating
	 *    observables.
	 */
	virtual std::size_t get_num_parameters() const = 0;

	/**
	 * \brief Gets a function that calculates an observable, given a set of
	 *    model parameters.
	 *
	 * This default function verifies that the model and observable are
	 * compatible and then returns a function that casts `this` to the
	 * observable type and calls the proper function.
	 *
	 * \throw molstat::IncompatibleObservable if the desired observable is
	 *    incompatible with this model.
	 *
	 * \param[in] obs The type_index of the class for the observable.
	 * \return A function that calculates the observable.
	 */
	virtual ObservableFunction getObservableFunction(
		const ObservableIndex &obs) const;

	/**
	 * \brief Generates a set of model parameters using the specified random
	 *    distributions.
	 *
	 * \param[in] r Handle to the GSL random number generator.
	 * \return A set of model parameters.
	 */
	virtual std::valarray<double> generateParameters(gsl_rng_ptr &r) const;

	// the factory needs to get at the internal details
	friend class SimulateModelFactory;
};

#if 0
/**
 * \brief Shortcut for a function that constructs a molstat::SimulateModel.
 */
using SimulateModelFactory =
	std::function<std::shared_ptr<SimulateModel>()>;

/**
 * \brief Gets a molstat::SimulateModelFactory for the given model type.
 *
 * \tparam T The type of molstat::SimulateModel to construct.
 * \return A function that creates the model when invoked.
 */
template<typename T>
inline SimulateModelFactory GetSimulateModelFactory() {
	return [] () -> std::shared_ptr<SimulateModel> {
		return std::make_shared<T>();
	};
}
#endif

template<typename T>
inline ObservableIndex GetObservableIndex() {
	return std::type_index{ typeid(T) };
}

} // namespace molstat

#endif