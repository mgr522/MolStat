/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model.h
 * \brief Defines the molstat::SimulateModel class for simulating histograms,
 *    the molstat::SimulateModelFactory class for creating models at runtime,
 *    and other auxiliary types/aliases for the simulator interface.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __simulate_model_h__
#define __simulate_model_h__

#include <memory>
#include <functional>
#include <valarray>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <typeinfo>
#include <typeindex>
#include <general/random_distributions/rng.h>
#include <general/string_tools.h>

namespace molstat {

// forward declarations
class SimulateModel;
class RandomDistribution;

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

/**
 * \brief Factory class for creating molstat::SimulateModels at runtime.
 *
 * This class instantiates models at runtime, making sure the model is in a
 * valid state before giving the user access to it.
 */
class SimulateModelFactory {
private:
	SimulateModelFactory() = default;

	/**
	 * \brief Pointer to the model being constructed.
	 */
	std::shared_ptr<SimulateModel> model;

	/**
	 * \brief The set of distribution names for the model that still need to
	 *    be specified.
	 */
	std::set<std::string> remaining_names;

	/**
	 * \brief Cache of the names of distributions required by the model.
	 */
	std::vector<std::string> model_names;

public:
	SimulateModelFactory(SimulateModelFactory &&) = default;
	SimulateModelFactory &operator=(SimulateModelFactory &&) = default;

	/**
	 * \brief Creates a molstat::SimulateModelFactory with the underlying model
	 *    of type T.
	 *
	 * \tparam T The type of molstat::SimulateModel to build with this factory.
	 * \return The factory.
	 */
	template<typename T>
	static SimulateModelFactory makeFactory();

	/**
	 * \brief Adds a random distribution to the model.
	 *
	 * \param[in] name The name of the distribution being added.
	 * \param[in] dist The distribution being added.
	 * \return The factory.
	 */
	 SimulateModelFactory &setDistribution(std::string name,
	 	std::shared_ptr<const RandomDistribution> dist);

	 /**
	  * \brief Returns the constructed model.
	  *
	  * Some runtime error checking, such as making sure all distributions are
	  * specified is performed here.
	  *
	  * \throw molstat::MissingDistribution if one of the required distributions
	  *    has not been specified.
	  *
	  * \return Pointer to the constructed model.
	  */
	 std::shared_ptr<SimulateModel> getModel();
};

/**
 * \brief Gets a function that produces a molstat::ObservableIndex for the
 *    given observable type.
 *
 * \tparam T The type of molstat::Observable to use.
 * \return A function that gives the observable's molstat::ObservableIndex.
 */
template<typename T>
inline ObservableIndex GetObservableIndex() {
	return std::type_index{ typeid(T) };
}

/**
 * \brief Shortcut for a function that constructs a
 *    molstat::SimulateModelFactory.
 */
using SimulateModelFactoryFactory =
	std::function<SimulateModelFactory()>;

/**
 * \brief Gets a function that produces a molstat::SimulateModelFactory
 *    for the given model type.
 *
 * \tparam T The type of molstat::SimulateModel to construct.
 * \return A function that creates the model when invoked.
 */
template<typename T>
inline SimulateModelFactoryFactory GetSimulateModelFactory() {
	return [] () -> SimulateModelFactory {
		return SimulateModelFactory::makeFactory<T>();
	};
}

// templated definitions
template<typename T>
SimulateModelFactory SimulateModelFactory::makeFactory() {
	using namespace std;

	SimulateModelFactory factory;

	factory.model = make_shared<T>();

	// get the names of required distributions
	factory.model_names = factory.model->get_names();

	// convert the ordered vector to just a set
	for(const std::string &iter : factory.model_names)
		factory.remaining_names.emplace(to_lower(iter));

	// set the size of the model's vector of distributions
	factory.model->dists.resize(factory.model->get_num_parameters());

	return factory;
}

} // namespace molstat

#endif