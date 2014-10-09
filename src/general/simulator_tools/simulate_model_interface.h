/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface.h
 * \brief Defines abstract classes for simulating histograms.
 *
 * \todo Write a better description of how all these models work together.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __simulate_model_interface_h__
#define __simulate_model_interface_h__

#include <memory>
#include <array>
#include <string>
#include <map>
#include <stdexcept>
#include <utility>
#include <functional>
#include <general/random_distributions/rng.h>

namespace molstat {

/**
 * \brief Base class for simulating data.
 *
 * This model requires the number of observables being calculated with each
 * set of random numbers.
 *
 * Subclasses must implement the molstat::Simulator::simulate function
 * to generate a set of random numbers (given a generator) and calculate the
 * requested observables.
 *
 * The subsequent class molstat::ModelSimulator does most of the heavy
 * lifting, but is also templated over the number of model parameters. This
 * class provides an easy way to type a generic model.
 *
 * \tparam OBS The number of observables.
 */
template<std::size_t OBS>
class Simulator {
public:
	virtual ~Simulator() = default;

	/**
	 * \brief Calculates the desired observables using the random number
	 *    generator.
	 *
	 * \param[in] r The GSL random number generator handle.
	 * \return The simulated observables.
	 */
	virtual std::array<double, OBS> simulate(gsl_rng_ptr &r) const = 0;
};

/**
 * \brief Base class for an observable.
 *
 * An observable must implement the operator() function (see below). As
 * demonstrated in examples, the operator() function should immediately
 * redirect to another function whose name can be specified by the
 * observable. In this way, a particular model can inherit from multiple
 * observables and implement each of their \"observable\" functions.
 *
 * Each molstat::Observable is essentially an interface that is meant to be
 * inherited, along with molstat::SimulateModel.
 *
 * \tparam MPs The number of model parameters needed to calculate this
 *    observable.
 */
template<std::size_t MPs>
class Observable {
public:
	virtual ~Observable() = default;

	/**
	 * \brief Function that calculates the observable, given a set of
	 *    model parameters.
	 *
	 * \param[in] params The model parameters.
	 * \return The value of the observable for these parameters.
	 */
	virtual double operator()(const std::array<double, MPs> &params) const = 0;
};

/**
 * \brief Base class for a model that uses MPs model parameters to calculate
 *    observables.
 *
 * This class stores all of the random number distributions for simulating
 * data and provides a function for generating a random set of model
 * parameters.
 *
 * All models for simulating data should derive from this class. Derived classes
 * should pass in the names of parameters they require during construction
 * and subsequently process given parameters to calculate an observable.
 *
 * \tparam MPs The number of model parameters needed to calculate an
 *    observable.
 */
template<std::size_t MPs>
class SimulateModel {
private:
	/**
	 * \brief The random distributions, ordered as specified during
	 *    construction.
	 */
	std::array<std::shared_ptr<const RandomDistribution>, MPs> dists;

protected:
	/**
	 * \brief Converts a map laying out the desired order of parameters to
	 *    an array with the parameters so ordered.
	 *
	 * \throw out_of_range If a name for one of the indices O, 1, ..., MPs-1
	 *    is missing.
	 *
	 * \param[in] param_order Map that orders the parameter names.
	 * \return An array with the names in order.
	 */
	std::array<std::string, MPs> order_from_map(
		const std::map<std::size_t, std::string> &param_order);

public:
	SimulateModel() = delete;
	virtual ~SimulateModel() = default;

	/**
	 * \brief The number of model parameters, for reference at compile/runtime.
	 */
	constexpr static std::size_t numModelParameters = MPs;

	/**
	 * \brief Constructor requiring a list of available distributions and an
	 *    array of required distributions.
	 *
	 * The order of parameters in the list `names` will be preserved when
	 * randomly generating a set of parameters.
	 *
	 * \throw runtime_error If a required distribution is not found among the
	 *    available distributions.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 * \param[in] names The names of required distributions.
	 */
	SimulateModel(
		const std::map<std::string,
		               std::shared_ptr<RandomDistribution>> &avail,
		const std::array<std::string, MPs> &names);

	/**
	 * \brief Produces a random set of parameters from the specified
	 *    distributions.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \return The randomly-generated model parameters, in the order specified
	 *    during construction.
	 */
	std::array<double, MPs> generateParameters(gsl_rng_ptr &r) const;
};

// forward declarations of classes used to construct molstat::ModelSimulator
// objects. They need to be friends of molstat::ModelSimulator
template<std::size_t OBS>
class SimulatorFactory;

namespace SimulatorFactoryHelper {
template<std::size_t MPs>
class ObservableSetter;
} // namespace molstat::SimulatorFactoryHelper

/**
 * \brief Base class for a simulator that employs a model with MPs parameters
 *    to calculate OBS observables.
 *
 * This class stores and, when requested, invokes the appropriate observable
 * functions for simulating observables.
 *
 * \tparam OBS The number of observables.
 * \tparam MPs The number of model parameters needed to calculate an
 *    observable.
 */
template<std::size_t OBS, std::size_t MPs>
class ModelSimulator : public Simulator<OBS> {
private:
	/**
	 * \brief Shortcut for a function that calculates an observable.
	 */
	using ObservableFunction =
		std::function<double(const std::array<double, MPs> &)>;

	/**
	 * \brief The model used to calculate observables.
	 */
	std::shared_ptr<SimulateModel<MPs>> model;

	/**
	 * \brief The observable functions.
	 */
	std::array<ObservableFunction, OBS> observables;

	/**
	 * \brief Default constructor.
	 *
	 * Hidden so that objects must be created using the factory function.
	 */
	ModelSimulator() = default;

public:
	virtual ~ModelSimulator() = default;

	/**
	 * \brief Function for a \"zero\" observable.
	 *
	 * \param[in] params The set of model parameters.
	 * \return The observable; in this case, 0.
	 */
	constexpr static double ZeroObs(const std::array<double, MPs> &params)
		noexcept {

		return 0.;
	}

	/**
	 * \brief Calculates the desired observables using the random number
	 *    generator.
	 *
	 * \param[in] r The GSL random number generator handle.
	 * \return The simulated observables.
	 */
	virtual std::array<double, OBS>
		simulate(gsl_rng_ptr &r) const override final;

	// Give the SimulatorFactory class access to the internals so that it can
	// produce a Simulator.
	friend class SimulatorFactory<OBS>;
	friend class SimulatorFactoryHelper::ObservableSetter<MPs>;
};

/**
 * \brief Class used to aid in the construction of molstat::Simulator objects
 *    (as implemented by molstat::ModelSimulator).
 *
 * \tparam OBS The number of observables.
 */
template<std::size_t OBS>
class SimulatorFactory {
private:
	/**
	 * \brief the molstat::ModelSimulator object that we're building.
	 */
	std::unique_ptr<Simulator<OBS>> model;

	/**
	 * \brief The maximum number of model parameters a molstat::SimulateModel
	 *    can use.
	 *
	 * I haven't figured out a way to regain the template parameter MPs from
	 * here. (It's only needed when constructing molstat::ModelSimulator
	 * objects). Thus, the implementation of this Factory does a hard-coded
	 * dynamic_cast through the options up to (and including) the following
	 * value. If a higher value is needed, it simply must be increased here.
	 */
	constexpr static std::size_t MAX_MPs = 6;

public:
	SimulatorFactory() : model(nullptr) {}
	SimulatorFactory(SimulatorFactory<OBS> &&) = default;
	~SimulatorFactory() = default;
	SimulatorFactory &operator=(SimulatorFactory<OBS> &&) = default;

	/**
	 * \brief Constructs a molstat::SimulatorFactory that uses a model of type
	 *    T.
	 *
	 * Construction of the SimulateModel underlying the Simulator requires a
	 * map of available molstat::RandomDistributions.
	 *
	 * \throw runtime_error If a required distribution is not found among the
	 *    available distributions.
	 *
	 * \tparam T The type of model (derived from molstat::SimulateModel) to use.
	 * \param[in] avail The available distributions, keyed by name.
	 */
	template<typename T>
	static SimulatorFactory<OBS> makeFactory(
		const std::map<std::string,
		               std::shared_ptr<RandomDistribution>> &avail);

	/**
	 * \brief Sets the `i`th observable for the simulator.
	 *
	 * This function will dynamic_cast the model to make sure the specified
	 * function is available; that is, the model implements the desired
	 * observable. If so, the observable is set. Otherwise, a runtime_error
	 * exception is thrown.
	 *
	 * \throw out_of_range If j is not in 0, ..., OBS-1.
	 * \throw runtime_error If the model and observable are incompatible; that is,
	 *    the model does not implement the observable.
	 *
	 * \tparam T The class (interface) for the observable.
	 * \param[in] j The index for this observable (j = 0, ..., OBS-1).
	 */
	template<template<std::size_t> class T>
	SimulatorFactory<OBS> &setObservable(std::size_t j);

	/**
	 * \brief Transfers the simulator being created.
	 *
	 * \return The simulator.
	 */
	std::unique_ptr<Simulator<OBS>> create() noexcept;
};

namespace SimulatorFactoryHelper {
/**
 * \brief This auxiliary class determines (recursively) the correct number of
 *    model parameters (MPs) for the model being constructed and, once
 *    determined, adds the desired observable to the molstat::Simulator
 *    object (via molstat::ModelSimulator).
 *
 * This is the other part of the inelegant solution for needing access to
 * MPs when adding observables.
 *
 * This class is recursive in the template parameter MPs. The invoking
 * function should call the member function with an upper bound for MPs. If
 * this value works, the observable is assigned and we return. If not, we
 * call the function with one less MPs and repeat. A base case (MPs=0) is
 * provided as a partial specialization.
 *
 * \throw runtime_error If, once cast to a molstat::ModelSimulator, the
 *    model is incompatible with the specified observable.
 *
 * \tparam MPs The number of model parameters to try.
 */
template<std::size_t MPs>
class ObservableSetter {
public:
	/**
	 * \brief Function that attempts to cast the molstat::Simulator object
	 *    to a molstat::ModelSimulator object with MPs model parameters, in
	 *    order to add an observable.
	 *
	 * \tparam OBS The number of observables in the molstat::Simulator.
	 * \tparam T The class (interface) for the observable.
	 * \param[in,out] ptr Raw pointer to the molstat::Simulator we're creating.
	 *    This pointer will be dynamically cast to a molstat::ModelSimulator;
	 *    if successful, the observable will be added in the correct place.
	 *    ptr will not be deleted.
	 * \param[in] j The index for this observable (j = 0, ..., OBS-1).
	 * \return True if the assignment was successful (for this or any other
	 *    value of MPs); false otherwise.
	 */
	template<std::size_t OBS, template<std::size_t> class T>
	inline static bool setObservableMPs(Simulator<OBS> *ptr, std::size_t j);
};

/**
 * \brief Partial specialization of
 *    molstat::SimulatorFactoryHelper::ObservableSetter for the case of MPS=0.
 */
template<>
class ObservableSetter<0> {
public:
	/**
	 * \brief There's nothing to do at MPs=0. We shouldn't ever be here. But, if
	 *    we are, something went wrong and the cast failed.
	 *
	 * \tparam OBS The number of observables in the molstat::Simulator.
	 * \tparam T The class (interface) for the observable.
	 * \param[in,out] ptr Raw pointer to the molstat::Simulator we're creating.
	 *    This pointer will be dynamically cast to a molstat::ModelSimulator;
	 *    if successful, the observable will be added in the correct place.
	 *    ptr will not be deleted.
	 * \param[in] j The index for this observable (j = 0, ..., OBS-1).
	 * \return True if the assignment was successful (for this or any other
	 *    value of MPs); false otherwise.
	 */
	template<std::size_t OBS, template<std::size_t> class T>
	inline constexpr static bool setObservableMPs(Simulator<OBS> *ptr,
		std::size_t j) noexcept;
};
} // namespace molstat::SimulatorFactoryHelper

// Helper types for use in the simulator
/**
 * \brief Shortcut for a function that creates a molstat::SimulatorFactory.
 *
 * \tparam OBS The number of observables.
 */
template<std::size_t OBS>
using SimulateModelFunction = std::function<
	SimulatorFactory<OBS>(
		const std::map<std::string,
		               std::shared_ptr<RandomDistribution>> &)>;

/**
 * \brief Gets a function that creates a molstat::SimulatorFactory for the
 *    given model and number of observables.
 *
 * \tparam OBS The number of observables.
 * \tparam T The model.
 */
template<std::size_t OBS, typename T>
SimulateModelFunction<OBS> GetSimulateModelFunction();

/**
 * \brief Shortcut for a function that sets an observable in
 *    molstat::SimulatorFactory.
 *
 * \tparam OBS The number of observables.
 */
template<std::size_t OBS>
using ObservableFunction = std::function<
	void(SimulatorFactory<OBS> &, std::size_t)>;

/**
 * \brief Gets a function that sets an observable in the
 *    molstat::Simulator.
 *
 * \tparam OBS The number of observables in the molstat::Simulator.
 * \tparam T The class (interface) of the observable.
 */
template<std::size_t OBS, template<std::size_t> class T>
ObservableFunction<OBS> GetObservableFunction();

} // namespace molstat

// Owing to the use of templates, include the *.cc implementation file.
#include "simulate_model_interface.cc"

#endif
