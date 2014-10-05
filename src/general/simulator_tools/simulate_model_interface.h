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
 * The subsequent class molstat::SimulateObservables does most of the heavy
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
 * In this sense, each molstat::Observable is essentially an interface that
 * is meant to be inherited, along with molstat::SimulateModel.
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
 * \brief Base class for a model that using MPs model parameters for
 *    calculating observables.
 *
 * This class stores all of the random number distributions for simulating
 * data and provides a function for generating a random set of model
 * parameters.
 *
 * All models for simulating data should derive from this class. Derived classes
 * should pass in the names of parameters they require (during construction)
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

// declaration
template<std::size_t OBS>
class SimulatorFactory;

/**
 * \brief Base class for a model that has MPs parameters and calculates OBS
 *    observables.
 *
 * This class invokes the appropriate observable functions for the
 * molstat::SimulatorInterface.
 *
 * \tparam OBS The number of observables.
 * \tparam MPs The number of model parameters needed to calculate an
 *    observable.
 */
template<std::size_t OBS, std::size_t MPs>
class SimulateObservables : public Simulator<OBS> {
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
	 * \brief The Observable functions.
	 */
	std::array<ObservableFunction, OBS> observables;

	/**
	 * \brief Default constructor.
	 *
	 * Hidden so that objects must be created using the factory function.
	 */
	SimulateObservables() = default;

public:
	virtual ~SimulateObservables() = default;

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
};

/**
 * \brief Class used to aid in the construction of molstat::Simulator objects
 *    (as implemented by molstat::SimulateObservables).
 *
 * \tparam OBS The number of observables.
 */
template<std::size_t OBS>
class SimulatorFactory {
private:
	/**
	 * \brief the molstat::SimulateObservables object that we're building.
	 */
	std::unique_ptr<Simulator<OBS>> model;

	/**
	 * \brief The maximum number of model parameters a molstat::SimulateModel
	 *    can use.
	 *
	 * I haven't figured out a way to regain the template parameter MPs from
	 * here. (It's only needed when constructing SimulatorObservables...).
	 * Thus, the implementation of this Factory does a hard-coded dynamic_cast
	 * through the options up to (and including) the following value. If a
	 * higher value is needed, it must be increased here (for error checking)
	 * and also in the molstat::SimulatorFactory::setObservable function.
	 */
	constexpr static std::size_t MAX_MPs = 6;

	/**
	 * \brief Casts the model for the specified number of model parameters.
	 *
	 * This is the other part of the inelegant solution for needing access to
	 * MPs when adding observables.
	 *
	 * \throw runtime_error If, once cast to a molstat::SimulateObservables,
	 *    the model is incompatible with the specified observable.
	 *
	 * \tparam MPs The number of model parameters.
	 * \tparam T The class (interface) for the observable.
	 * \param[in,out] ptr Raw pointer to the molstat::Simulator we're creating.
	 *    This pointer will be dynamically cast to a
	 *    molstat::SimulateObservables; if successful, the observable will be
	 *    added in the correct place. ptr will not be deleted.
	 * \param[in] j The index for this observable (j = 0, ..., OBS-1).
	 * \return True if the assignment was successful; false otherwise.
	 */
	template<std::size_t MPs, template<std::size_t> class T>
	static bool setObservableMPs(Simulator<OBS> *ptr, std::size_t j);

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
	 * map of available molstat::RandomDistributions; these distributions.
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
