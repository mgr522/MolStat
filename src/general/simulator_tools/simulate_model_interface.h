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
public:
	/**
	 * \brief Shortcut for a function that calculates an observable.
	 */
	using ObservableFunction =
		std::function<double(const std::array<double, MPs> &)>;

private:
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
	 * \brief Factory function for creating a molstat::SimulateObservables
	 *    object.
	 *
	 * SimulateObservables objects are created by passing in a map of available
	 * molstat::RandomDistributions; these distributions are needed for
	 * constructing the underlying molstat::SimulateModel.
	 *
	 * \throw runtime_error If there a required distribution is not found among
	 *    the available distributions.
	 *
	 * \tparam T The type of model (derived from molstat::SimulateModel) to use.
	 * \param[in] avail The available distributions, keyed by name,
	 */
	template<typename T>
	static std::unique_ptr<SimulateObservables<OBS, MPs>> Factory(
		const std::map<std::string,
		               std::shared_ptr<RandomDistribution>> &avail) {

		using namespace std;

		unique_ptr<SimulateObservables<OBS, MPs>> ret(
			new SimulateObservables<OBS, MPs>());

		// make the model from the specified type
		ret->model = make_shared<T>(avail);

		// initialize all observables to the zero function
		for(std::size_t j = 0; j < OBS; ++j)
			ret->observables[j] = ZeroObs;

		return ret;
	}

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
	 * \brief Sets the `i`th observable for the simulation.
	 *
	 * This function will dynamic_cast the model to make sure the specified
	 * function is available; that is, the model implements the desired
	 * observable. If this check is passed, the observable is set. Otherwise,
	 * a runtime_error exception is thrown.
	 *
	 * \throw runtime_error If the model and observable are incompatible; that is,
	 *    the model does not implement the observable.
	 *
	 * \tparam T The class (interface) for the observable.
	 * \param[in] i The index for this observable (i = 0... OBS-1).
	 */
	template<template<std::size_t> class T>
	void setObservable(std::size_t i) {

		std::shared_ptr<T<MPs>> cast = std::dynamic_pointer_cast<T<MPs>>(model);

		if(cast == nullptr)
			throw std::runtime_error("Incompatible model and observable.");

		observables[i] =
			[cast] (const std::array<double, MPs> &params)
				-> double {

				return (cast->operator())(params);
			};
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
};

} // namespace molstat

// Owing to templates, include the *.cc implementation file.
#include "simulate_model_interface.cc"

#endif
