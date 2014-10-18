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
#include <typeindex>
#include <general/random_distributions/rng.h>

namespace molstat {

/**
 * \brief The signature of a function that calculates an observable.
 */
using ObservableFunction =
	std::function<double(const std::valarray<double> &)>;

/**
 * \brief Base class for a model that uses model parameters to calculate
 *    observables.
 *
 * All models for simulating data should derive from this class. Derived
 * classes must implement functions that provide the number of model parameters
 * needed, as well as a name for each parameter. Derived classes should also
 * probably derive from molstat::Observable so that the MolStat simulator knows
 * that this model is compatible with that observable.
 */
class SimulateModel : public std::enable_shared_from_this<SimulateModel> {
protected:
	/**
	 * \brief Set of observables compatible with this model, keyed by
	 *    the type_index of the class for the observable.
	 */
	std::map<std::type_index, ObservableFunction> compatible_observables;

	/**
	 * \brief Ordered vector of random number distributions for the various
	 *    model parameters.
	 */
	std::vector<std::shared_ptr<const RandomDistribution>> dists;

	/**
	 * \brief Tells whether or not all required distributions have been
	 *    supplied.
	 */
	bool distributions_specified = false;

	/**
	 * \brief Gets a map of parameter name to index.
	 *
	 * \return The ordered list of names of distributions.
	 */
	virtual std::vector<std::string> get_names() const = 0;

public:
	SimulateModel() = default;
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
	 * compatible and then returns a function that casts *this to the
	 * observable type and calls the proper function.
	 *
	 * \throw molstat::IncompatibleObservable if the desired observable is
	 *    incompatible with this model.
	 *
	 * \param[in] obs The type_index of the class for the observable.
	 * \return A function that calculates the observable.
	 */
	virtual ObservableFunction getObservableFunction(
		const std::type_index &obs) const;

	/**
	 * \brief Generates a set of model parameters using the specified random
	 *    distributions.
	 *
	 * \param[in] r Handle to the GSL random number generator.
	 * \return A set of model parameters.
	 */
	virtual std::valarray<double> generateParameters(gsl_rng_ptr &r) const;

	/**
	 * \brief Sets the distribution for the specified parameter.
	 *
	 * \param[in] name The name of the parameter.
	 * \param[in] dist The distribution.
	 */
	void setDistribution(const std::string &name,
		std::shared_ptr<const RandomDistribution> dist);
};

} // namespace molstat

#endif