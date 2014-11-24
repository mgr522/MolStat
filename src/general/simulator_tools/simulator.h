/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulator.h
 * \brief Defines the molstat::Simulator class for simulating histograms.
 *
 * A fairly high-level introduction to the classes for simulating histograms
 * can be found in \ref sec_add_simulate. Here we expand on this a little bit.
 *
 * The molstat::Simulator class is the interface for the whole package. A
 * molstat::Simulator contains the model (and the possible hierarchy of
 * submodels therein) and the list of functions for calculating the desired
 * observables. The molstat::SimulateModel (top-level model) must be specified
 * when constructing the molstat::Simulator. The molstat::Simulator::simulate
 * function then simulates a random set of model parameters and calculates the
 * requested observables.
 *
 * Details for the molstat::SimulateModel class and molstat::Observable
 * template are found in the documentation of simulate_model.h.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __simulator_h__
#define __simulator_h__

#include <memory>
#include <valarray>
#include <vector>
#include <typeindex>
#include <general/random_distributions/rng.h>
#include "simulate_model.h"

namespace molstat {

/**
 * \brief Class for simulating data.
 *
 * This class contains the actual model used to simulate data and functions
 * for constructing the model, setting observables, etc.
 */
class Simulator
{
private:
	/// The actual model used to simulate data.
	std::shared_ptr<SimulateModel> model;

	/**
	 * \brief The functions that calculate observables, given a set of model
	 *    parameters.
	 */
	std::vector<ObservableFunction> obs_functions;

public:
	Simulator() = delete;

	/**
	 * \brief Constructs molstat::Simulator using the specified model.
	 *
	 * \throw runtime_error if a nullptr is provided.
	 * \throw molstat::FullModelRequired if the specified model is a
	 *    submodel type.
	 *
	 * \param[in] model_ The model to be used.
	 */
	Simulator(std::shared_ptr<SimulateModel> model_);

	/**
	 * \brief Calculates the desired observables using the random number
	 *    generator.
	 *
	 * \throw molstat::NoObservables if no observables have been set.
	 * \throw molstat::MissingDistribution if any of the model's required
	 *    distributions is unspecified.
	 *
	 * \param[in] engine The C++11 random number engine.
	 * \return The simulated observables.
	 */
	std::valarray<double> simulate(Engine &engine) const;

	/**
	 * \brief Sets the `j`th observable for the simulator.
	 *
	 * First verify that the observable is compatible with the specified model.
	 * If compatible, get a function for calculating the observable and store
	 * it for use later.
	 *
	 * \throw molstat::IncompatibleObservable If our models is incompatible
	 *    with the desired observable.
	 * \throw out_of_range If `j` is out of range (not between 0 and
	 *    obs_functions.size()).
	 *
	 * \param[in] j The output index for this observable.
	 * \param[in] obs The identifier of the observable.
	 */
	void setObservable(std::size_t j, const ObservableIndex &obs);
};

} // namespace MolStat

#endif