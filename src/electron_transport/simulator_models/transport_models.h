/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_models.h
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __transport_models_h__
#define __transport_models_h__

#include <memory>
#include <string>
#include <map>
#include <functional>
#include <general/random_distributions/rng.h>
#include <general/simulator_tools/simulate_model_interface.h>

using std::shared_ptr;

/**
 * \brief Loads the transport models into the MolStat "database".
 *
 * Models are stored as a map from a string (the name of the model) to a
 * function that creates an instance of the model (of type
 * SimulateModelInstantiator).
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for transport have been added to it.
 */
void load_transport_models(
	std::map<std::string, SimulateModelInstantiator> &models);

/**
 * \brief Loads the transport observables into the MolStat "database".
 *
 * Observables are stored as a map from a string (the name of the observable)
 * to a function that checks the model/observable pair and returns an
 * Observable<2>. 1D observables have a 0. padded onto them for programming
 * convenience. The observable function takes the GSL random number generator,
 * simulates parameters, and calculates the observable using the specified
 * model.
 *
 * The extra function from shared_ptr<SimulateModel> -> Observable is needed
 * to typecheck the model/observable combination.
 *
 * \param[in,out] observables The map of observables in MolStat. On output, the
 *    observables for transport have been added to it.
 */
void load_transport_observables(std::map<std::string, std::function<
	Observable<2>(const shared_ptr<SimulateModel>)>> &observables);

#endif
