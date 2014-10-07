/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_simulate_models.h
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __transport_simulate_models_h__
#define __transport_simulate_models_h__

#include <memory>
#include <string>
#include <map>
#include <functional>
#include <general/random_distributions/rng.h>
#include <general/simulator_tools/simulate_model_interface.h>

namespace molstat {

/**
 * \brief Loads the transport models into the MolStat "database".
 *
 * Models are stored as a map from a string (the name of the model) to a
 * function that creates a molstat::SimulatorFactory.
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for transport have been added to it.
 */
void load_transport_models(
	std::map<std::string,
	         SimulateModelFunction<2>> &models);

/**
 * \brief Loads the transport observables into the MolStat "database".
 *
 * Observables are stored as a map from a string (the name of the observable)
 * to a function that adds the observable to a molstat:SimulatorFactory.
 *
 * \param[in,out] observables The map of observables in MolStat. On output, the
 *    observables for transport have been added to it.
 */
void load_transport_observables(
	std::map<std::string,
	         ObservableFunction<2>> &observables);

} // namespace molstat

#endif
