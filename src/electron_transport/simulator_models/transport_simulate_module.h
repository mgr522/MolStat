/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_simulate_module.h
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __transport_simulate_module_h__
#define __transport_simulate_module_h__

#include <string>
#include <map>
#include <general/simulator_tools/simulate_model.h>

namespace molstat {
namespace transport {

/**
 * \brief Loads the transport models into the MolStat "database".
 *
 * Models are stored as a map from a string (the name of the model) to a
 * function that creates the desired molstat::SimulateModelFactory.
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for transport have been added to it.
 */
void load_models(
	std::map<std::string,
	         SimulateModelFactoryFunction> &models);

/**
 * \brief Loads the transport observables into the MolStat "database".
 *
 * Observables are stored as a map from a string (the name of the observable)
 * to an index that can be used to add the observable to a molstat:Simulator.
 *
 * \param[in,out] observables The map of observables in MolStat. On output, the
 *    observables for transport have been added to it.
 */
void load_observables(
	std::map<std::string,
	         ObservableIndex> &observables);

} // namespace molstat::transport
} // namespace molstat

#endif
