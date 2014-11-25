/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file echem_simulate_module.h
 * \brief Functions that load the electrochemistry models and observables
 *    into MolStat.
 *
 * \author Bo Fu, Matthew G.\ Reuter
 * \date November 2014
 */

#ifndef __echem_simulate_module_h__
#define __echem_simulate_module_h__

#include <string>
#include <map>
#include <general/simulator_tools/simulate_model.h>

namespace molstat {
namespace echem {

/**
 * \brief Loads the electrochemistry models into the MolStat "database".
 *
 * Models are stored as a map from a string (the name of the model) to a
 * function that creates the desired molstat::SimulateModelFactory.
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for electrochemistry have been added to it.
 */
void load_models(
	std::map<std::string,
	         SimulateModelFactoryFunction> &models);

/**
 * \brief Loads the electrochemistry observables into the MolStat "database".
 *
 * Observables are stored as a map from a string (the name of the observable)
 * to an index that can be used to add the observable to a molstat:Simulator.
 *
 * \param[in,out] observables The map of observables in MolStat. On output, the
 *    observables for electrochemistry have been added to it.
 */
void load_observables(
	std::map<std::string,
	         ObservableIndex> &observables);

} // namespace molstat::echem
} // namespace molstat

#endif
