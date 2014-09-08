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
 * \brief Loads the transport models into the MolStat ``database''.
 *
 * Models are stored as a map from a string (the name of the model) to a
 * function that creates an instance of the model. The function has signature
 * function< shared_ptr<SimulatorModel>
 *   ( map<string, shared_ptr<RandomDistribution>> ) >.
 * The argument is a map of available distributions, as required by the
 * constructors.
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for transport have been added to it.
 */
void load_transport_models(std::map<std::string,
	std::function<shared_ptr<SimulateModel>(
		const std::map<std::string, shared_ptr<RandomDistribution>>&)>> &models);

#endif
