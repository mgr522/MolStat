/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_fit_module.h
 * \brief Loads transport models for fitting electron transport data.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __transport_fit_models_h__
#define __transport_fit_models_h__

#include <general/fitter_tools/fit_model_interface.h>
#include <map>
#include <string>

namespace molstat {
namespace transport {

/**
 * \brief Loads the transport models into the MolStat "database".
 *
 * Models are stored as a map from a string (the name of the model) to a
 * factory for the model (of type FitModelFactory).
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for transport have been added to it.
 */
void load_models(
	std::map<std::string, FitModelFactory<1>> &models);

} // namespace molstat::transport
} // namespace molstat

#endif
