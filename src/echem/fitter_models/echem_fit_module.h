/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file echem_fit_module.h
 * \brief Loads electrochemistry models for fitting experimental data.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#ifndef __echem_fit_models_h__
#define __echem_fit_models_h__

#include <general/fitter_tools/fit_model_interface.h>
#include <map>
#include <string>

namespace molstat {
namespace echem {

/**
 * \brief Loads the electrochemistry models into the MolStat "database".
 *
 * Models are stored as a map from a string (the name of the model) to a
 * factory for the model (of type FitModelFactory).
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for transport have been added to it.
 */
void load_models(
	std::map<std::string, FitModelFactory<1>> &models);

} // namespace molstat::echem
} // namespace molstat

#endif
