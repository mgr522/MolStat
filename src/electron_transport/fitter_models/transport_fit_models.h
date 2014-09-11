/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_fit_models.h
 * \brief Specializations of the FitModel template for the 1-dimensional fit
 *    functions used to fit electron transport data.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __transport_fit_models_h__
#define __transport_fit_models_h__

#include <general/fitter_tools/fit_model_interface.h>
#include <map>
#include <string>

/**
 * \brief Loads the transport models into the MolStat "database".
 *
 * Models are stored as a map from a string (the name of the model) to a
 * function that creates an instance of the model
 * (of type FitModelInstantiator).
 *
 * \param[in,out] models The map of models in MolStat. On output, the models
 *    for transport have been added to it.
 */
void load_transport_models(
	std::map<std::string, FitModelInstantiator<1>> &models);

#endif
