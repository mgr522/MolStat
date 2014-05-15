/**
 * \file fit_model_1d.h
 * \brief Specializations of the FitModel template for 1-dimensional fit
 *    functions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __fit_model_1d_h__
#define __fit_model_1d_h__

#include "fit_model.h"

/**
 * \brief Gets a FitModel<1> object from a tokenized string (name of the model)
 *    and a data set. These types of models can fit 1D conductance histograms.
 *
 * \exception std::invalid_argument if the model name is not found.
 *
 * \param[in] name The name of the FitModel.
 * \param[in] data The data points to fit against.
 * \return Pointer to the FitModel.
 */
std::shared_ptr<FitModel<1>> get_fit_model(const std::string &name,
	std::list<std::pair<std::array<double, 1>, double>> &data);

#endif
