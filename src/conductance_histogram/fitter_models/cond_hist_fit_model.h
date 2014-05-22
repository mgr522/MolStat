/**
 * \file cond_hist_fit_model.h
 * \brief Specializations of the FitModel template for the 1-dimensional fit
 *    functions used to fit conductance histograms.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __cond_hist_fit_model_h__
#define __cond_hist_fit_model_h__

#include <general/fitter_tools/fit_model_interface.h>

/**
 * \brief Gets a FitModel<1> object from a tokenized string (name of the model)
 *    and a data set for fitting a 1D conductance histograms.
 *
 * \exception std::invalid_argument if the model name is not found.
 *
 * \param[in] name The name of the FitModel.
 * \param[in] data The data points to fit against.
 * \return Pointer to the FitModel.
 */
std::shared_ptr<FitModel<1>> get_cond_hist_fit_model(const std::string &name,
	std::list<std::pair<std::array<double, 1>, double>> &data);

#endif
