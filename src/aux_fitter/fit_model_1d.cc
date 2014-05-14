/**
 * \file fit_model_1d.cc
 * \brief Implementation of the FitModel interface for its 1-dimensional
 *    (i.e., the fit function maps \f$\mathbb{R}\f$ to \f$\mathbb{R}\f$)
 *    specialization.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "fit_model_1d.h"

/**
 * \brief Gets a FitModel<1> object from a tokenized string (name of the model)
 *    and a data set. These types of models can fit 1D conductance histograms.
 *
 * \exception std::invalid_argument if the model name is not found.
 *
 * \param[in] name The name of the FitModel.
 * \param[in] data The data that will be fit against.
 * \return Pointer to the FitModel.
 */
std::shared_ptr<FitModel<1>> get_fit_model(const std::string &name,
	std::vector<std::pair<std::array<double, 1>, double>> &data) {

	return nullptr;
}
