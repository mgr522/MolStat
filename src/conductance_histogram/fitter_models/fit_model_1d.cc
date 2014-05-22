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
#include "symmetric_resonant.h"
#include "symmetric_nonresonant.h"
#include "asymmetric_resonant.h"
#include <exception>

using namespace std;

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
	std::list<std::pair<std::array<double, 1>, double>> &data) {

	shared_ptr<FitModel<1>> ret;

	if(name == "symmetricresonant")
		ret = make_shared<SymmetricResonantFitModel>(data);
	else if(name == "symmetricnonresonant")
		ret = make_shared<SymmetricNonresonantFitModel>(data);
	else if(name == "asymmetricresonant")
		ret = make_shared<AsymmetricResonantFitModel>(data);
	else {
		throw invalid_argument("Error: unrecognized fit model: '%s'.\n" \
			"Recognized options are\n" \
			"   - 'SymmetricResonant'\n" \
			"     resonant tunneling through a single site, symmetric coupling.\n" \
			"   - 'SymmetricNonresonant'\n" \
			"     nonresonant tunneling through a single site, symmetric coupling.\n" \
			"   - 'AsymmetricResonant'\n" \
			"     resonant tunneling through a single site, asymmetric coupling."
		);
	}

	return ret;
}
