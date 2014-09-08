/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file cond_hist_fit_model.cc
 * \brief Implementation of the FitModel interface for its 1-dimensional
 *    (i.e., the fit function maps \f$\mathbb{R}\f$ to \f$\mathbb{R}\f$)
 *    specialization.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 * \endinternal
 */

#include "cond_hist_fit_model.h"
#include "symmetric_resonant.h"
#include "symmetric_nonresonant.h"
#include "asymmetric_resonant.h"
#include <exception>

using namespace std;

std::shared_ptr<FitModel<1>> get_cond_hist_fit_model(const std::string &name,
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
