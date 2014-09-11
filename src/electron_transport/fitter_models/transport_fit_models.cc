/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file transport_fit_models.cc
 * \brief Implementation of the FitModel interface for its 1-dimensional
 *    (i.e., the fit function maps \f$\mathbb{R}\f$ to \f$\mathbb{R}\f$)
 *    specialization.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 * \endinternal
 */

#include "transport_fit_models.h"
#include "symmetric_resonant.h"
#include "symmetric_nonresonant.h"
#include "asymmetric_resonant.h"
#include <exception>

using namespace std;

void load_transport_models(
	std::map<std::string, FitModelInstantiator<1>> &models) {

	models["symmetricresonant"] =
		FitModelAdd<SymmetricResonantFitModel, 1>();

	models["symmetricnonresonant"] =
		FitModelAdd<SymmetricNonresonantFitModel, 1>();

	models["asymmetricresonant"] =
		FitModelAdd<AsymmetricResonantFitModel, 1>();
}
