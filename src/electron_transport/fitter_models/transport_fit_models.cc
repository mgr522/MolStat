/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file transport_fit_models.cc
 * \brief Function that loads transport fit models.
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

namespace molstat {

void load_transport_models(
	std::map<std::string, FitModelFactory<1>> &models)
{
	models["symmetricresonant"] =
		GetFitModelFactory<SymmetricResonantFitModel, 1>();

	models["symmetricnonresonant"] =
		GetFitModelFactory<SymmetricNonresonantFitModel, 1>();

	models["asymmetricresonant"] =
		GetFitModelFactory<AsymmetricResonantFitModel, 1>();
}

} // namespace molstat