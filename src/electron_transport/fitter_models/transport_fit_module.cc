/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_fit_module.cc
 * \brief Function that loads transport fit models.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "transport_fit_module.h"
#include "symmetric_resonant.h"
#include "symmetric_nonresonant.h"
#include "asymmetric_resonant.h"

using namespace std;

namespace molstat {
namespace transport {

void load_models(
	std::map<std::string, FitModelFactory<1>> &models)
{
	models["symmetricresonant"] =
		GetFitModelFactory<SymmetricResonantFitModel, 1>();

	models["symmetricnonresonant"] =
		GetFitModelFactory<SymmetricNonresonantFitModel, 1>();

	models["asymmetricresonant"] =
		GetFitModelFactory<AsymmetricResonantFitModel, 1>();
}

} // namespace molstat::transport
} // namespace molstat
