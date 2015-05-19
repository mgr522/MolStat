/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file echem_fit_module.cc
 * \brief Function that loads electrochemistry fit models.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include "echem_fit_module.h"

using namespace std;

namespace molstat {
namespace echem {

void load_models(
	std::map<std::string, FitModelFactory<1>> &models)
{
#if 0
// preserve the code for ease of adding a model, once one's ready
	models["symmetricresonant"] =
		GetFitModelFactory<SymmetricResonantFitModel, 1>();
#endif
}

} // namespace molstat::echem
} // namespace molstat
