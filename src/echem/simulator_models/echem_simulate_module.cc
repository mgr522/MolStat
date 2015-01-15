/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file echem_simulate_module.cc
 * \brief Functions that load the electrochemistry models and observables into
 *    MolStat.
 *
 * \author Bo Fu, Matthew G.\ Reuter
 * \date November 2014
 */

#include <config.h>
#include <general/string_tools.h>
#include "observables.h"
 
#include "non-nernstian.h"

namespace molstat {
namespace echem {

void load_models(
	std::map<std::string,
	         SimulateModelFactoryFunction> &models)
{
#if HAVE_CVODE
	// models that require the optional CVODE package
	models.emplace(
		to_lower("NonNernstianReaction"),
		GetSimulateModelFactory<NonNernstianReaction>() );
#endif
}

void load_observables(
	std::map<std::string,
	         ObservableIndex> &observables)
{
	observables.emplace(
		to_lower("ForwardETPotential"),
		GetObservableIndex<ForwardETPotential>() );

	observables.emplace(
		to_lower("BackwardETPotential"),
		GetObservableIndex<BackwardETPotential>() );
}

} // namespace molstat::echem
} // namespace molstat
