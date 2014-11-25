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

#include <general/string_tools.h>
//#include "transport_simulate_module.h"
#include "observables.h"

//#include "asym_two_site_channel.h"

namespace molstat {
namespace echem {

void load_models(
	std::map<std::string,
	         SimulateModelFactoryFunction> &models)
{
#if 0
	models.emplace(
		to_lower("TransportJunction"),
		GetSimulateModelFactory<TransportJunction>() );
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
