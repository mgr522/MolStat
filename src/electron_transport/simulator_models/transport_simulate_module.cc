/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_simulate_module.cc
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <utility>
#include <general/string_tools.h>
#include "transport_simulate_module.h"
#include "observables.h"
#include "junction.h"

#include "sym_one_site_channel.h"
#include "asym_one_site_channel.h"
#include "sym_two_site_channel.h"
#include "asym_two_site_channel.h"

namespace molstat {
namespace transport {

/**
 * \brief Helper function that adds items the map.
 *
 * This function uses emplace to sidestep use of the default constructor.
 *
 * \param[in,out] collection The map. On output, has the item added to it.
 * \param[in] name The name of the item (key in the map).
 * \param[in] item The item to insert.
 */
template<typename T>
inline static void add_entry(std::map<std::string, T> &collection,
                             std::string name,
                             T &&item) {

	// need to use emplace in case T doesn't have a default constructor
	collection.emplace(
		make_pair(to_lower(name), std::forward<T>(item))
	);
}

void load_models(
	std::map<std::string,
	         SimulateModelFactoryFunction> &models) {

	add_entry(models,
		to_lower("TransportJunction"),
		GetSimulateModelFactory<TransportJunction>() );

	add_entry(models,
		to_lower("SymmetricOneSiteChannel"),
		GetSimulateModelFactory<SymOneSiteChannel>() );

	add_entry(models,
		to_lower("AsymmetricOneSiteChannel"),
		GetSimulateModelFactory<AsymOneSiteChannel>() );

	add_entry(models,
		to_lower("SymmetricTwoSiteChannel"),
		GetSimulateModelFactory<SymTwoSiteChannel>() );

	add_entry(models,
		to_lower("AsymmetricTwoSiteChannel"),
		GetSimulateModelFactory<AsymTwoSiteChannel>() );
}

void load_observables(
	std::map<std::string,
	         ObservableIndex> &observables) {

	add_entry(observables,
		to_lower("AppliedBias"),
		GetObservableIndex<AppliedBias>() );

	add_entry(observables,
		to_lower("StaticConductance"),
		GetObservableIndex<StaticConductance>() );

	add_entry(observables,
		to_lower("DifferentialConductance"),
		GetObservableIndex<DifferentialConductance>() );
}

} // namespace molstat::transport
} // namespace molstat
