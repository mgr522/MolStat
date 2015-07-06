/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file transport_simulate_module.cc
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <general/string_tools.h>
#include "transport_simulate_module.h"
#include "observables.h"
#include "junction.h"

#include "sym_one_site_channel.h"
#include "asym_one_site_channel.h"
#include "sym_two_site_channel.h"
#include "asym_two_site_channel.h"
#include "rectangular_barrier.h"
#include "sym_interference.h"

namespace molstat {
namespace transport {

void load_models(
	std::map<std::string,
	         SimulateModelFactoryFunction> &models)
{
	models.emplace(
		to_lower("TransportJunction"),
		GetSimulateModelFactory<TransportJunction>() );

	models.emplace(
		to_lower("SymmetricOneSiteChannel"),
		GetSimulateModelFactory<SymOneSiteChannel>() );

	models.emplace(
		to_lower("AsymmetricOneSiteChannel"),
		GetSimulateModelFactory<AsymOneSiteChannel>() );

	models.emplace(
		to_lower("SymmetricTwoSiteChannel"),
		GetSimulateModelFactory<SymTwoSiteChannel>() );

	models.emplace(
		to_lower("AsymmetricTwoSiteChannel"),
		GetSimulateModelFactory<AsymTwoSiteChannel>() );

	models.emplace(
		to_lower("RectangularBarrierChannel"),
		GetSimulateModelFactory<RectangularBarrier>() );

	models.emplace(
		to_lower("InterferenceChannel"),
		GetSimulateModelFactory<SymInterferenceChannel>() );
}

void load_observables(
	std::map<std::string,
	         ObservableIndex> &observables)
{
	observables.emplace(
		to_lower("AppliedBias"),
		GetObservableIndex<AppliedBias>() );

	observables.emplace(
		to_lower("ElectricCurrent"),
		GetObservableIndex<ElectricCurrent>() );

	observables.emplace(
		to_lower("StaticConductance"),
		GetObservableIndex<StaticConductance>() );

	observables.emplace(
		to_lower("ZeroBiasConductance"),
		GetObservableIndex<ZeroBiasConductance>() );

	observables.emplace(
		to_lower("DifferentialConductance"),
		GetObservableIndex<DifferentialConductance>() );

	observables.emplace(
		to_lower("SeebeckCoefficient"),
		GetObservableIndex<SeebeckCoefficient>() );
}

} // namespace molstat::transport
} // namespace molstat
