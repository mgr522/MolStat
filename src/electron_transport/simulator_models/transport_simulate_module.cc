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

#include <general/string_tools.h>
#include "transport_simulate_module.h"
#include "observables.h"
#include "junction.h"
#if 0
#include "sym_one_site_simulate_model.h"
#include "asym_one_site_simulate_model.h"
#include "ind_two_chan_simulate_model.h"
#include "sym_two_site_simulate_model.h"
#include "asym_two_site_simulate_model.h"
#endif

namespace molstat {
namespace transport {

void load_models(
	std::map<std::string,
	         SimulateModelFactory> &models) {

#if 0
	models.emplace(make_pair(to_lower("SymmetricOneSite"),
		GetSimulateModelFactory<SymOneSiteSimulateModel>() ));

	models["symmetriconesite"] =
		GetSimulatorFactory<2, SymOneSiteSimulateModel>();

	models["asymmetriconesite"] =
		GetSimulatorFactory<2, AsymOneSiteSimulateModel>();

	models["independenttwochannel"] =
		GetSimulateModelFunction<2, IndTwoChanSimulateModel>();

	models["symmetrictwosite"] =
		GetSimulatorFactory<2, SymTwoSiteSimulateModel>();

	models["asymmetrictwosite"] =
		GetSimulatorFactory<2, AsymTwoSiteSimulateModel>();
#endif
}

void load_observables(
	std::map<std::string,
	         ObservableIndex> &observables) {

	observables.emplace(make_pair(to_lower("AppliedBias"),
		GetObservableIndex<AppliedBias>()) );

	observables.emplace(make_pair(to_lower("StaticConductance"),
		GetObservableIndex<StaticConductance>()) );

	observables.emplace(make_pair(to_lower("DifferentialConductance"),
		GetObservableIndex<DifferentialConductance>()) );
}

} // namespace molstat::transport
} // namespace molstat
