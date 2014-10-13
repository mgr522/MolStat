/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_simulate_models.cc
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "transport_simulate_models.h"
#include "transport_observables.h"
#include "sym_one_site_simulate_model.h"
#include "asym_one_site_simulate_model.h"
 #if 0
#include "ind_two_chan_simulate_model.h"
#endif
#include "sym_two_site_simulate_model.h"
#include "asym_two_site_simulate_model.h"

namespace molstat {
namespace transport {

void load_models(
	std::map<std::string,
	         SimulatorFactory<2>> &models) {

	models["symmetriconesite"] =
		GetSimulatorFactory<2, SymOneSiteSimulateModel>();

	models["asymmetriconesite"] =
		GetSimulatorFactory<2, AsymOneSiteSimulateModel>();

#if 0
	models["independenttwochannel"] =
		GetSimulateModelFunction<2, IndTwoChanSimulateModel>();
#endif

	models["symmetrictwosite"] =
		GetSimulatorFactory<2, SymTwoSiteSimulateModel>();

	models["asymmetrictwosite"] =
		GetSimulatorFactory<2, AsymTwoSiteSimulateModel>();
}

void load_observables(
	std::map<std::string,
	         ObservableSetter<2>> &observables) {

	observables["appliedbias"] =
		GetObservableSetter<2, AppliedBias>();

	observables["staticconductance"] =
		GetObservableSetter<2, StaticConductance>();

	observables["differentialconductance"] =
		GetObservableSetter<2, DifferentialConductance>();
}

} // namespace molstat::transport
} // namespace molstat