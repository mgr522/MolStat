/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_simulate_models.cc
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "transport_simulate_models.h"
#include "transport_observables.h"
#include "sym_one_site_simulate_model.h"
#include "asym_one_site_simulate_model.h"
#include "ind_two_chan_simulate_model.h"
#include "sym_two_site_simulate_model.h"
#include "asym_two_site_simulate_model.h"

namespace molstat {

void load_transport_models(
	std::map<std::string, SimulateModelFactory> &models) {

	models["symmetriconesite"] =
		GetSimulateModelFactory<SymOneSiteSimulateModel>();

	models["asymmetriconesite"] =
		GetSimulateModelFactory<AsymOneSiteSimulateModel>();

	models["independenttwochannel"] =
		GetSimulateModelFactory<IndTwoChanSimulateModel>();

	models["symmetrictwosite"] =
		GetSimulateModelFactory<SymTwoSiteSimulateModel>();

	models["asymmetrictwosite"] =
		GetSimulateModelFactory<AsymTwoSiteSimulateModel>();
}

void load_transport_observables(
	std::map<std::string,
	         std::function<Observable<2>(const std::shared_ptr<SimulateModel>)>>
		&observables) {

	// load all observables here.
	// If the observable is 1D, call Obs2(ObservableCheck(...))
	// If the observable is 2D, call ObservableCheck(...)

	observables["staticconductance"] =
		ObservableCheck(&StaticConductance::StaticG);

	observables["differentialconductance"] =
		ObservableCheck(&DifferentialConductance::DiffG);
}

} // namespace molstat