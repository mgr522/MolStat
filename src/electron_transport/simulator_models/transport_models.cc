/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file transport_models.cc
 * \brief Functions that load the transport models and observables into
 *    MolStat.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "transport_models.h"
#include "transport_observables.h"
#include "sym_one_site_simulate_model.h"

void load_transport_models(
	std::map<std::string, SimulateModelInstantiator> &models) {

	models["symmetriconesite"] =
		SimulateModelInstance<SymOneSiteSimulateModel>();
}

void load_transport_observables(std::map<std::string, std::function<
	Observable<2>(const shared_ptr<SimulateModel>)>> &observables) {

	// load all observables here.
	// If the observable is 1D, call Obs2(ObservableCheck(...))
	// If the observable is 2D, call ObservableCheck(...)

	observables["staticconductance"] =
		ObservableCheck(&StaticConductance::StaticG);

	observables["differentialconductance"] =
		ObservableCheck(&DifferentialConductance::DiffG);
}
