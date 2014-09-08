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
#include "sym_one_site_simulate_model.h"

using namespace std;

/**
 * \brief Creates the wrapper function for instantiating a model.
 *
 * \return The wrapper to a function for instantiating the model.
 */
template<typename T>
inline function<shared_ptr<SimulateModel>
	(const map<string, shared_ptr<RandomDistribution>> &)> model_maker() {

	return [] (const map<string, shared_ptr<RandomDistribution>> &avail)
		-> shared_ptr<SimulateModel> {

		return make_shared<T>(avail);
	};
}

void load_transport_models(std::map<std::string,
	std::function<shared_ptr<SimulateModel>(
		const std::map<std::string, shared_ptr<RandomDistribution>>&)>>
	&models) {

	models["symmetriconesite"] = model_maker<SymOneSiteSimulateModel>();
}
