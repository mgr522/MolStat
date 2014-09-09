/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface.cc
 * \brief Implementation of an abstract class for simulating histograms.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "simulate_model_interface.h"
#include <array>

using namespace std;

SimulateModel::SimulateModel(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail,
	const std::vector<std::string> &names)
	: dists(names.size()) {

	int j, n;

	// go through each of the requested parameters
	n = names.size();
	for(j = 0; j < n; ++j) {
		try {
			dists[j] = avail.at(names[j].c_str());
		}
		catch(const out_of_range &e) {
			throw runtime_error(names[j].c_str());
		}
	}
}

void SimulateModel::sample(shared_ptr<gsl_rng> r, vector<double> &vals) const {
	int j, n;

	// sample each distribution, in order
	n = dists.size();
	for(j = 0; j < n; ++j)
		vals[j] = dists[j]->sample(r);
}

std::function<Observable<2>(const shared_ptr<SimulateModel>)> Obs2(
	const std::function<Observable<1>(const shared_ptr<SimulateModel>)> &f) {

	return [=] (const shared_ptr<SimulateModel> model) -> Observable<2> {
		Observable<1> obs1 = f(model);

		return [=] (shared_ptr<gsl_rng> r) -> array<double, 2> {
			return { { 0., obs1(r)[0] } };
		};
	};
}
