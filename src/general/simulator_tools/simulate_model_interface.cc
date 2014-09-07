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

using namespace std;

SimulateModel::SimulateModel(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail,
	const std::vector<std::string> &names) {

	int j, n;

	// go through each of the requested parameters
	n = names.size();
	for(j = 0; j < n; ++j) {
		try {
			dists[j] = avail.at(names[j].c_str());
		}
		catch(const out_of_range &e) {
			throw runtime_error(("A distribution for \"" + names[j] + "\" must " \
				"be specified.").c_str());
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
