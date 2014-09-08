/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_one_site_simulate_model.cc
 * \brief Tight-binding model with one site that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "sym_one_site_simulate_model.h"

using namespace std;

SymOneSiteSimulateModel::SymOneSiteSimulateModel(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, { "ef", "v" }) {
}

std::array<double, 1> SymOneSiteSimulateModel::ZeroBiasG(shared_ptr<gsl_rng> r)
	const {

	return {0.};
}
