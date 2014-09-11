/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file ind_two_chan_simulate_model.cc
 * \brief Generic model for two independent channels that each couple
 *    symmetrically to both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "ind_two_chan_simulate_model.h"
#include "sym_one_site_simulate_model.h"

using namespace std;

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const vector<string> IndTwoChanSimulateModel::parameters =
	{"ef", "v", "epsilon1", "gamma1", "a1", "epsilon2", "gamma2", "a2"};

void IndTwoChanSimulateModel::unpack_parameters(const std::vector<double> &vec,
	double &ef, double &v, double &epsilon1, double &gamma1, double &a1,
	double &epsilon2, double &gamma2, double &a2) {

	ef = vec[0];
	v = vec[1];
	epsilon1 = vec[2];
	gamma1 = vec[3];
	a1 = vec[4];
	epsilon2 = vec[5];
	gamma2 = vec[6];
	a2 = vec[7];
}

IndTwoChanSimulateModel::IndTwoChanSimulateModel(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}

std::array<double, 2> IndTwoChanSimulateModel::DiffG(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(8);
	double ef, v, eps1, gamma1, a1, eps2, gamma2, a2;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps1, gamma1, a1, eps2, gamma2, a2);

	return {{ v, diff_conductance(params) }};
}

std::array<double, 2> IndTwoChanSimulateModel::StaticG(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(8);
	double ef, v, eps1, gamma1, a1, eps2, gamma2, a2;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps1, gamma1, a1, eps2, gamma2, a2);

	return {{ v, static_conductance(params) }};
}

double IndTwoChanSimulateModel::transmission(const double e, const double v,
	const double eps1, const double gamma1, const double a1, const double eps2,
	const double gamma2, const double a2) {

	return SymOneSiteSimulateModel::transmission(e, v, eps1, gamma1, a1) +
		SymOneSiteSimulateModel::transmission(e, v, eps2, gamma2, a2);
}

double IndTwoChanSimulateModel::static_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps1, gamma1, a1, eps2, gamma2, a2;
	vector<double> vec1(5), vec2(5);

	// unpack the model parameters
	unpack_parameters(vec, ef, v, eps1, gamma1, a1, eps2, gamma2, a2);

	vec1[0] = vec2[0] = ef;
	vec1[1] = vec2[1] = v;
	vec1[2] = eps1;
	vec1[3] = gamma1;
	vec1[4] = a1;
	vec2[2] = eps2;
	vec2[3] = gamma2;
	vec2[4] = a2;

	return SymOneSiteSimulateModel::static_conductance(vec1) +
		SymOneSiteSimulateModel::static_conductance(vec2);
}

double IndTwoChanSimulateModel::diff_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps1, gamma1, a1, eps2, gamma2, a2;
	vector<double> vec1(5), vec2(5);

	// unpack the parameters
	unpack_parameters(vec, ef, v, eps1, gamma1, a1, eps2, gamma2, a2);

	vec1[0] = vec2[0] = ef;
	vec1[1] = vec2[1] = v;
	vec1[2] = eps1;
	vec1[3] = gamma1;
	vec1[4] = a1;
	vec2[2] = eps2;
	vec2[3] = gamma2;
	vec2[4] = a2;

	return SymOneSiteSimulateModel::diff_conductance(vec1) +
		SymOneSiteSimulateModel::diff_conductance(vec2);
}
