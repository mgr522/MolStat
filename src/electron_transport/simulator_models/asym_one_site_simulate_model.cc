/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_one_site_simulate_model.cc
 * \brief Tight-binding model with one site that couples asymmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "asym_one_site_simulate_model.h"
#include <cmath>

namespace molstat {

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const std::vector<std::string> AsymOneSiteSimulateModel::parameters =
	{"ef", "v", "epsilon", "gammal", "gammar", "a"};

void AsymOneSiteSimulateModel::unpack_parameters(const std::vector<double> &vec,
	double &ef, double &v, double &epsilon, double &gammal, double &gammar,
	double &a) {

	ef = vec[0];
	v = vec[1];
	epsilon = vec[2];
	gammal = vec[3];
	gammar = vec[4];
	a = vec[5];
}

AsymOneSiteSimulateModel::AsymOneSiteSimulateModel(
	const std::map<std::string, std::shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}

std::array<double, 2> AsymOneSiteSimulateModel::DiffG(gsl_rng_ptr &r) const {
	std::vector<double> params(6);
	double ef, v, eps, gammal, gammar, a;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gammal, gammar, a);

	return {{ v, diff_conductance(params) }};
}

std::array<double, 2> AsymOneSiteSimulateModel::StaticG(gsl_rng_ptr &r) const {
	std::vector<double> params(6);
	double ef, v, eps, gammal, gammar, a;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gammal, gammar, a);

	return {{ v, static_conductance(params) }};
}

double AsymOneSiteSimulateModel::transmission(const double e, const double v,
	const double eps, const double gammal, const double gammar,
	const double a) {

	return 4.*gammal*gammar / (4.*(e - eps - a*v)*(e - eps - a*v) +
		(gammal + gammar)*(gammal + gammar));
}

double AsymOneSiteSimulateModel::static_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gammal, gammar, a;

	// unpack the model parameters
	unpack_parameters(vec, ef, v, eps, gammal, gammar, a);

	return 2.*gammal*gammar / (v*(gammal + gammar)) *
		(atan(2. * (ef-eps+(0.5-a)*v) / (gammal + gammar))
		- atan(2. * (ef-eps-(0.5+a)*v) / (gammal + gammar)));
}

double AsymOneSiteSimulateModel::diff_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gammal, gammar, a;

	// unpack the parameters
	unpack_parameters(vec, ef, v, eps, gammal, gammar, a);

	return (0.5 - a) * transmission(ef+0.5*v, v, eps, gammal, gammar, a) +
		(0.5 + a) * transmission(ef-0.5*v, v, eps, gammal, gammar, a);
}

} // namespace molstat
