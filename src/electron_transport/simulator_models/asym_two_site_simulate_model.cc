/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_two_site_simulate_model.cc
 * \brief Tight-binding model of a two-site chain that couples asymmetrically
 *    to both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#include "asym_two_site_simulate_model.h"
#include <cmath>
#include <complex>

using namespace std;

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const vector<string> AsymTwoSiteSimulateModel::parameters =
	{"ef", "v", "epsilon", "gammal", "gammar", "beta"};

void AsymTwoSiteSimulateModel::unpack_parameters(
	const std::vector<double> &vec, double &ef, double &v, double &epsilon,
	double &gammal, double &gammar, double &beta) {

	ef = vec[0];
	v = vec[1];
	epsilon = vec[2];
	gammal = vec[3];
	gammar = vec[4];
	beta = vec[5];
}

AsymTwoSiteSimulateModel::AsymTwoSiteSimulateModel(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}

std::array<double, 2> AsymTwoSiteSimulateModel::DiffG(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(6);
	double ef, v, eps, gammal, gammar, beta;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gammal, gammar, beta);

	return {{ v, diff_conductance(params) }};
}

std::array<double, 2> AsymTwoSiteSimulateModel::StaticG(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(6);
	double ef, v, eps, gammal, gammar, beta;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gammal, gammar, beta);

	return {{ v, static_conductance(params) }};
}

double AsymTwoSiteSimulateModel::transmission(const double e, const double v,
	const double eps, const double gammal, const double gammar,
	const double beta) {

	double temp = 4.*(e-eps)*(e-eps) - 4.*beta*beta - gammal*gammar;

	return 16.*gammal*gammar*beta*beta /
		(temp*temp + 4.*(gammal+gammar)*(gammal+gammar)*(e-eps)*(e-eps));
}

double AsymTwoSiteSimulateModel::static_c_integral(const double z,
	const double eps, const double gammal, const double gammar,
	const double beta) {

	complex<double> bgg = sqrt(complex<double>((gammal-gammar)*(gammal-gammar)
		- 16.*beta*beta, 0.));

	complex<double> denom1 = sqrt(-8.*beta*beta + gammal*gammal + gammar*gammar
		- (gammal + gammar) * bgg);

	complex<double> denom2 = sqrt(-8.*beta*beta + gammal*gammal + gammar*gammar
		+ (gammal + gammar) * bgg);

	return sqrt(128.)*gammal*gammar*beta*beta / (gammal + gammar) * real(
		(atan(sqrt(8.)*(z-eps) / denom1) / denom1 -
		 atan(sqrt(8.)*(z-eps) / denom2) / denom2) / bgg
		);
}

double AsymTwoSiteSimulateModel::static_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gammal, gammar, beta;

	// unpack the model parameters
	unpack_parameters(vec, ef, v, eps, gammal, gammar, beta);

	return (static_c_integral(ef + 0.5*v, eps, gammal, gammar, beta) -
		static_c_integral(ef - 0.5*v, eps, gammal, gammar, beta)) / v;
}

double AsymTwoSiteSimulateModel::diff_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gammal, gammar, beta;

	// unpack the parameters
	unpack_parameters(vec, ef, v, eps, gammal, gammar, beta);

	return 0.5*transmission(ef + 0.5*v, v, eps, gammal, gammar, beta) +
		0.5*transmission(ef - 0.5*v, v, eps, gammal, gammar, beta);
}
