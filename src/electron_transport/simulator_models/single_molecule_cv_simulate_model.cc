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

#include "single_molecule_cv_simulate_model.h"
#include <cmath>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>

using namespace std;

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const vector<string> SingleMoleculeCV::parameters =
	{"e0", "eref", "lambda", "af", "ab", "v", "n", "poinitial", "temperature","tlimit"};

void SingleMoleculeCV::unpack_parameters(const std::vector<double> &vec,
	double &e0, double &eref, double &lambda, double &af, double &ab,
	double &v, double &n, double &poinitial, double &temperature, double &tlimit) {

	e0 = vec[0];
	eref = vec[1];
	lambda = vec[2];
	af = vec[3];
	ab = vec[4];
	v = vec[5];
    n = vec[6];
    poinitial = vec[7];
    temperature = vec[8];
    tlimit = vec[9];
}

SingleMoleculeCV::SingleMoleculeCV(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}

std::array<double, 2> SingleMoleculeCV::DiffG(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(10);
	double ef, v, eps, gammal, gammar, a;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gammal, gammar, a);

	return {{ v, diff_conductance(params) }};
}

std::array<double, 2> SingleMoleculeCV::PeakPotentials(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(10);
	double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);

	return {{ v, static_conductance(params) }};/*Modification needed later for new function.*/
}

double SingleMoleculeCV::kf(const std::vector<double> &vecl) {

    double e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit;

    // umpack the model parameters
    unpack_parameters(params, e0, eref, lambda, af, ab, v, n, poinitial, temperature, tlimit);
}
double SingleMoleculeCV::static_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gammal, gammar, a;

	// unpack the model parameters
	unpack_parameters(vec, ef, v, eps, gammal, gammar, a);

	return 2.*gammal*gammar / (v*(gammal + gammar)) *
		(atan(2. * (ef-eps+(0.5-a)*v) / (gammal + gammar))
		- atan(2. * (ef-eps-(0.5+a)*v) / (gammal + gammar)));
}

double SingleMoleculeCV::diff_conductance(
	const std::vector<double> &vec) {

	double ef, v, eps, gammal, gammar, a;

	// unpack the parameters
	unpack_parameters(vec, ef, v, eps, gammal, gammar, a);

	return (0.5 - a) * transmission(ef+0.5*v, v, eps, gammal, gammar, a) +
		(0.5 + a) * transmission(ef-0.5*v, v, eps, gammal, gammar, a);
}
