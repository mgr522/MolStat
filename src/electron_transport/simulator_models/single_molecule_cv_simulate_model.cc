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

using namespace std;

// if the order of the following list is changed, the unpack_parameters
// function MUST also be updated
const vector<string> SingleMoleculeCV::parameters =
	{"ef", "v", "epsilon", "gammal", "gammar", "a", "b"};

void SingleMoleculeCV::unpack_parameters(const std::vector<double> &vec,
	double &ef, double &v, double &epsilon, double &gammal, double &gammar,
	double &a, double &b) {

	ef = vec[0];
	v = vec[1];
	epsilon = vec[2];
	gammal = vec[3];
	gammar = vec[4];
	a = vec[5];
    b = vec[6];
}

SingleMoleculeCV::SingleMoleculeCV(
	const std::map<std::string, shared_ptr<RandomDistribution>> &avail)
	: SimulateModel(avail, parameters) {
}



std::array<double, 2> SingleMoleculeCV::PeakPotentials(shared_ptr<gsl_rng> r)
	const {

	vector<double> params(7);
	double ef, v, eps, gammal, gammar, a, b;

	// generate and unpack the parameters
	sample(r, params);
	unpack_parameters(params, ef, v, eps, gammal, gammar, a, b);

	return {{ v, E_applied(0, params) }};
}



double SingleMoleculeCV::kf( double t,
	const std::vector<double> &vec) {

	double ef, v, eps, gammal, gammar, a, b;

	// unpack the model parameters
	unpack_parameters(vec, ef, v, eps, gammal, gammar, a, b);

	return 2.*gammal*gammar / (v*(gammal + gammar)) *
		(atan(2. * (ef-eps+(0.5-a)*v) / (gammal + gammar))
		- atan(2. * (ef-eps-(0.5+a)*v) / (gammal + gammar)));
}



double SingleMoleculeCV::E_applied(double t,
    const std::vector<double> &vec) {

    double ef, v, eps, gammal, gammar, a, b;

    //upack the model paramters
    unpack_parameters(vec, ef, v, eps, gammal, gammar, a, b);

    double E;
    double e0 = 1.0;
    double tlimit = 1.0;
    
    if (t >= 0 && t <= tlimit)
        E = e0 + v * t;
    if (t > tlimit && t <= 2.0 * tlimit)
        E = e0 + 2.0 * v * tlimit - v * t;
    return E;
}
