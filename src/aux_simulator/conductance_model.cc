/**
 * \file conductance_model.cc
 * \brief Function for constructing a ConductanceModel from the name of the
 *    model and a list of named random number distributions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#include "conductance_model.h"
#include "symmetric_one_site.h"
#include "asymmetric_one_site.h"
#include "symmetric_voltage_one_site.h"
#include "symmetric_voltage_two_site.h"

using namespace std;

shared_ptr<ConductanceModel> make_model(const std::string str,
	const std::map<std::string, shared_ptr<RandomDistribution>> &parameters) {

	// create the actual model and process the random number distributions
	if(str == "symmetriconesitemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps;

		// populate the gamma and epsilon distributions
		try {
			dist_gamma = parameters.at("gamma");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gamma\" must be " \
				"specified.");
		}

		try {
			dist_eps = parameters.at("epsilon");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"epsilon\" must be " \
				"specified.");
		}

		return make_shared<SymmetricOneSiteModel>(dist_eps, dist_gamma);
	}
	else if(str == "symmetricvoltageonesitemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps;

		// populate the gamma and epsilon distributions
		try {
			dist_gamma = parameters.at("gamma");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gamma\" must be " \
				"specified.");
		}

		try {
			dist_eps = parameters.at("epsilon");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"epsilon\" must be " \
				"specified.");
		}

		return make_shared<SymmetricVoltageOneSiteModel>(dist_eps, dist_gamma);
	}
	else if(str == "asymmetriconesitemodel") {
		shared_ptr<RandomDistribution> dist_gammaL, dist_gammaR, dist_eps;

		// populate the gamma and epsilon distributions
		try {
			dist_gammaL = parameters.at("gammal");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gammaL\" must be " \
				"specified.");
		}

		try {
			dist_gammaR = parameters.at("gammar");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"gammaR\" must be " \
				"specified.");
		}

		try {
			dist_eps = parameters.at("epsilon");
		}
		catch(const out_of_range &e) {
			throw invalid_argument("A distribution for \"epsilon\" must be " \
				"specified.");
		}

		return make_shared<AsymmetricOneSiteModel>(dist_eps, dist_gammaL,
			dist_gammaR);
	}
	else
		throw invalid_argument("Unrecognized model. Options are:\n" \
			"   SymmetricOneSiteModel - " \
				"Symmetric-Coupling, One-Site Model\n" \
			"   SymmetricVoltageOneSiteModel - " \
				"Symmetric-Coupling, Voltage-Dependent One-Site Model\n" \
			"   AsymmetricOneSiteModel - " \
				"Asymmetric-Coupling, One-Site Model\n");

	// should never be here
	throw invalid_argument("Shouldn't be here.");
	return nullptr;
}
