/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
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
#include "symmetric_two_site.h"
#include "asymmetric_two_site.h"
#include "symmetric_voltage_two_site.h"

using namespace std;

/**
 * \brief Gets a distribution with the desired name from the map, throwing an
 *    exception if none is found.
 *
 * \throw invalid_argument if there is not a distribution with the name `name`.
 *
 * \param[in] name The name of the distribution.
 * \param[in] parameters The map of distributions.
 * \return The distribution.
 */
static shared_ptr<RandomDistribution> find_distribution(const string name,
	const map<string, shared_ptr<RandomDistribution>> &parameters) {

	shared_ptr<RandomDistribution> ret;

	try {
		ret = parameters.at(name.c_str());
	}
	catch(const out_of_range &e) {
		throw invalid_argument(("A distribution for \"" + name + "\" must be " \
			"specified.").c_str());
	}

	return ret;
}

shared_ptr<ConductanceModel> make_model(const std::string str,
	const std::map<std::string, shared_ptr<RandomDistribution>> &parameters) {

	// create the actual model and process the random number distributions
	if(str == "symmetriconesitemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps;

		dist_gamma = find_distribution("gamma", parameters);
		dist_eps = find_distribution("epsilon", parameters);

		return make_shared<SymmetricOneSiteModel>(dist_eps, dist_gamma);
	}
	else if(str == "asymmetriconesitemodel") {
		shared_ptr<RandomDistribution> dist_gammaL, dist_gammaR, dist_eps;

		dist_gammaL = find_distribution("gammal", parameters);
		dist_gammaR = find_distribution("gammar", parameters);
		dist_eps = find_distribution("epsilon", parameters);

		return make_shared<AsymmetricOneSiteModel>(dist_eps, dist_gammaL,
			dist_gammaR);
	}
	else if(str == "symmetricvoltageonesitemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps;

		dist_gamma = find_distribution("gamma", parameters);
		dist_eps = find_distribution("epsilon", parameters);

		return make_shared<SymmetricVoltageOneSiteModel>(dist_eps, dist_gamma);
	}
	else if(str == "symmetrictwositemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps, dist_beta;

		dist_gamma = find_distribution("gamma", parameters);
		dist_eps = find_distribution("epsilon", parameters);
		dist_beta = find_distribution("beta", parameters);

		return make_shared<SymmetricTwoSiteModel>(dist_eps, dist_gamma,
			dist_beta);
	}
	else if(str == "asymmetrictwositemodel") {
		shared_ptr<RandomDistribution> dist_gammaL, dist_gammaR, dist_eps,
			dist_beta;

		dist_gammaL = find_distribution("gammal", parameters);
		dist_gammaR = find_distribution("gammar", parameters);
		dist_eps = find_distribution("epsilon", parameters);
		dist_beta = find_distribution("beta", parameters);

		return make_shared<AsymmetricTwoSiteModel>(dist_eps, dist_gammaL,
			dist_gammaR, dist_beta);
	}
	else if(str == "symmetricvoltagetwositemodel") {
		shared_ptr<RandomDistribution> dist_gamma, dist_eps, dist_beta;

		dist_gamma = find_distribution("gamma", parameters);
		dist_eps = find_distribution("epsilon", parameters);
		dist_beta = find_distribution("beta", parameters);

		return make_shared<SymmetricVoltageTwoSiteModel>(dist_eps, dist_gamma,
			dist_beta);
	}
	else
		throw invalid_argument("Unrecognized model. Options are:\n" \
			"   SymmetricOneSiteModel - " \
				"Symmetric-Coupling, One-Site Model\n" \

			"   AsymmetricOneSiteModel - " \
				"Asymmetric-Coupling, One-Site Model\n" \

			"   SymmetricVoltageOneSiteModel - " \
				"Symmetric-Coupling, Voltage-Dependent One-Site Model\n" \

			"   SymmetricTwoSiteModel - " \
				"Symmetric-Coupling, Two-Site Model\n" \

			"   AsymmetricTwoSiteModel - " \
				"Asymmetric-Coupling, Two-Site Model\n" \

			"   SymmetricVoltageTwoSiteModel - " \
				"Symmetric-Coupling, Voltage-Dependent Two-Site Model\n");

	// should never be here
	throw invalid_argument("Shouldn't be here.");
	return nullptr;
}
