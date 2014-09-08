/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file conductance_model.h
 * \brief Interface for the various tight-binding models and their transmission
 *        functions/conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#ifndef __conductance_model_h__
#define __conductance_model_h__

#include <memory>
#include <string>
#include <map>
#include <general/random_distributions/rng.h>

using std::shared_ptr;

/**
 * \brief Abstract base class a conductance model.
 */
class ConductanceModel {
public:
	/**
	 * \internal
	 * \brief Default constructor.
	 * \endinternal
	 */
	ConductanceModel() = default;

	/**
	 * \internal
	 * \brief Destructor.
	 * \endinternal
	 */
	virtual ~ConductanceModel() = default;

	/**
	 * \brief Gets the static conductance for a random set of model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \param[in] eta The relative voltage drop.
	 * \param[in] V The voltage.
	 * \return The static conductance.
	 */
	virtual double static_conductance(shared_ptr<gsl_rng> r,
		const double EF, const double eta, const double V) const = 0;

	/**
	 * \brief Gets the differential conductance for a random set of model
	 *        parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \param[in] eta The relative voltage drop.
	 * \param[in] V The voltage.
	 * \return The differential conductance.
	 */
	virtual double diff_conductance(shared_ptr<gsl_rng> r,
		const double EF, const double eta, const double V) const = 0;

	/**
	 * \brief Gets the zero bias (differential) conductance for a random set of
	 *        model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \return The zero-bias (differential) conductance.
	 */
	virtual double zero_bias_conductance(shared_ptr<gsl_rng> r,
		const double EF) const = 0;
};

/**
 * \brief Uses the distributions from the input deck to constructs a model of
 *    the specified type.
 *
 * \throw invalid_argument if any required distributions are not specified.
 *
 * \param[in] str The name of the model type, assumed to be lowercase.
 * \param[in] parameters Map of distribution names (lowercase) to the
 *    distributions.
 * \return A shared pointer to the model object.
 */
shared_ptr<ConductanceModel> make_model(const std::string str,
	const std::map<std::string, shared_ptr<RandomDistribution>> &parameters);

#endif
