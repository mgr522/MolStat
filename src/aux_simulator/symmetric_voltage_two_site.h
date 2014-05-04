/**
 * \file symmetric_voltage_two_site.h
 * \brief The symmetric-coupling, voltage-dependent two-site tight-binding
 *        model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __symmetric_voltage_two_site_h__
#define __symmetric_voltage_two_site_h__

#include <memory>
#include "../aux_random_distributions/rng.h"
#include "conductance_model.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the voltage-dependent, two-site model (symmetric
 *    coupling).
 */
class SymmetricVoltageTwoSiteModel : public ConductanceModel {
protected:
	/**
	 * \brief Random distribution for epsilon, the channel energy.
	 */
	shared_ptr<const RandomDistribution> dist_eps;

	/**
	 * \brief Random distribution for gamma, the channel-lead coupling.
	 */
	shared_ptr<const RandomDistribution> dist_gamma;

	/**
	 * \brief Random distribution for beta, the site-site coupling.
	 */
	shared_ptr<const RandomDistribution> dist_beta;

	/**
	 * \brief Calculates the antiderivative needed for the static
	 *    conductance (fixed values of the model parameters).
	 *
	 * \param[in] z The limit of integration.
	 * \param[in] V The voltage.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The antiderivative needed for the static conductance.
	 */
	static double static_c_integral(const double z, const double V,
		const double eps, const double gamma, const double beta);

	/**
	 * \brief Calculates the antiderivative needed for the differential
	 *    conductance (fixed values of the model parameters).
	 *
	 * \param[in] z The limit of integration.
	 * \param[in] V The voltage.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The antiderivative needed for the differential conductance.
	 */
	static double diff_c_integral(const double z, const double V,
		const double gamma, const double beta);

public:
	SymmetricVoltageTwoSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for epsilon.
	 * \param[in] gamma The distribution for gamma.
	 * \param[in] beta The distribution for beta.
	 */
	SymmetricVoltageTwoSiteModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gamma,
		const shared_ptr<const RandomDistribution> &beta);

	/**
	 * \brief Destructor.
	 */
	virtual ~SymmetricVoltageTwoSiteModel() = default;

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
		const double EF, const double eta, const double V) const;

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
		const double EF, const double eta, const double V) const;

	/**
	 * \brief Gets the zero-bias (differential) conductance for a random set of
	 *    model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \return the zero-bias (differential) conductance.
	 */
	virtual double zero_bias_conductance(shared_ptr<gsl_rng> r,
		const double EF) const;

	/**
	 * \brief Calculates the transmission for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] E The incident energy of the electron.
	 * \param[in] V The voltage.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double V, const double eps,
		const double gamma, const double beta);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma,
		const double beta);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma,
		const double beta);
};

#endif
