/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file symmetric_two_site.h
 * \brief The symmetric-coupling, two-site tight-binding model for calculating
 *    conductances.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __symmetric_two_site_h__
#define __symmetric_two_site_h__

#include <memory>
#include <general/random_distributions/rng.h>
#include "conductance_model.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the voltage-dependent, two-site model (symmetric
 *    coupling).
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling,
 * - `beta` (\f$\beta\f$), the inter-site coupling.
 *
 * Starting from \f[ \hat{H} = \left[ \begin{array}{cc} \varepsilon & \beta \\ \beta & \varepsilon \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{L} = \left[ \begin{array}{cc} -i\Gamma/2 & 0 \\ 0 & 0 \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{R} = \left[ \begin{array}{cc} 0 & 0 \\ 0 & -i\Gamma/2 \end{array} \right], \f]
 * the transmission function is
 * \f[ T(E) = \frac{16 \Gamma^2 \beta^2}{\left[ 4(E-\varepsilon)^2-4\beta^2-\Gamma^2\right]^2 + 16 \Gamma^2(E-\varepsilon)^2}. \f]
 * - Differential conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \left[ \eta T(E_\mathrm{F} + \eta eV) + (1-\eta) T(E_\mathrm{F} + (\eta-1)eV) \right]. \f]
 * - Static conductance:
 *   \f{eqnarray*}{ g(V) & = & \frac{2e^2}{h} \frac{2\beta\Gamma}{eV(4\beta^2+e^2V^2)} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon + \eta eV)}{2\beta + i\Gamma} \right) \right] \\
 *   && -\frac{2e^2}{h} \frac{2\beta\Gamma}{eV(4\beta^2+e^2V^2)} \mathrm{Re} \left[ (\Gamma + i2\beta) \mathrm{arctanh}\left( \frac{2(E_\mathrm{F} - \varepsilon + (\eta-1) eV)}{2\beta + i\Gamma} \right) \right]. \f}
 */
class SymmetricTwoSiteModel : public ConductanceModel {
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
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The antiderivative needed for the static conductance.
	 */
	static double static_c_integral(const double z, const double eps,
		const double gamma, const double beta);

public:
	SymmetricTwoSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for epsilon.
	 * \param[in] gamma The distribution for gamma.
	 * \param[in] beta The distribution for beta.
	 */
	SymmetricTwoSiteModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gamma,
		const shared_ptr<const RandomDistribution> &beta);

	/**
	 * \brief Destructor.
	 */
	virtual ~SymmetricTwoSiteModel() = default;

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
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double eps,
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
