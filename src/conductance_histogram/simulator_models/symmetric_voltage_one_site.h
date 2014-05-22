/**
 * \file symmetric_voltage_one_site.h
 * \brief The symmetric-coupling, voltage-dependent one-site tight-binding
 *        model for calculating conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#ifndef __symmetric_voltage_one_site_h__
#define __symmetric_voltage_one_site_h__

#include <memory>
#include <general/random_distributions/rng.h>
#include "conductance_model.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the voltage-dependent, one-site model (symmetric
 *    coupling).
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon + eV \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma/2 \right] \f$, the
 * transmission function is
 * \f[ T(E) = \frac{\Gamma^2}{(E-\varepsilon-eV)^2 + \Gamma^2}. \f]
 * - Differential conductance (with intermediate quantities):
 *   \f[ \frac{\partial}{\partial V}T(E) = \frac{2e\Gamma^2(E - \varepsilon-eV)}{[(E-\varepsilon-eV)^2+\Gamma^2]^2}; \f]
 *   \f[ \frac{2e}{h} \int\limits_{E_\mathrm{F}+(\eta-1)eV}^{E_\mathrm{F}+\eta eV} \mathrm{d}E \frac{\partial}{\partial V} T(E) = \frac{2e^2}{h} T(E_\mathrm{F}+(\eta-1)eV) - \frac{2e^2}{h}T(E_\mathrm{F}+\eta eV); \f]
 *   \f[ g(V) = \frac{2e^2}{h} \left[ (\eta-1) T(E_\mathrm{F} + \eta eV) + (2-\eta) T(E_\mathrm{F} + (\eta-1)eV) \right]. \f]
 * - Static conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \frac{\Gamma}{eV} \left[ \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (\eta-1) eV}{\Gamma} \right) - \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (\eta-2) eV}{\Gamma} \right) \right]. \f]
 */
class SymmetricVoltageOneSiteModel : public ConductanceModel {
protected:
	/**
	 * \brief Random distribution for epsilon, the channel energy.
	 */
	shared_ptr<const RandomDistribution> dist_eps;

	/**
	 * \brief Random distribution for gamma, the channel-lead coupling.
	 */
	shared_ptr<const RandomDistribution> dist_gamma;

public:
	SymmetricVoltageOneSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for epsilon.
	 * \param[in] gamma The distribution for gamma.
	 */
	SymmetricVoltageOneSiteModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gamma);

	/**
	 * \brief Destructor.
	 */
	virtual ~SymmetricVoltageOneSiteModel() = default;

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
	 * \brief Calculates the transmission for fixed values of epsilon and gamma.
	 *
	 * \param[in] E The incident energy of the electron.
	 * \param[in] V The voltage.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double V, const double eps,
		const double gamma);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gamma The channel-lead coupling, gamma.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gamma);
};

#endif
