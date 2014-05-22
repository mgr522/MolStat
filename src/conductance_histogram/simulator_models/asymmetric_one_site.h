/**
 * \file asymmetric_one_site.h
 * \brief The asymmetric-coupling, tight-binding model for calculating
 *        conductances.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __asymmetric_one_site_h__
#define __asymmetric_one_site_h__

#include <memory>
#include <general/random_distributions/rng.h>
#include "conductance_model.h"

using std::shared_ptr;

/**
 * \brief Class encapsulating the asymmetric-coupling, one-site tight-binding
 *    model.
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gammaL` (\f$\Gamma_\mathrm{L}\f$), the site/lead coupling for one electrode.
 * - `gammaR` (\f$\Gamma_\mathrm{R}\f$), the site/lead coupling for the other electrode.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma_\mathrm{L,R}/2 \right] \f$,
 * the transmission function is
 * \f[ T(E) = \frac{4\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{4(E-\varepsilon)^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})^2}. \f]
 * - Differential conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \left[ \eta T(E_\mathrm{F} + \eta eV) + (1-\eta) T(E_\mathrm{F} + (\eta-1)eV) \right]. \f]
 * - Static conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \frac{2\Gamma_\mathrm{L} \Gamma_\mathrm{R}}{eV(\Gamma_\mathrm{L} +\Gamma_\mathrm{R})} \left[ \arctan\left( \frac{2(E_\mathrm{F} - \varepsilon + \eta eV)}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) - \arctan\left( \frac{2(E_\mathrm{F} - \varepsilon + (\eta-1)eV)}{\Gamma_\mathrm{L} + \Gamma_\mathrm{R}} \right) \right]. \f]
 */
class AsymmetricOneSiteModel : public ConductanceModel {
protected:
	/**
	 * \brief Random distribution for epsilon, the channel energy.
	 */
	shared_ptr<const RandomDistribution> dist_eps;

	/**
	 * \brief Random distribution for gammaL, one channel-lead coupling.
	 */
	shared_ptr<const RandomDistribution> dist_gammaL;

	/**
	 * \brief Random distribution for gammaR, one channel-lead coupling.
	 */
	shared_ptr<const RandomDistribution> dist_gammaR;

public:
	AsymmetricOneSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for epsilon.
	 * \param[in] gammaL The distribution for gammaL.
	 * \param[in] gammaR The distribution for gammaR.
	 */
	AsymmetricOneSiteModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gammaL,
		const shared_ptr<const RandomDistribution> &gammaR);

	/**
	 * \brief Destructor.
	 */
	virtual ~AsymmetricOneSiteModel() = default;

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
	 *        model parameters.
	 *
	 * \param[in] r The handle for GSL random number generation.
	 * \param[in] EF The Fermi energy.
	 * \return The differential conductance.
	 */
	virtual double zero_bias_conductance(shared_ptr<gsl_rng> r,
		const double EF) const;

	/**
	 * \brief Calculates the transmission for fixed values of epsilon and the.
	 *     gammas.
	 *
	 * \param[in] E The incident energy of the electron.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL The left channel-lead coupling, gammaL.
	 * \param[in] gammaR The right channel-lead coupling, gammaR.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double eps,
		const double gammaL, const double gammaR);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL The left channel-lead coupling, gammaL.
	 * \param[in] gammaR The right channel-lead coupling, gammaR.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL The left channel-lead coupling, gammaL.
	 * \param[in] gammaR The right channel-lead coupling, gammaR.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR);
};

#endif
