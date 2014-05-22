/**
 * \file asymmetric_two_site.h
 * \brief The asymmetric-coupling, two-site tight-binding model for calculating
 *    conductances.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __asymmetric_two_site_h__
#define __asymmetric_two_site_h__

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
 * - `gammaL` (\f$\Gamma_\mathrm{L}\f$), the site/lead coupling for one
 *   electrode,
 * - `gammaR` (\f$\Gamma_\mathrm{R}\f$), the site/lead coupling for the other
 *   electrode,
 * - `beta` (\f$\beta\f$), the inter-site coupling.
 *
 * Starting from \f[ \hat{H} = \left[ \begin{array}{cc} \varepsilon & \beta \\ \beta & \varepsilon \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{L} = \left[ \begin{array}{cc} -i\Gamma_\mathrm{L}/2 & 0 \\ 0 & 0 \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{R} = \left[ \begin{array}{cc} 0 & 0 \\ 0 & -i\Gamma_\mathrm{R}/2 \end{array} \right], \f]
 * the transmission function is
 * \f[ T(E) = \frac{16 \Gamma_\mathrm{L} \Gamma_\mathrm{R} \beta^2}{\left[ 4(E-\varepsilon)^2-4\beta^2-\Gamma_\mathrm{L} \Gamma_\mathrm{R} \right]^2 + 4 (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})^2(E-\varepsilon)^2}. \f]
 * - Differential conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \left[ \eta T(E_\mathrm{F} + \eta eV) + (1-\eta) T(E_\mathrm{F} + (\eta-1)eV) \right]. \f]
 * - Indefinite integral for the static conductance:
 *   \f{eqnarray*}{ \int \mathrm{d}E T(E) & = & \frac{8\sqrt{2}\Gamma_\mathrm{L} \Gamma_\mathrm{R} \beta^2}{(\Gamma_\mathrm{L}+\Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}} \left\{ \frac{\mathrm{arctan}\left[ \frac{\sqrt{8}(E-\varepsilon)}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 - 8\beta^2 - (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}}} \right]}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 -8\beta^2 - (\Gamma_\mathrm{L} + \Gamma_\mathrm{R}) \sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2-16\beta^2}}} \right. \\
 *   && \left. - \frac{\mathrm{arctan}\left[ \frac{\sqrt{8}(E-\varepsilon)}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 - 8\beta^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}}} \right]}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 -8\beta^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R}) \sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2-16\beta^2}}} \right\}. \f}
 */
class AsymmetricTwoSiteModel : public ConductanceModel {
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
	 * \param[in] gammaL One channel-lead coupling, gammaL.
	 * \param[in] gammaR One channel-lead coupling, gammaR.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The antiderivative needed for the static conductance.
	 */
	static double static_c_integral(const double z, const double eps,
		const double gammaL, const double gammaR, const double beta);

public:
	AsymmetricTwoSiteModel() = delete;

	/**
	 * \brief Constructor requiring the random distributions.
	 *
	 * \param[in] eps The distribution for epsilon.
	 * \param[in] gammaL The distribution for gammaL.
	 * \param[in] gammaR The distribution for gammaR.
	 * \param[in] beta The distribution for beta.
	 */
	AsymmetricTwoSiteModel(
		const shared_ptr<const RandomDistribution> &eps,
		const shared_ptr<const RandomDistribution> &gammaL,
		const shared_ptr<const RandomDistribution> &gammaR,
		const shared_ptr<const RandomDistribution> &beta);

	/**
	 * \brief Destructor.
	 */
	virtual ~AsymmetricTwoSiteModel() = default;

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
	 * \param[in] gammaL One channel-lead coupling, gammaL.
	 * \param[in] gammaR One channel-lead coupling, gammaR.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The transmission.
	 */
	static double transmission(const double E, const double eps,
		const double gammaL, const double gammaR, const double beta);

	/**
	 * \brief Calculates the static conductance for fixed values of the model
	 *    parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL One channel-lead coupling, gammaL.
	 * \param[in] gammaR One channel-lead coupling, gammaR.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The static conductance.
	 */
	static double static_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR, const double beta);

	/**
	 * \brief Calculates the differential conductance for fixed values of the
	 *    model parameters.
	 *
	 * \param[in] EF The Fermi energy.
	 * \param[in] V The voltage.
	 * \param[in] eta The relative voltage drops at the leads.
	 * \param[in] eps The channel energy, epsilon.
	 * \param[in] gammaL One channel-lead coupling, gammaL.
	 * \param[in] gammaR One channel-lead coupling, gammaR.
	 * \param[in] beta The site-site coupling, beta.
	 * \return The differential conductance.
	 */
	static double diff_conductance(const double EF, const double V,
		const double eta, const double eps, const double gammaL,
		const double gammaR, const double beta);
};

#endif
