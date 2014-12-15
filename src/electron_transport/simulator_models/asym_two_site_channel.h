/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file asym_two_site_channel.h
 * \brief Tight-binding channel with a two-site chain that couples
 *    asymmetrically to both electrodes. The chain does not drop voltage.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __asym_two_site_channel_h__
#define __asym_two_site_channel_h__

#include "observables.h"
#include "junction.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator submodel for transport through a two-site site that couples
 *    asymmetrically to both electrodes.
 *
 * Inherited model parameters (from molstat::transport::TransportJunction) are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias.
 *
 * Model parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gammal` (\f$\Gamma_\mathrm{L}\f$), the site/lead coupling for one
 *   electrode,
 * - `gammar` (\f$\Gamma_\mathrm{R}\f$), the site/lead coupling for the other
 *   electrode,
 * - `beta` (\f$\beta\f$), the inter-site coupling.
 *
 * Starting from \f[ \hat{H} = \left[ \begin{array}{cc} \varepsilon & \beta \\ \beta & \varepsilon \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{L} = \left[ \begin{array}{cc} -i\Gamma_\mathrm{L}/2 & 0 \\ 0 & 0 \end{array} \right], \;\;\;\;
 * \hat{\Sigma}_\mathrm{R} = \left[ \begin{array}{cc} 0 & 0 \\ 0 & -i\Gamma_\mathrm{R}/2 \end{array} \right], \f]
 * the transmission function is
 * \f[ T(E) = \frac{16 \Gamma_\mathrm{L} \Gamma_\mathrm{R} \beta^2}{\left[ 4(E-\varepsilon)^2-4\beta^2-\Gamma_\mathrm{L} \Gamma_\mathrm{R} \right]^2 + 4 (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})^2(E-\varepsilon)^2}. \f]
 * - Differential conductance:
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \frac{1}{2} \left[ T(E_\mathrm{F} + eV/2) + T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Indefinite integral for the static conductance:
 *   \f{eqnarray*}{ \int \mathrm{d}E T(E) & = & \frac{8\sqrt{2}\Gamma_\mathrm{L} \Gamma_\mathrm{R} \beta^2}{(\Gamma_\mathrm{L}+\Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}} \left\{ \frac{\mathrm{arctan}\left[ \frac{\sqrt{8}(E-\varepsilon)}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 - 8\beta^2 - (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}}} \right]}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 -8\beta^2 - (\Gamma_\mathrm{L} + \Gamma_\mathrm{R}) \sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2-16\beta^2}}} \right. \\
 *   && \left. - \frac{\mathrm{arctan}\left[ \frac{\sqrt{8}(E-\varepsilon)}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 - 8\beta^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R})\sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2 - 16\beta^2}}} \right]}{\sqrt{\Gamma_\mathrm{L}^2 + \Gamma_\mathrm{R}^2 -8\beta^2 + (\Gamma_\mathrm{L} + \Gamma_\mathrm{R}) \sqrt{(\Gamma_\mathrm{L} - \Gamma_\mathrm{R})^2-16\beta^2}}} \right\}. \f}
 */
class AsymTwoSiteChannel : public Channel,
	public ElectricCurrent,
	public DifferentialConductance,
	public StaticConductance
{
private:
	/**
	 * \internal
	 * \brief Calculates the antiderivative needed for the electric current
	 *    (fixed values of the model parameters).
	 *
	 * \param[in] z The limit of integration.
	 * \param[in] eps The channel energy, \f$\varepsilon\f$.
	 * \param[in] gammal The left channel-lead coupling,
	 *    \f$\Gamma_\mathrm{L}\f$.
	 * \param[in] gammar The right channel-lead coupling,
	 *    \f$\Gamma_\mathrm{R}\f$.
	 * \param[in] beta The site-site coupling, \f$\beta\f$.
	 * \return The antiderivative needed for the static conductance.
	 * \endinternal
	 */
	static double current_integral(const double z, const double eps,
		const double gammal, const double gammar, const double beta);

public:
	/// Container index for the Fermi energy.
	static const std::size_t Index_EF;

	/// Container index for the applied bias.
	static const std::size_t Index_V;

	/// Container index for the site energy.
	static const std::size_t Index_epsilon;

	/// Container index for the left site-lead coupling.
	static const std::size_t Index_gammaL;

	/// Container index for the right site-lead coupling.
	static const std::size_t Index_gammaR;

	/// Container index for the inter-site coupling.
	static const std::size_t Index_beta;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	virtual ~AsymTwoSiteChannel() = default;

	/**
	 * \brief Calculates the transmission for a set of model parameters.
	 *
	 * \param[in] e The energy of the incident electron.
	 * \param[in] v The applied bias.
	 * \param[in] eps The level energy.
	 * \param[in] gammal The left level-electrode coupling.
	 * \param[in] gammar The right level-electrode coupling.
	 * \param[in] beta The site-site coupling.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double v, const double eps,
		const double gammal, const double gammar, const double beta);

	double ECurrent(const std::valarray<double> &params) const override;
	double StaticG(const std::valarray<double> &params) const override;
	double DiffG(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
