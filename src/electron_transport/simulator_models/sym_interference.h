/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_interference.h
 * \brief Tight-binding channel with two sites that have an interference 
 *    feature.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#ifndef __sym_interference_h__
#define __sym_interference_h__

#include "observables.h"
#include "junction.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator submodel (molstat::transport::Channel) for transport
 *    through two sites that yield a destructive interference feature.
 *
 * This model uses two degenerate sites with one site coupling to both
 * electrodes.
 *
 * Inherited model parameters (from molstat::transport::TransportJunction) are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias.
 *
 * Submodel parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling,
 * - `beta` (\f$\beta\f$), the inter-site coupling.
 *
 * Starting from \f$ \hat{H} = \left[ \begin{array}{cc} \varepsilon & \beta \\ \beta & \varepsilon \end{array} \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ \begin{array}{cc} -i\Gamma/2 & 0 \\ 0 & 0 \end{array} \right] \f$, the
 * transmission function is
 * \f[ T(E) = \frac{\Gamma^2(E-\varepsilon)^2}{[(E-\varepsilon)^2 - \beta^2] + (E-\varepsilon)^2 \Gamma^2}. \f]
 */
class SymInterferenceChannel : public Channel,
	public ZeroBiasConductance
{
public:
	/// Container index for the Fermi energy.
	static const std::size_t Index_EF;

	/// Container index for the applied bias.
	static const std::size_t Index_V;

	/// Container index for the site energy.
	static const std::size_t Index_epsilon;

	/// Container index for the site-lead coupling.
	static const std::size_t Index_gamma;

	/// Container index for the inter-site coupling.
	static const std::size_t Index_beta;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	virtual ~SymInterferenceChannel() = default;

	/**
	 * \brief Calculates the transmission for a set of model parameters.
	 *
	 * \param[in] e The energy of the incident electron.
	 * \param[in] eps The level energy.
	 * \param[in] gamma The level-electrode coupling.
	 * \param[in] beta The inter-site coupling.
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double eps,
		const double gamma, const double beta);
	
	virtual double ZeroBiasG(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
