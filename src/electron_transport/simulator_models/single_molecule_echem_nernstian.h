/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file single_molecule_echem_nernstian.h
 * \brief Nernstian limit of single molecule electrochemistry model.
 *
 * \author Bo Fu
 * \date December 2014
 */

#ifndef __single_molecule_echem_nernstian_h__
#define __single_molecule_echem_nernstian_h__

#include "observables.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator model (molstat::transport) for single molecule electrochemistry
 *    under Nernstian limit with a slow enough scanning rate.
 *
 * Inherited model parameters (from molstat::transport::TransportJunction) are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias.
 *
 * Submodel parameters are
 * - `epsilon` (\f$\varepsilon\f$), the site-energy,
 * - `gamma` (\f$\Gamma\f$), the site/lead coupling,
 * - `a` (\f$a\f$), the scaling factor for the voltage drop.
 *
 * Starting from \f$ \hat{H} = \left[ \varepsilon + aeV \right] \f$ and
 * \f$ \hat{\Sigma}_\mathrm{L/R} = \left[ -i\Gamma/2 \right] \f$, the
 * transmission function is
 * \f[ T(E) = \frac{\Gamma^2}{(E-\varepsilon-aeV)^2 + \Gamma^2}. \f]
 * - Differential conductance (with intermediate quantities):
 *   \f[ \frac{\partial}{\partial V}T(E) = \frac{2ea\Gamma^2(E - \varepsilon-aeV)}{[(E-\varepsilon-aeV)^2+\Gamma^2]^2}; \f]
 *   \f[ \frac{2e}{h} \int\limits_{E_\mathrm{F}-eV/2}^{E_\mathrm{F}+eV/2} \mathrm{d}E \frac{\partial}{\partial V} T(E) = \frac{2e^2a}{h} T(E_\mathrm{F}-eV/2) - \frac{2e^2a}{h}T(E_\mathrm{F}+eV/2); \f]
 *   \f[ G_\mathrm{d}(V) = \frac{2e^2}{h} \left[ (1/2-a) T(E_\mathrm{F} + eV/2) + (1/2+a) T(E_\mathrm{F} - eV/2) \right]. \f]
 * - Static conductance:
 *   \f[ G_\mathrm{s}(V) = \frac{2e^2}{h} \frac{\Gamma}{eV} \left[ \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (1/2-a) eV}{\Gamma} \right) - \arctan\left( \frac{E_\mathrm{F} - \varepsilon - (1/2+a) eV}{\Gamma} \right) \right]. \f]
 */
class SingleMoleculeEchemNernstian : 
	public PeakPotential,
	public virtual molstat::SimulateModel
{
public:
	/// Container index for the reference potential.
	static const std::size_t Index_E0;

	/// Container index for the prefactor of forward half-reaction rate constant.
	static const std::size_t Index_Af;

	/// Container index for the prefactor of backward half-reaction rate constant.

	static const std::size_t Index_Ab;

	/// Container index for the temperature.
	static const std::size_t Index_T;

	/// Container index for the number of electrons involved in the half-reaction.
	static const std::size_t Index_n;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	virtual ~SingleMoleculeEchemNernstian() = default;

	virtual double PeakV(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
