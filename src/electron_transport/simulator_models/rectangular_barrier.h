/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file rectangular_barrier.h
 * \brief Channel with transmission given by a rectangular barrier.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#ifndef __rectangular_barrier_h__
#define __rectangular_barrier_h__

#include "observables.h"
#include "junction.h"

namespace molstat {
namespace transport {

/**
 * \brief Simulator submodel (molstat::transport::Channel) for transport
 *    through a rectangular barrier. This submodel is designed to represent
 *    ``vacuum'' or background (direct electrode-electrode) tunneling.
 *
 * Inherited model parameters (from molstat::transport::TransportJunction) are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias.
 *
 * Submodel parameters are
 * - `height` (\f$\h\f$), the energy height of the barrier (in eV),
 * - `width` (\f$w\f$), the width of the barrier (in nm).
 *
 * Note that, if this submodel is used, all energies in other submodels (and
 * the Fermi energy) should be reported in eV, as well.
 *
 * The rectangular barrier model is discussed in most standard quantum
 * mechanics texts. As a brief overview, the potential energy is 0 for
 * \f$x<0\f$ and for \f$x>w\f$ (\f$w>0\f$). In between (\f$0<x<w\f$), the
 * potential is \f$h\f$. The Schr\"odinger equation is then solved in each of 
 * the three segments, and both the wavefunction and its derivative are matched
 * at the boundaries between the three regions. Finally, we assume there is no
 * incoming wave from the right and use the magnitude of the wavefunction in
 * the right as the transmission amplitude. This, in square modulus, is the
 * transmission for Landauer-B\"uttiker theory. The equation is
 * \f[
 * T(E) = \frac{1}{ 1 + \frac{1}{\sinh^2( \sqrt{2m(h-E)} w / \hbar ) h^2} {4E(h-E)} }.
 *
 * This model does not take the applied bias directly into account, and should
 * only be used for zero-bias simulations.
 * \f]
 */
class RectangularBarrier : public Channel,
	public DifferentialConductance
{
public:
	/// Container index for the Fermi energy.
	static const std::size_t Index_EF;

	/// Container index for the applied bias.
	static const std::size_t Index_V;

	/// Container index for the height (energy) of the barrier.
	static const std::size_t Index_h;

	/// Container index for the width (distance) of the barrier.
	static const std::size_t Index_w;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	virtual ~RectangularBarrier() = default;

	/**
	 * \brief Calculates the transmission for a set of model parameters.
	 *
	 * \param[in] e The energy of the incident electron.
	 * \param[in] V The applied bias.
	 * \param[in] h The height of the barrier (energy, in eV).
	 * \param[in] w The width of the barrier (distance, in nm).
	 * \return The transmission for this set of parameters.
	 */
	static double transmission(const double e, const double V, const double h,
		const double w);
	
	virtual double DiffG(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

#endif
