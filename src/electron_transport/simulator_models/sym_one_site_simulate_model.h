/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file sym_one_site_simulate_model.h
 * \brief Tight-binding model with one site that couples symmetrically to
 *    both electrodes.
 *
 * \author Matthew G.\ Reuter
 * \date September 2014
 */

#ifndef __sym_one_site_simulate_model_h__
#define __sym_one_site_simulate_model_h__

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <general/random_distributions/rng.h>
#include <general/simulator_tools/simulate_model_interface.h>
#include "transport_observables.h"

using std::shared_ptr;

/**
 * \brief Simulator model for transport through a single site that couples
 *    symmetrically to both electrodes.
 *
 * Model parameters are
 * - `ef` (\f$E_\mathrm{F}\f$), the Fermi energy,
 * - `v` (\f$V\f$), the applied bias,
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
 *   \f[ \frac{2e}{h} \int\limits_{E_\mathrm{F}+(\eta-1)eV}^{E_\mathrm{F}+\eta eV} \mathrm{d}E \frac{\partial}{\partial V} T(E) = \frac{2e^2a}{h} T(E_\mathrm{F}+(\eta-1)eV) - \frac{2e^2a}{h}T(E_\mathrm{F}+\eta eV); \f]
 *   \f[ g(V) = \frac{2e^2}{h} \left[ (\eta-a) T(E_\mathrm{F} + \eta eV) + (a+1-\eta) T(E_\mathrm{F} + (\eta-1)eV) \right]. \f]
 * - Static conductance:
 *   \f[ g(V) = \frac{2e^2}{h} \frac{\Gamma}{eV} \left[ \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (\eta-a) eV}{\Gamma} \right) - \arctan\left( \frac{E_\mathrm{F} - \varepsilon + (\eta-1-a) eV}{\Gamma} \right) \right]. \f]
 */
class SymOneSiteSimulateModel : public SimulateModel,
	public ZeroBiasConductance {

public:
	SymOneSiteSimulateModel() = delete;

	/**
	 * \internal
	 * \brief Default constructor.
	 *
	 * The constructor specifies the required parameters.
	 *
	 * \param[in] avail The available distributions, keyed by name.
	 * \endinternal
	 */
	SymOneSiteSimulateModel(
		const std::map<std::string, shared_ptr<RandomDistribution>> &avail);

	virtual ~SymOneSiteSimulateModel() = default;

	/**
	 * \brief Returns the zero-bias conductance for a randomly-generated set of
	 *    model parameters.
	 * 
	 * \param[in] r The GSL random number generator handle.
	 * \return The zero-bias conductance.
	 */
	virtual std::array<double, 1> ZeroBiasG(shared_ptr<gsl_rng> r) const
		override;
};

#endif
