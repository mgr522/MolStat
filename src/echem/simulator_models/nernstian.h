/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file nernstian.h
 * \brief Simulator model for electron transfer with a Nernstian reaction.
 *
 * \author Bo Fu
 * \date January 2015
 */

#ifndef __nernstian_reaction_h__
#define __nernstian_reaction_h__

#include "observables.h"

namespace molstat {
namespace echem {

/**
 * \brief Simulator model for single molecule electrochemistry
 *    under Nernstian limit with a slow enough scanning rate.
 *
 *  Model parameters are
 * - `Eref` (\f$E_{ref}\f$), the standard redox potential,
 * - `Af` (\f$A_f\f$), the prefactor of forward half-reaction rate constant,
 * - `Ab` (\f$A_b\f$), the prefactor of backward half-reaction rate constant,
 * - `T`  (\f$T\f$),   the temperature,
 * - `n`  {\f$n\f$),   the number of electrons involved in the half-reaction.  
 *
 * The probability \f$P_{mathrm{O}}\f$ for Nernstian reaction is
 * \f[ P_\mathrm{O}(E)_{\mathrm{Nernstian}} = \frac{k_b(E)}{k_f(E)+k_b(E)}\f]
 * Then the peak potential can be solved from \f$ P_\mathrm{O}(\epsilon_P)_{\mathrm{Nernstian}}=0.5\f$. 
 * \f[ \epsilon_P=E^0 - \frac{k_B T}{ne}\log\frac{A_b}{A_f}\f]
 * - Peak Potential:
 *   \f[ \epsilon_P=E^0 - \frac{k_B T}{ne}\log\frac{A_b}{A_f}\f]
 */
class NernstianReaction : 
	public RedoxETPotential,
	public virtual molstat::SimulateModel
{
public:
	/// Container index for the reference potential.
	static const std::size_t Index_Eref;

	/// Container index for the prefactor of forward half-reaction rate constant.
	static const std::size_t Index_Af;

	/// Container index for the prefactor of backward half-reaction rate constant.

	static const std::size_t Index_Ab;

protected:
	virtual std::vector<std::string> get_names() const override;

public:
	/// Destructor.
	virtual ~NernstianReaction() = default;

	/**
	 * \brief The redox potential for Nernstian reaction.
	 * \param[in] params The set of model parameters.
	 * \return The redox potential.
	 */
	virtual double RedoxETP(const std::valarray<double> &params) const override;
};

} // namespace molstat::echem
} // namespace molstat

#endif
