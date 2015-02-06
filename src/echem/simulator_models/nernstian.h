/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file nernstian.h
 * \brief Simulator model for electron transfer with a Nernstian reaction.
 *
 * \author Bo Fu, Matthew G.\ Reuter
 * \date January 2015
 */

#ifndef __nernstian_reaction_h__
#define __nernstian_reaction_h__

#include "observables.h"

namespace molstat {
namespace echem {

/**
 * \brief Simulator model for single molecule electron transfer using a
 *    Nernstian reacation with \f$n\f$ electrons transferred.
 *
 * Physical parameters include
 * - \f$T\f$, the temperature,
 * - \f$n\f$, the number of electrons transferred in the reaction.
 * .
 * Note that these physical parameters are not directly used by the model
 * (i.e., they do not need to be specified in the input file), but they are
 * used in the following reduced unit system employed by the model and model
 * parameters.
 * - Electric potential is measured in \f$ k_\mathrm{B} T / (ne)\f$.
 *
 * Model parameters are
 * - `Eref` (\f$E_\mathrm{ref}\f$), the reference potential,
 * - `Af` (\f$A_\mathrm{f}\f$), the prefactor of forward half-reaction rate
 *   constant,
 * - `Ab` (\f$A_\mathrm{b}\f$), the prefactor of backward half-reaction rate
 *   constant.
 *
 * The redox potential is determined for the half-reaction
 * \f[ O +ne^-\rightleftharpoons R. \f]
 * Because we are interested in single-molecule properties, the probability
 * that the molecule is in the oxidized (reduced) state is used instead of
 * the concentration used in conventional theories. This model is similar to
 * that used by the non-Nernstian reaction
 * \if fullref 
 * (molstat::echem::NonNernstianReaction)
 * \endif
 * -- this model is simply the limit where the potential changes very slowly.
 *
 * Ultimately, the probability \f$ P_\mathrm{O} \f$ is
 * \f[ P_\mathrm{O}(E) = \frac{k_b(E)}{k_f(E)+k_b(E)}, \f]
 * where \f$E\f$ is the potential. Then the redox potential is then obtained
 * from \f$ P_\mathrm{O}(E) = 1/2 \f$:
 * \f[ E = E_\mathrm{ref} - \frac{k_\mathrm{B} T}{ne} \ln\left(\frac{A_\mathrm{b}}{A_\mathrm{f}}\right). \f]
 * In the reduced unit system this becomes
 * \f[ E = E_\mathrm{ref} - \ln\left( \frac{A_\mathrm{b}}{A_\mathrm{f}} \right). \f]
 *
 * Note that the forward and backward electron transfer potentials are identical
 * in this model.
 */
class NernstianReaction : 
	public ForwardETPotential,
	public BackwardETPotential
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
	virtual ~NernstianReaction() = default;

	virtual double ForwardETP(const std::valarray<double> &params) const
		override;
	virtual double BackwardETP(const std::valarray<double> &params) const
		override;
};

} // namespace molstat::echem
} // namespace molstat

#endif
