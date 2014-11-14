/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file bin_linear.h
 * \brief Implements linear binning.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 */

#ifndef __bin_linear_h__
#define __bin_linear_h__

#include "bin_style.h"

namespace molstat
{

/**
 * \brief Linear binning style. This is the \"identity\" mask.
 *
 * \f$u = f(x) = x\f$. Trivially, \f$x = f^{-1}(u) = u\f$ and \f$f'(x) = 1\f$.
 */
class BinLinear : public BinStyle
{
public:
	using BinStyle::BinStyle;
	virtual ~BinLinear() = default;

	/**
	 * \brief The mask function, \f$u = f(x) = x\f$.
	 *
	 * \param[in] x The unmasked data value.
	 * \return The transformed (masked) data value.
	 */
	virtual double mask(const double x) const override;

	/**
	 * \brief The inverse mask function, \f$x=f^{-1}(u) = u\f$.
	 *
	 * \param[in] u The transformed (masked) data value.
	 * \return The unmasked data value.
	 */
	virtual double invmask(const double u) const override;

	/**
	 * \brief The derivative \f$\mathrm{d}f / \mathrm{d}x = 1\f$.
	 *
	 * \param[in] x The unmasked value, \f$x\f$.
	 * \return The derivative evaluated at \f$x\f$, \f$f'(x)\f$.
	 */
	virtual double dmaskdx(const double x) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
