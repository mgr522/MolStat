/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file bin_log.h
 * \brief Implements logarithmic binning.
 *
 * \author Matthew G.\ Reuter
 * \date June 2014
 */

#ifndef __bin_log_h__
#define __bin_log_h__

#include "bin_style.h"

namespace molstat
{

/**
 * \brief Logarithmic binning style.
 *
 * \f$u = f(x) = \log_b(x)\f$. Then, \f$x = f^{-1}(u) = b^u\f$ and
 * \f$f'(x) = [x \ln(b)]^{-1} \f$.
 */
class BinLog : public BinStyle
{
protected:
	/// The base of the logarithm.
	const double b;

public:
	BinLog() = delete;
	virtual ~BinLog() = default;

	/**
	 * \brief Constructor specifying the base of logarithm to use.
	 *
	 * \param[in] nbin_ The number of bins.
	 * \param[in] b_ The base.
	 */
	BinLog(const std::size_t nbin_, const double b_);

	/**
	 * \brief The mask function, \f$u = f(x) = \log_b(x)\f$.
	 *
	 * \param[in] x The unmasked data value.
	 * \return The transformed (masked) data value.
	 */
	virtual double mask(const double x) const override;

	/**
	 * \brief The inverse mask function, \f$x=f^{-1}(u) = b^u\f$.
	 *
	 * \param[in] u The transformed (masked) data value.
	 * \return The unmasked data value.
	 */
	virtual double invmask(const double u) const override;

	/**
	 * \brief The derivative \f$\mathrm{d}f / \mathrm{d}x = [x \ln(b)]^{-1}\f$.
	 *
	 * \param[in] x The unmasked data value, \f$x\f$.
	 * \return The derivative evaluated at \f$x\f$, \f$f'(x)\f$.
	 */
	virtual double dmaskdx(const double x) const override;

	virtual std::string info() const override;
};

} // namespace molstat

#endif
