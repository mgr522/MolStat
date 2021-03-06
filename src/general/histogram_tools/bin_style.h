/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file bin_style.h
 * \brief Provides an interface for various binning styles.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __bin_style_h__
#define __bin_style_h__

#include <string>
#include <memory>
#include <vector>
#include <general/string_tools.h>

namespace molstat
{

/**
 * \brief Interface for an arbitrary binning style.
 *
 * The purpose of histograms (at least for determining single-molecule
 * conductances) is to estimate the probability density function of the
 * conductance. In the simplest case, we bin on a \"linear\" scale; that is,
 * bin in \f$x\f$ to directly estimate \f$P_{\hat{x}}(x)\f$.
 *
 * But, this is not the only way to perform the binning, and other venues may
 * produce better results. This class interfaces all of the mathematical
 * quantities needed for a binning style:
 * -# The mask function for binning in \f$u\f$: \f$ u = f(x) \f$.
 * -# The inverse of the mask function: \f$ x = f^{-1}(u) \f$.
 * -# The derivative of the mask function: \f$ \mathrm{d}f / \mathrm{d}x \f$.
 * .
 * Full details about binning styles can be found in \ref sec_histograms.
 */
class BinStyle
{
public:
	/// The number of bins.
	const std::size_t nbins;

	BinStyle() = delete;
	virtual ~BinStyle() = default;

	/**
	 * \brief Constructor requires the number of bins.
	 *
	 * \param[in] nbins_ The number of bins.
	 */
	BinStyle(const std::size_t nbins_);

	/**
	 * \brief The mask function, \f$u=f(x)\f$.
	 *
	 * \param[in] x The unmasked data value.
	 * \return The transformed (masked) data value.
	 */
	virtual double mask(const double x) const = 0;

	/**
	 * \brief The inverse mask function, \f$x=f^{-1}(u)\f$.
	 *
	 * \param[in] u The transformed (masked) data value.
	 * \return The unmasked data value.
	 */
	virtual double invmask(const double u) const = 0;

	/**
	 * \brief The derivative \f$\mathrm{d}f / \mathrm{d}x \f$.
	 *
	 * \param[in] x The unmasked data value, \f$x\f$.
	 * \return The derivative evaluated at \f$x\f$, \f$f'(x)\f$.
	 */
	virtual double dmaskdx(const double x) const = 0;

	/**
	 * \brief Create a string summary of this binning style.
	 *
	 * \return The string representation.
	 */
	virtual std::string info() const = 0;
};

/**
 * \brief Gets a binning style from a vector of string tokens.
 *
 * The first token is the name of the binning style. All subsequent tokens, if
 * any, are options for that binning style.
 *
 * This function destroys the tokens.
 *
 * \throw invalid_argument if a BinStyle object cannot be constructed from the
 *    provided arguments.
 *
 * \param[in] tokens The vector of input tokens for creating a BinStyle.
 * \return The BinStyle object.
 */
std::unique_ptr<BinStyle> BinStyleFactory(TokenContainer &&tokens);

} // namespace molstat

#endif
