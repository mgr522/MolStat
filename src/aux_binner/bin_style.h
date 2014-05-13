/**
 * \file bin_style.h
 * \brief Provides an interface for various binning styles.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __bin_style_h__
#define __bin_style_h__

/**
 * \brief Interface for an arbitrary binning style.
 *
 * The purpose of histograms (at least for determining single-molecule
 * conductances) is to estimate the probability density function of the
 * conductance. In the most simplest case, we bin on a \"linear\" scale; that
 * is, bin in \f$g\f$ to directly estimate \f$P(g)\f$.
 *
 * But, this is not the only way to perform the binning, and other venues may
 * produce better results. For example, it is common to bin logarithmically
 * when the measured conductance values span an order of magnitude (or more).
 * In this case we bin in \f$\log(g)\f$, but still try to estimate \f$P(g)\f$
 * rather than \f$P(\log(g))\f$.
 *
 * Suppose we bin in a variable \f$u(g)\f$ that, as written, is a function of
 * \f$g\f$. The obtained histogram estimates \f$P(u(g))\f$, which, using the
 * change-of-variable formula for probability density functions \cite
 * bk:ghahramani-2001, is related to \f$P(g)\f$ by
 * \f[ P(g) = P(u(g)) \frac{\mathrm{d}u}{\mathrm{d}g}. \f]
 *
 * This class provides a general framework for specifying \f$u(g)\f$,
 * \f$u^{-1}(u)\f$, and \f$\mathrm{d}u / \mathrm{d}g\f$.
 */
class BinStyle {
public:
	/**
	 * \brief Default constructor.
	 */
	BinStyle() = default;

	/**
	 * \brief Destructor.
	 */
	virtual ~BinStyle() = default;

	/**
	 * \brief The conductance mask, \f$u(g)\f$.
	 *
	 * \param[in] g The conductance value.
	 * \return The transformed conductance.
	 */
	virtual double gmask(const double g) const = 0;

	/**
	 * \brief The inverse conductance mask, \f$g=u^{-1}(u)\f$.
	 *
	 * \param[in] u The transformed conductance value.
	 * \return The conductance.
	 */
	virtual double invgmask(const double u) const = 0;

	/**
	 * \brief The derivative \f$\mathrm{d}u / \mathrm{d}g \f$.
	 *
	 * \param[in] g The conductance value, \f$g\f$.
	 * \return The derivative evaluated at \f$g\f$, \f$u'(g)\f$.
	 */
	virtual double dudg(const double g) const = 0;
};

#endif
