/**
 * \file bin_linear.h
 * \brief Implements linear binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __bin_linear_h__
#define __bin_linear_h__

#include "bin_style.h"

/**
 * \brief Linear binning style.
 *
 * \f$u(g) = g\f$. Trivially, \f$u^{-1}(u) = u\f$ and \f$u'(g) = 1\f$.
 */
class BinLinear : public BinStyle {
public:
	/**
	 * \brief Default constructor.
	 */
	BinLinear() = default;

	/**
	 * \brief Destructor.
	 */
	virtual ~BinLinear() = default;

	/**
	 * \brief The conductance mask, \f$u(g) = g\f$.
	 *
	 * \param[in] g The conductance value.
	 * \return The transformed conductance.
	 */
	virtual double gmask(const double g) const;

	/**
	 * \brief The inverse conductance mask, \f$g=u^{-1}(u) = u\f$.
	 *
	 * \param[in] u The transformed conductance value.
	 * \return The conductance.
	 */
	virtual double invgmask(const double u) const;

	/**
	 * \brief The derivative \f$\mathrm{d}u / \mathrm{d}g = 1\f$.
	 *
	 * \param[in] g The conductance value, \f$g\f$.
	 * \return The derivative evaluated at \f$g\f$, \f$u'(g)\f$.
	 */
	virtual double dudg(const double g) const;
};

#endif
