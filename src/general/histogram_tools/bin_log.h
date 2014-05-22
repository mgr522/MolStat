/**
 * \file bin_log.h
 * \brief Implements logarithmic binning.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __bin_log_h__
#define __bin_log_h__

#include "bin_style.h"

/**
 * \brief Logarithmic binning style.
 *
 * \f$u(g) = \log_b(g)\f$. Then, \f$u^{-1}(u) = b^u\f$ and
 * \f$u'(g) = [g \ln(b)]^{-1} \f$.
 */
class BinLog : public BinStyle {
protected:
	/**
	 * \brief The base of the logarithm.
	 */
	const double b;

public:
	BinLog() = delete;

	/**
	 * \brief Constructor specifying the base of logarithm to use.
	 *
	 * \param[in] b_ The base.
	 */
	BinLog(const double b_);

	/**
	 * \brief Destructor.
	 */
	virtual ~BinLog() = default;

	/**
	 * \brief The conductance mask, \f$u(g) = \log(g)\f$.
	 *
	 * \param[in] g The conductance value.
	 * \return The transformed conductance.
	 */
	virtual double gmask(const double g) const;

	/**
	 * \brief The inverse conductance mask, \f$g=u^{-1}(u) = 10^u\f$.
	 *
	 * \param[in] u The transformed conductance value.
	 * \return The conductance.
	 */
	virtual double invgmask(const double u) const;

	/**
	 * \brief The derivative \f$\mathrm{d}u / \mathrm{d}g = [g \ln(10)]^{-1}\f$.
	 *
	 * \param[in] g The conductance value, \f$g\f$.
	 * \return The derivative evaluated at \f$g\f$, \f$u'(g)\f$.
	 */
	virtual double dudg(const double g) const;
};

#endif
