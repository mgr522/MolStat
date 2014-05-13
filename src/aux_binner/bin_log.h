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
 * \f$u(g) = \log(g)\f$. Then, \f$u^{-1}(u) = 10^u\f$ and
 * \f$u'(g) = [g \ln(10)]^{-1} \f$.
 */
class BinLog : public BinStyle {
public:
	/**
	 * \brief Default constructor.
	 */
	BinLog() = default;

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
	virtual double gmask(const double g);

	/**
	 * \brief The inverse conductance mask, \f$g=u^{-1}(u) = 10^u\f$.
	 *
	 * \param[in] u The transformed conductance value.
	 * \return The conductance.
	 */
	virtual double invgmask(const double u);

	/**
	 * \brief The derivative \f$\mathrm{d}u / \mathrm{d}g = [g \ln(10)]^{-1}\f$.
	 *
	 * \param[in] g The conductance value, \f$g\f$.
	 * \return The derivative evaluated at \f$g\f$, \f$u'(g)\f$.
	 */
	virtual double dudg(const double g);
};

#endif
