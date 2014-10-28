/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \internal
 * \file identity_tools.h
 * \brief Interfaces for an observable and model used to test the simulator.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 * \endinternal
 */

#ifndef __identity_tools_h_
#define __identity_tools_h_

#include <valarray>
#include "observable.h"

namespace molstat {

/**
 * \internal
 * \brief Dummy observable for testing the simulator.
 * \endinternal
 */
class IdentityObservable : public Observable<IdentityObservable>
{
public:
	IdentityObservable()
		: Observable<IdentityObservable>(&IdentityObservable::identity) {}

	virtual ~IdentityObservable() = default;

	/**
	 * \internal
	 * \brief Returns the \"identity\", as defined by the model.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The identity.
	 * \endinternal
	 */
	virtual double identity(const std::valarray<double> &params) const = 0;
};

/**
 * \internal
 * \brief Dummy model that implements the dummy identity observable.
 * \endinternal
 */
class IdentityModel :
	public IdentityObservable
{
protected:
	/**
	 * \internal
	 * \brief List of required distributions.
	 *
	 * \return The list of required distributions.
	 * \endinternal
	 */
	virtual std::vector<std::string> get_names() const override;

public:
	/**
	 * \internal
	 * \brief Gets the identity observable.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The identity observable.
	 * \endinternal
	 */
	virtual double identity(const std::valarray<double> &params) const override;
};

} // namespace molstat

 #endif