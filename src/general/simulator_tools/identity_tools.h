/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file identity_tools.h
 * \brief Interfaces for an observable and model used to test the simulator.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __identity_tools_h_
#define __identity_tools_h_

#include <valarray>
#include "observable.h"

namespace molstat {

/// Dummy observable for testing the simulator.
class IdentityObservable : public Observable<IdentityObservable>
{
public:
	IdentityObservable()
		: Observable<IdentityObservable>(&IdentityObservable::identity)
	{
	}

	virtual ~IdentityObservable() = default;

	/**
	 * \brief Returns the \"identity\", as defined by the model.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The identity.
	 */
	virtual double identity(const std::valarray<double> &params) const = 0;
};

/// Dummy model that implements the dummy identity observable.
class IdentityModel : public IdentityObservable
{
protected:
	virtual std::vector<std::string> get_names() const override;

public:
	virtual double identity(const std::valarray<double> &params) const override;
};

} // namespace molstat

 #endif
