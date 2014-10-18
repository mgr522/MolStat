/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file observable.h
 * \brief Defines the molstat::Observable class, and helper functions.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __observable_h__
#define __observable_h__

#include <memory>
#include <valarray>
#include <typeinfo>
#include <typeindex>
#include "simulate_model.h"

namespace molstat {

template<typename T>
class Observable : public virtual SimulateModel {
public:
	Observable() = delete;
	virtual ~Observable() = default;

	/**
	 * \brief Constructor that processes observable information.
	 *
	 * Adds the type information and a function for calculating the observable
	 * to the model.
	 *
	 * \throw molstat::IncompatibleObservable if the dynamic cast fails. This
	 *    should not happen if the Observable class is used as intended.
	 *
	 * \param[in] obsfunc Member pointer to the observable function.
	 */
	Observable(double (T::*obsfunc)(const std::valarray<double> &) const) {
		// cast this to the inherited class
		std::shared_ptr<const T> cast
			= std::dynamic_pointer_cast<const T>(shared_from_this());

		if(cast == nullptr)
			throw IncompatibleObservable();

		// add the function to the list of compatible observables
		compatible_observables[std::typeindex(typeid(T))] =
			[cast] (const std::valarray<double> &params) -> double {
				return cast->*obsfunc(params);
			};
	}
};

} // namespace MolStat

#endif