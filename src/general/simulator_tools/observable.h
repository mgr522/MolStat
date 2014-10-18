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
#include "simulator_exceptions.h"

namespace molstat {

template<typename T>
class Observable : public virtual SimulateModel {
public:
	Observable() = delete;
	virtual ~Observable() = default;

	/**
	 * \brief Constructor that processes observable information.
	 *
	 * Adds the type information and a factory for the observable to the model
	 * information.
	 *
	 * The factory function needs a molstat::SimulateModel to produce the
	 * actual molstat::ObservableFunction. It should presumably be `*this`,
	 * but shared_from_this() can't be invoked in the constructor. The
	 * factory defers use of `shared_from_this()` to later, and the function
	 * it returns binds `*this` to the observable function.
	 *
	 * \param[in] obsfunc Member pointer to the observable function.
	 */
	Observable(double (T::*obsfunc)(const std::valarray<double> &) const) {
		using namespace std;

		// add the function to the list of compatible observables
		compatible_observables[GetObservableIndex<T>()] =
			[obsfunc] (shared_ptr<const SimulateModel> model)
				-> ObservableFunction {

				// cast this to the observable (derived) class
				std::shared_ptr<const T> cast
					= std::dynamic_pointer_cast<const T>(model);

				// make sure the observable is compatible
				if(cast == nullptr)
					throw IncompatibleObservable();

				// make the actual Observable function
				return [cast, obsfunc] (const std::valarray<double> &params)
					-> double {

					return (cast.get()->*obsfunc)(params);
				};
			};
	}
};

} // namespace molstat

#endif