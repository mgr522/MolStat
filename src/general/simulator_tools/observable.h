/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file observable.h
 * \brief Defines the molstat::Observable class for observables, the
 *    molstat::CompositeObservable class for observables reliant on submodels,
 *    and other helper functions.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __observable_h__
#define __observable_h__

#include <list>
#include <memory>
#include <utility>
#include <valarray>
#include <typeinfo>
#include <typeindex>
#include <functional>
#include "simulate_model.h"
#include "simulator_exceptions.h"

namespace molstat {

/**
 * \brief Base class for an observable.
 *
 * Any observable that can be simulated should, in some way, derive from this
 * class. The derived class will have a function of signature
 * \code{.cpp}
 * double DerivedClass::function_name(const std::valarray<double> &) const
 * \endcode
 * that calculates the observable, given a set of model parameters.
 *
 * The default constructor of molstat::Observable is intentionally deleted; the
 * deriving class must tell molstat::Observable the name of the observable
 * function. molstat::Observable then registers the observable function with
 * the molstat::SimulateModel, letting the runtime know that the model and
 * observable are \"compatible\".
 *
 * The observable function of the subclass may throw
 * molstat::NoObservableProduced if the specified set of model parameters does
 * not result in the observable being emitted by the system. (Not all models
 * will need this feature.) The main MolStat simulator will report the number
 * of trials that do not result in an observable.
 *
 * \note When dealing with composite simulator models, the
 *    `CompositeObservable` class may be preferential to `Observable`. More
 *    details for `CompositeObservable` are presented in its documentation.
 *
 * \tparam T The derived class for a specific observable. This is designed to
 *    be a \"curiously recurring template pattern\".
 */
template<typename T>
class Observable
	: public virtual SimulateModel
{
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
	 * but shared_from_this() can't be invoked here in the constructor. The
	 * factory defers use of `shared_from_this()` to later, and the function
	 * it returns binds the specified model to the observable function.
	 *
	 * The produced factory function will throw molstat::IncompativeObservable
	 * if the model and observable are incompatible. If used as intended
	 * (`this` is used; see above), this should not happen.
	 *
	 * \param[in] obsfunc Member pointer to the observable function.
	 */
	Observable(double (T::*obsfunc)(const std::valarray<double> &) const)
	{
		using namespace std;

		// add the function to the list of compatible observables
		compatible_observables[GetObservableIndex<T>()] =
			[obsfunc] (shared_ptr<const SimulateModel> model)
				-> ObservableFunction
			{
				// cast the model to the observable (derived) class
				shared_ptr<const T> cast = dynamic_pointer_cast<const T>(model);

				// make sure the observable is compatible
				if(cast == nullptr)
					throw IncompatibleObservable();

				// make the actual Observable function
				return [cast, obsfunc] (const std::valarray<double> &params)
						-> double
					{
						return (cast.get()->*obsfunc)(params);
					};
			};
	}
};

/**
 * \brief Base class for a composite observable; that is, an observable that
 *    is calculated from several submodels (used in conjunction with
 *    \c molstat::CompositeSimulateModel).
 *
 * This class is designed to simplify cases where the observable for the
 * composite model is simply calculated from the observables for each
 * submodel. Example: conductance -- conductance through the entire system
 * is just the sum of each channel's conductance. \c CompositeObservable
 * essentially provides the boilerplate code for such an operation; all the
 * composite model needs to specify is the operation used to combine the
 * observables from two submodels. (Ideally, this operation is associative
 * and commutative).
 *
 * If a composite model uses the observables from the submodels in a more
 * complicated way, it should derive from \c Observable and implement the
 * necessary observable function.
 *
 * The default constructor is intentionally deleted; the deriving class must
 * tell \c molstat::CompositeObservable the name (via the template) of the
 * underlying observable for each submodel. The deriving class must also
 * specify a function that combines two values from submodels into the
 * value of the composite observable.
 *
 * \tparam T The derived class for a specific observable. `T` should be a
 *    derived class of molstat::Observable.
 */
template<typename T>
class CompositeObservable
	: public virtual CompositeSimulateModel
{
protected:
	/**
	 * \brief Generate the \c ObservableFunction for the composite model.
	 *
	 * Create the function that calculates the observable \c T using the
	 * submodels. The observable from each submodel is calculated and all
	 * "sub-observables" are combined using the specified operation.
	 *
	 * \throw molstat::IncompatibleObservable if any of the underlying submodels
	 *    are incompatible with the observable.
	 * \throw molstat::NoSubmodels if the composite model has no submodels.
	 * \throw molstat::NotCompositeSimulateModel if the model is not a composite
	 *    model. If used as intended (`this` is used; see documentation on the
	 *    constructor), this should not happen.
	 *
	 * \param[in] oper Operation used to combine the observables from two
	 *    submodels. This operation should probably be associative and
	 *    commutative.
	 * \param[in] model Calculate the observable using this model.
	 * \return The function that calculates the observable, given a set of
	 *    model parameters.
	 */
	static ObservableFunction getCompositeObservableFunction(
		const std::function<double(double,double)> &oper,
		const std::shared_ptr<const SimulateModel> model)
	{
		// get the index for this observable
		const ObservableIndex oindex{ GetObservableIndex<T>() };

		// cast the model to a CompositeModel class so we can access
		// the submodels
		std::shared_ptr<const CompositeSimulateModel> cmodel
			= std::dynamic_pointer_cast<const CompositeSimulateModel>(model);

		// verify that the model is a composite model and has submodels
		if(cmodel == nullptr)
			throw NotCompositeSimulateModel();
		if(cmodel->submodels.size() == 0)
			throw NoSubmodels();

		// construct a list of submodel information; that is, a list of
		// parameters to pass to each submodel as well as the observable
		// function.
		std::list<std::pair<const std::valarray<size_t>,
		                    ObservableFunction>>
			subinfo;

		// go through all of the submodels
		for(const auto submodel : cmodel->submodels)
		{
			// getObservableFunction will throw IncompatibleObservable if
			// this submodel is incompatible with the observable. let this
			// exception pass up to the caller
			subinfo.emplace_back(make_pair(
				submodel.second,
				submodel.first->getObservableFunction(oindex)));
		}

		// make the actual Observable function for the composite observable.
		// this function goes through each submodel, calculates each
		// "sub-observable", and combines them using the specified operation
		return [oper, subinfo] (const std::valarray<double> &params)
				-> double
			{
				double ret{ 0. };
				bool isfirst{ true };

				// go through the submodels:
				// calculate the observable of each and combine them using
				// the specified operation
				for(const auto modelinfo : subinfo)
				{
					// modelinfo.first has a valarray that filters out the
					// correct model parameters to send to the submodel.
					// modelinfo.second is the function
					double obs = modelinfo.second(params[modelinfo.first]);

					if(isfirst)
					{
						// ret is uninitialized
						ret = obs;
						isfirst = false;
					}
					else
						ret = oper(ret, obs);
				}

				return ret;
			}; // end of the returned ObservableFunction
	}

public:
	CompositeObservable() = delete;
	virtual ~CompositeObservable() = default;

	/**
	 * \brief Constructor that processes composite observable information and,
	 *    if in order, adds the type information and a factory for the composite
	 *    observable to the model information.
	 *
	 * We need a `molstat::SimulateModel` to produce the actual
	 * `molstat::ObservableFunction`. It should presumably be `*this`, but
	 * `shared_from_this()` can't be invoked here in the constructor. We
	 * essentially defer the use of `shared_from_this()` to later by using
	 * `getCompositeObservableFunction`, which requires a `SimulateModel`.
	 *
	 * Note also that `CompositeObservable<T>` and `Observable<T>` both refer to
	 * the same observable (class `T`). A particular MolStat simulator model
	 * should not derive from both `Observable<T>` and `CompositeObservable<T>`.
	 *
	 * \param[in] oper Operation used to combine the observables from two
	 *    submodels. This operation should probably be associative and
	 *    commutative.
	 */
	CompositeObservable(const std::function<double(double, double)> &oper)
	{
		using namespace std::placeholders;

		// get the index for this observable
		const ObservableIndex oindex{ GetObservableIndex<T>() };

		// add the function to the list of compatible observables
		compatible_observables[oindex] = std::bind(
			getCompositeObservableFunction, oper, _1);
	}
};

} // namespace molstat

#endif
