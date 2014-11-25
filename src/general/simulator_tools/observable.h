/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
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
 *    molstat::CompositeSimulateModel).
 *
 * The default constructor is intentionally deleted; the deriving class must
 * tell molstat::CompositeObservable the name (via the template) of the
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
private:
	/**
	 * \brief Constructs a std::valarray from a std::list.
	 *
	 * std::valarray::operator[] can accept a std::valarray<std::size_t> to
	 * select desired indices, as will be needed when routing parameters from
	 * the composite model to the submodels. It is easiest to construct a list
	 * of the needed parameter indices and then convert the list to a valarray.
	 *
	 * \param[in] list The list of indices.
	 * \return A valarray containing the same indices, in the same order.
	 */
	inline static std::valarray<std::size_t> list2valarray(
		const std::list<std::size_t> &list);

public:
	CompositeObservable() = delete;
	virtual ~CompositeObservable() = default;

	/**
	 * \brief Constructor that processes composite observable information.
	 *
	 * Adds the type information and a factory for the composite observable to
	 * the model information.
	 *
	 * The factory function needs a molstat::SimulateModel to produce the
	 * actual molstat::ObservableFunction. It should presumably be `*this`,
	 * but shared_from_this() can't be invoked here in the constructor. The
	 * factory defers use of `shared_from_this()` to later, and the function
	 * it returns binds the specified model to the observable function.
	 *
	 * The produced factory function will throw molstat::IncompatibleObservable
	 * if any of the underlying submodels are incompatible with the observable.
	 * The factory function will also throw molstat::NoSubmodels is the
	 * composite model has no submodels. Finally, the factory function will
	 * throw molstat::NotCompositeSimulateModel if the model is not a composite
	 * model. If used as intended (`this` is used; see above), the latter should
	 * not happen.
	 *
	 * Note also that CompositeObservable<T> and Observable<T> both refer to
	 * the same observable (class T). A particular MolStat simulator model
	 * should not derive from both Observable<T> and CompositeObservable<T>.
	 *
	 * \param[in] oper Operation used to combine the observables from two
	 *    submodels. This operation should probably be associative and
	 *    commutative.
	 */
	CompositeObservable(const std::function<double(double, double)> &oper)
	{
		using namespace std;

		// get the index for this observable
		const ObservableIndex oindex{ GetObservableIndex<T>() };

		// add the function to the list of compatible observables
		compatible_observables[oindex] =
			[oper, oindex] (shared_ptr<const SimulateModel> model)
				-> ObservableFunction
			{
				// cast the model to a CompositeModel class so we can access
				// the submodels
				shared_ptr<const CompositeSimulateModel> cmodel
					= dynamic_pointer_cast<const CompositeSimulateModel>(model);

				// verify that the model is a composite model and has submodels
				if(cmodel == nullptr)
					throw NotCompositeSimulateModel();
				if(cmodel->getSubmodels().size() == 0)
					throw NoSubmodels();

				// construct a list of submodel information; that is, the list of
				// parameters to pass to each submodel as well as the observable
				// function for each observable.
				//
				// if any submodel does not have such a function,
				// IncompatibleObservable will be thrown, let the exception pass up
				//
				// ideally, this construction would be a subfunction, but there are
				// issues accessing the protected members of CompositeSimulateModel
				// because we're in a lambda
				list<pair<const valarray<size_t>, ObservableFunction>> subinfo;

				// get the number of model parameters explicitly required by the
				// composite model
				const size_t cparams{ cmodel->get_num_composite_parameters() };
				size_t tally { cparams };

				// go through all of the submodels
				for(const auto submodel : cmodel->getSubmodels())
				{
					const size_t sub_nparam{ submodel->get_num_parameters() };

					// create a list of the parameters indices that should be passed
					// to this submodel
					list<size_t> indices;
					// first add in the indices required by the composite model
					for(size_t j = 0; j < cparams; ++j)
						indices.emplace_back(j);
					// now add in the indices for the specific submodel
					for(size_t j = tally; j < tally + sub_nparam; ++j)
						indices.emplace_back(j);

					// getObservableFunction will throw IncompatibleObservable if
					// this submodel is incompatible with the observable. let this
					// exception pass up to the caller
					subinfo.emplace_back(make_pair(
						list2valarray(indices),
						submodel->getObservableFunction(oindex)));

					// now add the submodel's parameters to the tally for the next
					// offset
					tally += sub_nparam;
				}

				// make the actual Observable function for the composite observable
				return [oper, subinfo] (const std::valarray<double> &params)
						-> double
					{
						double ret{ 0. };
						bool isfirst{ true };

						// go through the submodels:
						// calculate the observable of each and combine them using
						// the specified operator
						for(auto modelinfo : subinfo)
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
					};
			};
	}
};

// template definitions
template<typename T>
std::valarray<std::size_t> CompositeObservable<T>::list2valarray(
	const std::list<std::size_t> &list)
{
	std::valarray<std::size_t> ret(list.size());
	std::size_t j = 0;

	for(auto index : list)
	{
		ret[j] = index;
		++j;
	}

	return ret;
}

} // namespace molstat

#endif