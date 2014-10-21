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
 * Any observable that can be simulated should, in some way derive from this
 * class. The derived class will have a function of signature
 * \verbatim
 * double DerivedClass::function_name(const std::valarray<double> &) const
 * \endverbatim
 * that calculates the observable, given a set of model parameters.
 *
 * The default constructor of molstat::Observable is intentionally deleted; the
 * deriving class must tell molstat::Observable the name of the function.
 * molstat::Observable then registers the observable function with the
 * molstat::SimulateModel, letting the runtime know that the model and
 * observable are \"compatible\".
 *
 * \tparam T The derived class for a specific observable.
 */
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
	Observable(double (T::*obsfunc)(const std::valarray<double> &) const) {
		using namespace std;

		// add the function to the list of compatible observables
		compatible_observables[GetObservableIndex<T>()] =
			[obsfunc] (shared_ptr<const SimulateModel> model)
				-> ObservableFunction {

				// cast the model to the observable (derived) class
				shared_ptr<const T> cast = dynamic_pointer_cast<const T>(model);

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
 * \tparam T The derived class for a specific observable.
 */
template<typename T>
class CompositeObservable : public virtual CompositeSimulateModel {
private:
	/**
	 * \brief Constructs a list of information about each of the submodels.
	 *
	 * The first element in the pair is a list of indices that gets the correct
	 * parameters to send to the submodel. (Remember, the total parameter set
	 * for the composite model contains the composite-specified parameters and
	 * then all the parameters for each submodel.) The second element is the
	 * function for calculating the observable for the submodel.
	 *
	 * \throw molstat::IncompatibleObservable If any of the submodels are
	 *    incompatible with the observable.
	 *
	 * \param[in] cmodel The composite model containing the submodels.
	 * \return The list of observable functions for the observable of type T.
	 */
	static std::list<std::pair<std::valarray<std::size_t>, ObservableFunction>>
		getSubmodelInfo(std::shared_ptr<const CompositeSimulateModel> cmodel);

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
	 * The factory function will also throw molstat::NotCompositeSimulateModel
	 * if the model has not submodels or is not a composite model. If used as
	 * intended (`this` is used; see above), the latter should not happen.
	 *
	 * Note also that CompositeObservable<T> and Observable<T> both refer to
	 * the same observable (class T). A particular MolStat simulator model
	 * should not derive from both Observable<T> and CompositeObservable<T>.
	 *
	 * \param[in] oper Operation used to combine the observables from two
	 *    submodels. This operation should probably be associative and
	 *    commutative.
	 */
	CompositeObservable(const std::function<double(double, double)> &oper) {
		using namespace std;

		// get the index for this observable
		const ObservableIndex oindex{ GetObservableIndex<T>() };

		// add the function to the list of compatible observables
		compatible_observables[oindex] =
			[oper, oindex] (shared_ptr<const SimulateModel> model)
				-> ObservableFunction {

				// cast the model to a CompositeModel class so we can access
				// the submodels
				shared_ptr<const CompositeSimulateModel> cast
					= dynamic_pointer_cast<const CompositeSimulateModel>(model);

				// verify that the model is a composite model and has submodels
				if(cast == nullptr || cast->submodels.size() == 0)
					throw NotCompositeSimulateModel();

				// construct a list of the observable functions for each of the
				// the submodels. if any submodel does not have such a function,
				// it will throw IncompatibleObservable. we will let that
				// exception pass upward
				const list<pair<valarray<size_t>, ObservableFunction>>
					subinfo{ getSubmodelInfo(cast) };

				// make the actual Observable function
				return [oper, subinfo] (const std::valarray<double> &params)
					-> double {

					double ret{ 0. };
					bool isfirst{ true };

					// go through the submodels:
					// calculate the observable of each and combine them using
					// the specified operator
					for(auto modelinfo : subinfo) {
						// modelinfo.first has a valarray the filters out the correct
						// model parameters to send to the submodel.
						// modelinfo.second is the function
						double obs = modelinfo.second(params[modelinfo.first]);

						if(isfirst) {
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
std::list<std::pair<std::valarray<std::size_t>, ObservableFunction>>
	CompositeObservable<T>::getSubmodelInfo(
		std::shared_ptr<const CompositeSimulateModel> cmodel) {

	// get the index for this observable
	const ObservableIndex oindex{ GetObservableIndex<T>() };

	// get the number of model parameters explicitly required by the
	// composite model
	const std::size_t cparams{ cmodel->get_composite_parameters() };
	std::size_t tally { cparams };

	// construct the information (parameters indices and observable function)
	// for each submodel
	list<ObservableFunction> ret;
	for(const auto submodel : cmodel->submodels) {
		const std::size_t submodel_nparam{ submodel->get_num_parameters() };

		// create a list of the parameters indices that should be passed to this
		// submodel
		std::list<std::size_t> indices;
		// first add in the indices required by the composite model
		for(std::size_t j = 0; j < cparams; ++j)
			indices.emplace_back(j);
		// now add in the indices for the specific submodel
		for(std::size_t j = tally; j < tally + submodel_nparam; ++j)
			indices.emplace_back(j);

		// getObservableFunction will throw IncompatibleObservable if this
		// submodel is incompatible with the observable. let this exception pass
		// up to the caller
		ret.emplace_back(std::make_pair(
			list2valarray(indices),
			submodel->getObservableFunction(oindex)));

		// now add the submodel's parameters to the tally for the next offset
		tally += submodel_nparam;
	}

	return ret;
}

template<typename T>
std::valarray<std::size_t> CompositeObservable<T>::list2valarray(
	const std::list<std::size_t> &list) {

	std::valarray<std::size_t> ret(list.size());
	std::size_t j = 0;

	for(auto index : list) {
		ret[j] = index;
		++j;
	}

	return ret;
}

} // namespace molstat

#endif