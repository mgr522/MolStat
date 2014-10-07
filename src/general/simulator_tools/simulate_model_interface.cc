/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface.cc
 * \brief Implements aspects of abstract classes for simulating histograms.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

// no #includes because this file is included by simulate_model_interface.h

namespace molstat {

template<std::size_t MPs>
SimulateModel<MPs>::SimulateModel(
	const std::map<std::string,
	               std::shared_ptr<RandomDistribution>> &avail,
	const std::array<std::string, MPs> &names)
	: dists{ { nullptr } } {

	// go through each of the requested parameters
	for(std::size_t j = 0; j < MPs; ++j) {
		try {
			dists[j] = avail.at(names[j]);
		}
		catch(const std::out_of_range &e) {
			throw std::runtime_error(names[j].c_str());
		}
	}
}

template<std::size_t MPs>
std::array<double, MPs>
	SimulateModel<MPs>::generateParameters(gsl_rng_ptr &r) const {

	std::array<double, MPs> ret;

	for(std::size_t j = 0; j < MPs; ++j)
		ret[j] = dists[j]->sample(r);

	return ret;
}

template<std::size_t OBS, std::size_t MPs>
std::array<double, OBS>
	ModelSimulator<OBS, MPs>::simulate(gsl_rng_ptr &r) const {

	std::array<double, OBS> ret;

	// generate the model parameters
	const std::array<double, MPs> params = model->generateParameters(r);

	// calculate each observable
	for(std::size_t j = 0; j < OBS; ++j)
		ret[j] = observables[j](params);

	return ret;
}

template<std::size_t OBS>
template<typename T>
SimulatorFactory<OBS> SimulatorFactory<OBS>::makeFactory(
	const std::map<std::string,
	               std::shared_ptr<RandomDistribution>> &avail) {

	static_assert(T::numModelParameters <= MAX_MPs,
		"The desired model has more parameters than supported. Increase " \
		"SimulatorFactory::MAX_MPs as noted in the code.");

	using SimulatorType = ModelSimulator<OBS, T::numModelParameters>;

	SimulatorFactory<OBS> ret;

	std::unique_ptr<SimulatorType> modelObs;

	// allocate the SimulateObservables object.
	modelObs.reset(new SimulatorType());

	// set the SimulateModel within the simulator
	modelObs->model = std::make_shared<T>(avail);

	// initialize all observables to the zero function
	for(std::size_t j = 0; j < OBS; ++j)
		modelObs->observables[j] = SimulatorType::ZeroObs;

	ret.model.reset(modelObs.release());

	return ret;
}

template<std::size_t OBS>
template<template<std::size_t> class T>
SimulatorFactory<OBS> &SimulatorFactory<OBS>::setObservable(std::size_t j) {
	if(j >= OBS)
		throw std::out_of_range("Observable index is too high.");

	// we need to perform some dynamic casting, but this isn't allowed on
	// unique_ptrs. So, get the raw pointer and use it (no worries about
	// ownership)
	Simulator<OBS> *ptr = model.get();
	
	// use the Factory's setObservableMPs function to find the correct number of
	// model parameters, cast to that molstat::ModelSimulator class, and add
	// the observable.
	//
	// setObservableMPs *should* return true (meaning the observable was
	// successfully assigned), but we check for complete prudence	
	if(!SimulatorFactoryHelper::ObservableSetter<MAX_MPs>::template
		setObservableMPs<OBS, T>(ptr, j))
		
		throw std::length_error("Should not be here.");

	return *this;
}

template<std::size_t OBS>
std::unique_ptr<Simulator<OBS>> SimulatorFactory<OBS>::create() noexcept {
	return std::move(model);
}

template<std::size_t MPs>
template<std::size_t OBS, template<std::size_t> class T>
bool SimulatorFactoryHelper::ObservableSetter<MPs>::setObservableMPs(
	Simulator<OBS> *ptr, std::size_t j) {

	ModelSimulator<OBS, MPs> *obs =
		dynamic_cast<ModelSimulator<OBS, MPs>*>(ptr);

	if(obs == nullptr) {
		// this is the wrong value of MPs; try one lower
		return ObservableSetter<MPs-1>::template
			setObservableMPs<OBS, T>(ptr, j);
	}
	else {
		std::shared_ptr<T<MPs>> cast =
			std::dynamic_pointer_cast<T<MPs>>(obs->model);

		if(cast == nullptr)
			throw std::runtime_error("Incompatible model and observable.");

		obs->observables[j] =
			[cast] (const std::array<double, MPs> &params)
				-> double {

				return (cast->operator())(params);
			};

		return true;
	}
}

template<std::size_t OBS, template<std::size_t> class T>
constexpr bool SimulatorFactoryHelper::ObservableSetter<0>::setObservableMPs(
	Simulator<OBS> *ptr, std::size_t j) noexcept {

	return false;
}

template<std::size_t OBS, typename T>
SimulateModelFunction<OBS> GetSimulateModelFunction() {
	return [] (const std::map<std::string,
		                       std::shared_ptr<RandomDistribution>> &avail)
			-> SimulatorFactory<OBS> {

		return SimulatorFactory<OBS>::template makeFactory<T>(avail);
	};
}

template<std::size_t OBS, template<std::size_t> class T>
ObservableFunction<OBS> GetObservableFunction() {
	return [] (SimulatorFactory<OBS> &sim, std::size_t j)
			-> void {

		sim.template setObservable<T>(j);
	};
}

} // namespace molstat
