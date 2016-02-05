/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file junction.cc
 * \brief Defines the junction interface for simulating electron trasport.
 *
 * A transport junction has several intrinsic parameters, such as the
 * the Fermi energy and the applied bias. Each junction then also some
 * discrete number of channels, which ultimately lead to conductance.
 *
 * The junction is thus modeled as a \"composite\" model, where each submodel
 * is a channel. This file defines these concepts.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "junction.h"
#include "rectangular_barrier.h"

namespace molstat {
namespace transport {

const std::size_t TransportJunction::Index_EF = 0;
const std::size_t TransportJunction::Index_V = 1;

std::vector<std::string> TransportJunction::get_names() const
{
	std::vector<std::string> ret(2);

	ret[Index_EF] = "ef";
	ret[Index_V] = "v";

	return ret;
}

/// \cond
TransportJunction::TransportJunction() :
	CompositeObservable<ElectricCurrent>(
		std::plus<double>()
	),
	CompositeObservable<StaticConductance>(
		std::plus<double>()
	),
	CompositeObservable<ZeroBiasConductance>(
		std::plus<double>()
	),
	CompositeObservable<DifferentialConductance>(
		std::plus<double>()
	)
{
}
/// \endcond

double TransportJunction::AppBias(const std::valarray<double> &params) const
{
	return params[Index_V];
}

double TransportJunction::SeebeckS(const std::valarray<double> &params) const
{
	const ObservableIndex zbg { GetObservableIndex<ZeroBiasConductance>() };
	const ObservableIndex s { GetObservableIndex<SeebeckCoefficient>() };

	// parse the parameters into submodels
	const SubmodelParameters subparams{ routeSubmodelParameters(params) };

	// the composite Seebeck coefficient is
	// S = sum_j (g_j S_j) / sum_k (g_k),
	// where g_j is the zero-bias conductance and s_j is the Seebeck coefficient
	// for the channel
	double sumgs { 0. }, sumg{ 0. };

	// go through each submodel
	for(const auto submodel : subparams) {
		// submodel.first is a pointer to the submodel
		// submodel.second is the parameters to pass to it

		// getObservableFunction will throw an exception if the submodel is
		// incompatible with the specified observable
		const double gj =
			submodel.first->getObservableFunction(zbg)(submodel.second);
		const double sj =
			submodel.first->getObservableFunction(s)(submodel.second);

		sumg += gj;
		sumgs += gj * sj;
	}

	return sumgs / sumg;
}

double TransportJunction::DispW(const std::valarray<double> &params) const
{
	// the displacement observable requires one channel to be a
	// rectangular barrier. search through the submodels for it.
	// if not found, throw an exception.

	for(const auto model : submodels)
	{
		const std::shared_ptr<RectangularBarrier> rbp =
			std::dynamic_pointer_cast<RectangularBarrier>(model.first);
		if(rbp != nullptr)
			return rbp->DispW(params);
	}

	throw IncompatibleObservable(
		"The displacement observable requires a rectangular barrier channel.");
}

} // namespace molstat::transport
} // namespace molstat
