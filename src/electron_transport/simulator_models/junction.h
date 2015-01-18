/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file junction.h
 * \brief Declares the junction interface for simulating electron trasport.
 *
 * A transport junction has several intrinsic parameters, such as the
 * the Fermi energy and the applied bias. Each junction then also some
 * discrete number of channels, which ultimately lead to conductance.
 *
 * The junction is thus modeled as a \"composite\" model; each submodel
 * is a channel. This file defines these concepts.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __transport_junction_h__
#define __transport_junction_h__

#include "observables.h"

namespace molstat {
namespace transport {

/// Type representing a channel for transport through junctions.
class Channel
	: public SimulateSubmodel<Channel>
{
};

/// Composite model representing a junction.
class TransportJunction :
	public AppliedBias,
	public CompositeObservable<ElectricCurrent>,
	public CompositeObservable<StaticConductance>,
	public CompositeObservable<ZeroBiasConductance>,
	public CompositeObservable<DifferentialConductance>
{
public:
	/// Container index for the Fermi energy.
	static const std::size_t Index_EF;

	/// Container index for the applied bias.
	static const std::size_t Index_V;

protected:
	virtual SimulateModelType getSubmodelType() const override;
	virtual std::vector<std::string> get_names() const override;

public:
	/**
	 * \brief Constructor that tells the MolStat framework to add conductances
	 *    from the channels (submodels).
	 */
	TransportJunction();

	virtual double AppBias(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

 #endif