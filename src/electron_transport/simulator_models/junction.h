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

/**
 * \brief Type representing the channel concept for transport through
 *    junctions.
 */
class Channel : public virtual SimulateModel
{
protected:
	/**
	 * \brief Identifies this model as a channel.
	 */
	virtual SimulateModelType getModelType() const override {
		return std::type_index{ typeid(Channel) };
	}
};

/**
 * \brief Composite model representing a junction.
 */
class TransportJunction :
	public AppliedBias,
	public CompositeObservable<StaticConductance>,
	public CompositeObservable<DifferentialConductance>
{

public:
	/**
	 * \brief Container index for the Fermi energy.
	 */
	static const std::size_t Index_EF;

	/**
	 * \brief Container index for the applied bias.
	 */
	static const std::size_t Index_V;

protected:
	/**
	 * \brief Junctions require channels for submodels.
	 *
	 * \return The ID for a Channel submodel.
	 */
	virtual SimulateModelType getSubmodelType() const override;

	/**
	 * \brief Junctions inherently require the Fermi energy and the applied
	 *    bias.
	 *
	 * \return A vector containing the names of these parameters.
	 */
	virtual std::vector<std::string> get_names() const override;

public:
	/**
	 * \brief Constructor that tells the MolStat framework to add conductances
	 *    from the channels (submodels).
	 */
	TransportJunction();

	/**
	 * \brief Calculates the applied bias, given a set of model parameters.
	 *
	 * \param[in] params A set of model parameters.
	 * \return The applied bias.
	 */
	virtual double AppBias(const std::valarray<double> &params) const override;
};

} // namespace molstat::transport
} // namespace molstat

 #endif