/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file echem/simulator_models/transport_observables.h
 * \brief Interfaces for the various observables related to electrochemistry.
 *
 * \author Bo Fu, Matthew G.\ Reuter
 * \date November 2014
 */

#ifndef __electrochemistry_observables_h__
#define __electrochemistry_observables_h__

#include <valarray>
#include <general/simulator_tools/observable.h>

namespace molstat {
namespace echem {

/**
 * \brief Observable class for the potential where electron transfer occurs in
 *    the forward sweep of a cyclic voltammogram.
 */
class ForwardETPotential : public Observable<ForwardETPotential>
{
public:
    ForwardETPotential()
    	: Observable<ForwardETPotential>(&ForwardETPotential::ForwardETP)
    {}

    virtual ~ForwardETPotential() = default;

    /**
     * \brief Returns the electron transfer potential for the forward sweep.
     *
     * \throw molstat::NoObservableProduced if the specified model parameters
     *    do not lead to such a potential.
     *
     * \param[in] params A set of model parameters.
     * \return The forward electron transfer potential.
     */
    virtual double ForwardETP(const std::valarray<double> &params) const = 0;
};

/**
 * \brief Observable class for the potential where electron transfer occurs in
 *    the backward sweep of a cyclic voltammogram.
 */
class BackwardETPotential : public Observable<BackwardETPotential>
{
public:
    BackwardETPotential()
    	: Observable<BackwardETPotential>(&BackwardETPotential::BackwardETP)
    {}

    virtual ~BackwardETPotential() = default;

    /**
     * \brief Returns the electron transfer potential for the backward sweep.
     *
     * \throw molstat::NoObservableProduced if the specified model parameters
     *    do not lead to such a potential.
     *
     * \param[in] params A set of model parameters.
     * \return The backward electron transfer potential.
     */
    virtual double BackwardETP(const std::valarray<double> &params) const = 0;
};

/**
 * \brief Observable class for the redox potential where electron transfer occurs in
 *    the Nernstian reaction.
 */
class RedoxETPotential : public Observable<RedoxETPotential>
{
public:
    RedoxETPotential()
    	: Observable<RedoxETPotential>(&RedoxETPotential::RedoxETP)
    {}

    virtual ~RedoxETPotential() = default;

    /**
     * \brief Returns the redox potential for Nernstian reaction.
     *
     * \param[in] params A set of model parameters.
     * \return The redox potentiali in a Nernstian reaction.
     */
    virtual double RedoxETP(const std::valarray<double> &params) const = 0;
};



} // namespace molstat::echem
} // namespace molstat

#endif
