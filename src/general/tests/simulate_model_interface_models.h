/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file simulate_model_interface_models.h
 * \brief Dummy models for testing the various MolStat classes for
 *    simulating data.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __test_simulate_model_interface_models_h__
#define __test_simulate_model_interface_models_h__

#include <general/simulator_tools/simulate_model.h>
#include "simulate_model_interface_observables.h"

/// Dummy model class.
class BasicTestModel :
	public BasicObs1,
	public BasicObs2,
	public BasicObs3
{
public:
	/// The set value for BasicObs2.
	constexpr static double obs2value = 4.;

	/// Value to throw for BasicObs3.
	constexpr static int obs3except = 4;

	virtual double Obs1(const valarray<double> &params) const override
	{
		return params[0];
	}

	virtual double Obs2(const valarray<double> &params) const override
	{
		return obs2value;
	}

	virtual double Obs3(const valarray<double> &params) const override
	{
		throw(obs3except);
		return 0.;
	}

	virtual vector<string> get_names() const override
	{
		return { "a" };
	}
};

/// Type of submodel used in the test.
class TestSubmodelType
	: public molstat::SimulateSubmodel<TestSubmodelType>
{
};

/**
 * \brief Dummy composite model class.
 *
 * This is intended to be a base class; derived classes must specify the
 * operation used to combine observables from the submodels.
 */
class CompositeTestModel :
	public molstat::UseSubmodelType<TestSubmodelType>,
	public BasicObs4,
	public molstat::CompositeObservable<BasicObs1>
{
protected:
	virtual vector<string> get_names() const override
	{
		return { "ef", "v" };
	}

public:
	/**
	 * \brief Constructor for a dummy composite model that requires a function
	 *    for combining the composite observables.
	 *
	 * \param[in] oper The operation for combining observables.
	 */
	CompositeTestModel(const std::function<double(double, double)> &oper) :
		molstat::CompositeObservable<BasicObs1>(oper) {}

	virtual double Obs4(const valarray<double> &params) const override
	{
		return params[0] + params[1];
	}
};

/// Dummy composite model class using addition to combine submodels.
class CompositeTestModelAdd
	: public CompositeTestModel
{
public:
	/// Constructor that specifies how to combine two submodels.
	CompositeTestModelAdd() :
		CompositeTestModel(
			std::plus<double>()
		) {}
};

/// Dummy composite model class using multiplication to combine submodels.
class CompositeTestModelMultiply
	: public CompositeTestModel
{
public:
	/// Constructor that specifies how to combine two submodels.
	CompositeTestModelMultiply() :
		CompositeTestModel(
			std::multiplies<double>()
		) {}
};

/// Dummy submodel type.
class CompositeSubModel
	: public TestSubmodelType,
	  public BasicObs1
{
protected:
	virtual vector<string> get_names() const override
	{
		return { "eps", "gamma" };
	}

public:
	virtual double Obs1(const valarray<double> &params) const override
	{	
		// v * (eps - gamma)
		return params[1] * (params[2] - params[3]);
	}
};

/**
 * \brief Failed submodel; used to test what happens when a submodel fails to
 *    implement a requested observable.
 */
class FailedSubModel
	: public TestSubmodelType
{
protected:
	virtual vector<string> get_names() const override
	{
		return { "eps" };
	}
};

#endif
