/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
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

/**
 * \internal
 * \brief Dummy model class.
 * \endinternal
 */
class TestModel : public virtual molstat::SimulateModel,
	public Observable1,
	public Observable2,
	public Observable3 {

public:
	double Obs1(const valarray<double> &vals) const override {
		return vals[0];
	}

	double Obs2(const valarray<double> &vals) const override {
		return constvalue;
	}

	double Obs3(const valarray<double> &vals) const override {
		throw(exceptvalue);
		return 0.;
	}

	vector<string> get_names() const override {
		return { "a" };
	}

	size_t get_num_parameters() const override {
		return 1;
	}
};

#endif
