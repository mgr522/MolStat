/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file test_simulate_model_interface.cc
 * \brief Test suite for the molstat::SimulateObservables,
 *    molstat::SimulateModel, and molstat::Observable templates.
 *
 * \test Tests the molstat::SimulateObservables, molstat::SimulateModel,
 *    and molstat::Observable templates.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include <memory>
#include <cstdio>
#include <map>
#include <string>
#include <functional>
#include <array>
#include <assert.h>

#include <general/random_distributions/constant.h>
#include "simulate_model_interface.h"

using namespace std;

const double distvalue = 7.5;
const double constvalue = 4.;
const int exceptvalue = 4;

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
template<size_t MPs>
class Observable1 : public molstat::Observable<MPs> {
public:
	virtual ~Observable1() = default;
	virtual double Obs1(const array<double, MPs> &vals) const = 0;

	virtual double operator()
		(const array<double, MPs> &vals) const override final {

		return Obs1(vals);
	}
};

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
template<size_t MPs>
class Observable2 : public molstat::Observable<MPs> {
public:
	virtual ~Observable2() = default;
	virtual double Obs2(const array<double, MPs> &vals) const = 0;

	virtual double operator()
		(const array<double, MPs> &vals) const override final {

		return Obs2(vals);
	}
};

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
template<size_t MPs>
class Observable3 : public molstat::Observable<MPs> {
public:
	virtual ~Observable3() = default;
	virtual double Obs3(const array<double, MPs> &vals) const = 0;

	virtual double operator()
		(const array<double, MPs> &vals) const override final {

		return Obs3(vals);
	}
};

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
template<size_t MPs>
class Observable4 : public molstat::Observable<MPs> {
public:
	virtual ~Observable4() = default;
	virtual double Obs4(const array<double, MPs> &vals) const = 0;

	virtual double operator()
		(const array<double, MPs> &vals) const override final {

		return Obs4(vals);
	}
};

/**
 * \internal
 * \brief Dummy model class.
 * \endinternal
 */
class TestModel : public molstat::SimulateModel<1>,
	public Observable1<1>, public Observable2<1>, public Observable3<1> {

public:
	TestModel(const map<string, shared_ptr<molstat::RandomDistribution>> &avail)
		: SimulateModel<1>(avail, { {"a"} }) {
	}

	double Obs1(const array<double, 1> &vals) const override {
		return vals[0];
	}

	double Obs2(const array<double, 1> &vals) const override {
		return constvalue;
	}

	double Obs3(const array<double, 1> &vals) const override {
		throw(exceptvalue);
		return 0.;
	}
};

/**
 * \internal
 * \brief Main function for testing the molstat::SimulateObservables,
 *    molstat::SimulateModel, and molstat::Observable templates.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 * \endinternal
 */
int main(int argc, char **argv) {
	molstat::SimulatorFactory<3> factory;
	shared_ptr<molstat::Simulator<3>> sim;
	map<string, shared_ptr<molstat::RandomDistribution>> parameters;
	molstat::gsl_rng_ptr r(nullptr, &gsl_rng_free);

	try {
		factory = molstat::SimulatorFactory<3>
			::makeFactory<TestModel>(parameters);

		// shouldn't be here because we didn't specify the parameter "a"
		assert(false);
	}
	catch(const runtime_error &e) {
		// we should be here!
	}

	parameters["a"] = make_shared<molstat::ConstantDistribution>(distvalue);

	try {
		// this time it should work!
		factory = molstat::SimulatorFactory<3>
			::makeFactory<TestModel>(parameters);
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// try to set an observable for Observable4.
	// This should fail (TestModel doesn't implement Observable4)
	try {
		factory.setObservable<Observable4>(0);
		assert(false);
	}
	catch(const runtime_error &e) {
		// should be here
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	// try to set an observable to a bad index
	try {
		factory.setObservable<Observable1>(3);
		assert(false);
	}
	catch(const out_of_range &e) {
		// should be here
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// now set observables for Observable1 and Observable2.
	try {
		factory.setObservable<Observable1>(0);
		factory.setObservable<Observable2>(2);
	}
	catch(const runtime_error &e) {
		// this should have worked...
		assert(false);
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	// cast to the simulator now that setup is complete
	sim = factory.create();

	// verify the set of observables generated...
	array<double, 3> data = sim->simulate(r);
	assert(abs(data[0] - distvalue) < 1.e-6);
	assert(abs(data[1] - 0.) < 1.e-6);
	assert(abs(data[2] - constvalue) < 1.e-6);

	// create a new simulator that uses Observable3
	try {
		// this time it should work!
		factory = molstat::SimulatorFactory<3>
			::makeFactory<TestModel>(parameters);
	}
	catch(const runtime_error &e) {
		assert(false);
	}

	// now load in Observable3
	try {
		factory.setObservable<Observable1>(0);
		factory.setObservable<Observable3>(1);
		factory.setObservable<Observable2>(2);
	}
	catch(const runtime_error &e) {
		// the set should have worked.
		assert(false);
	}
	catch(const out_of_range &e) {
		assert(false);
	}

	sim = factory.create();

	try {
		// Observable 3 should throw and exception; make sure we catch it
		data = sim->simulate(r);
		assert(false);
	}
	catch(int i) {
		assert(i == exceptvalue);
	}

	return 0;
}
