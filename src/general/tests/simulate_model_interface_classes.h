/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface_classes.h
 * \brief Test classes for testing the various MolStat classes for simulating
 *    data.
 *
 * \test Test classes for testing the various MolStat classes for simulating
 *    data.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __test_simulate_model_interface_classes_h__
#define __test_simulate_model_interface_classes_h__

#include <memory>
#include <valarray>
#include <vector>
#include <string>

#include <general/simulator_tools/observable.h>
#include <general/simulator_tools/simulate_model.h>

using namespace std;

const double distvalue = 7.5;
const double constvalue = 4.;
const int exceptvalue = 4;

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
class Observable1 : public molstat::Observable<Observable1> {
public:
	Observable1() :
		Observable(&Observable1::Obs1) {}

	virtual ~Observable1() = default;
	virtual double Obs1(const valarray<double> &vals) const = 0;
};

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
class Observable2 : public molstat::Observable<Observable2> {
public:
	Observable2() :
		Observable(&Observable2::Obs2) {}

	virtual ~Observable2() = default;
	virtual double Obs2(const valarray<double> &vals) const = 0;
};

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
class Observable3 : public molstat::Observable<Observable3> {
public:
	Observable3() :
		Observable(&Observable3::Obs3) {}

	virtual ~Observable3() = default;
	virtual double Obs3(const valarray<double> &vals) const = 0;
};

/**
 * \internal
 * \brief Dummy observable class.
 * \endinternal
 */
class Observable4 : public molstat::Observable<Observable4> {
public:
	Observable4() :
		Observable(&Observable4::Obs4) {}

	virtual ~Observable4() = default;
	virtual double Obs4(const valarray<double> &vals) const = 0;
};

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
