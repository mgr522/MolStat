/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_interface_classes.h
 * \brief Test classes for testing the molstat::SimulateObservables,
 *    molstat::SimulateModel, and molstat::Observable templates.
 *
 * \test Test classes for testing the molstat::SimulateObservables,
 *    molstat::SimulateModel, and molstat::Observable templates.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __test_simulate_model_interface_classes_h__
#define __test_simulate_model_interface_classes_h__

#include <memory>
#include <cstdio>
#include <map>
#include <string>
#include <functional>
#include <array>
#include <assert.h>

#include <general/simulator_tools/simulate_model_interface.h>

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
 * \brief Dummy class for testing the order of parameters and the
 *    order_from_map function.
 * \endinternal
 */
class FailedMapModel : public molstat::SimulateModel<3> {
public:
	FailedMapModel(
		const map<string, shared_ptr<molstat::RandomDistribution>> &avail)
		: SimulateModel<3>(avail,
		                   order_from_map({{0, "a"}, {1, "b"}, {3, "c"}})) {

	}
};

/**
 * \internal
 * \brief Dummy class for testing the order of parameters and the
 *    order_from_map function.
 * \endinternal
 */
class GoodMapModel : public molstat::SimulateModel<3> {
private:
	/**
	 * \internal
	 * \brief Map of indices to names of distributions.
	 * \endinternal
	 */
	const static map<size_t, string> map_of_order;

public:
	/**
	 * \internal
	 * \brief Order of the names of distributions.
	 * \endinternal
	 */
	const array<string, 3> order;

	GoodMapModel(
		const map<string, shared_ptr<molstat::RandomDistribution>> &avail)
		: SimulateModel<3>(avail,
		                   order_from_map(map_of_order)),
		  order(order_from_map(map_of_order)) {
	}
};

const map<size_t, string> GoodMapModel::map_of_order{{1, "c"}, {0, "a"},
	                                                  {2, "a"}, {3, "b"}};

#endif
