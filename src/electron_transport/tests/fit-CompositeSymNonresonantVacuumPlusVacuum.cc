/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/fit-CompositeSymNonresonantVacuumPlusVacuum.cc
 * \brief Test suite for the composite model (symmetric-coupling,
 *    nonresonant-tunneling and vacuum background) plus vacuum fit model.
 *
 * \test Test suite for the composite model (symmetric-coupling,
 *    nonresonant-tunneling and vacuum background) plus vacuum fit model.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#include <cassert>
#include <array>
#include <list>
#include <utility>
#include <cmath>

#include <electron_transport/fitter_models/composite_symmetric_nonresonant_vacuum_plus_vacuum.h>

using namespace std;

/// Shortcut for the model type in this test.
using ModelType = molstat::transport::CompositeSymmetricNonresonantVacuumPlusVacuumFitModel;

/// Alias for the type of data values.
using DataType = std::array<double, 1>;

/**
 * \brief Main function for testing the specified model.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status: 0 if the code passes the test, non-zero otherwise.
 */
int main(int argc, char **argv)
{
	constexpr double thresh{ 1.e-5 };
	constexpr double data1val{ 0.002 }, data2val{ 0.0011 }, data3val{ 0.0005 };
	constexpr double data1obs{ 20. }, data2obs{ 5. }, data3obs{ 4. };

	// set up a dummy data set
	DataType data1{{ data1val }}, data2{{ data2val }}, data3{{ data3val }};
	list<pair<DataType, double>> datalist;
	datalist.emplace_back(data1, data1obs);
	datalist.emplace_back(data2, data2obs);
	datalist.emplace_back(data3, data3obs);

	// make the model
	ModelType model(datalist);

	// fit parameters
	vector<double> fitparam(model.nfit);

	// return values
	pair<double, vector<double>> returnval;

	// check known values for several parameter sets
	fitparam[ModelType::C] = 100.;
	fitparam[ModelType::D] = 5.;
	fitparam[ModelType::GMINUS] = 1.e-4;
	fitparam[ModelType::NT] = 1.;
	fitparam[ModelType::NV] = 1.;
	fitparam[ModelType::NC] = 0.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(516.033 - returnval.first) / 516.033 < thresh);
	assert(abs(1.35887 - returnval.second[ModelType::C]) / 1.35887 < thresh);
	assert(abs(-34.9341 - returnval.second[ModelType::D]) / 34.9341 < thresh);
	assert(abs(-187829. - returnval.second[ModelType::GMINUS]) / 187829. < thresh);
	assert(abs(36.0327 - returnval.second[ModelType::NT]) / 36.0327 < thresh);
	assert(abs(500. - returnval.second[ModelType::NV]) / 500. < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(911.052 - returnval.first) / 911.052 < thresh);
	assert(abs(0.417719 - returnval.second[ModelType::C]) / 0.417719 < thresh);
	assert(abs(-14.4508 - returnval.second[ModelType::D]) / 14.4508 < thresh);
	assert(abs(-58689.4 - returnval.second[ModelType::GMINUS]) / 58689.4 < thresh);
	assert(abs(6.96077 - returnval.second[ModelType::NT]) / 6.96077 < thresh);
	assert(abs(909.091 - returnval.second[ModelType::NV]) / 909.091 < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(1996.39 - returnval.first) / 1996.39 < thresh);
	assert(abs(2.25534e-2 - returnval.second[ModelType::C]) / 2.25534e-2 < thresh);
	assert(abs(-1.26284 - returnval.second[ModelType::D]) / 1.26284 < thresh);
	assert(abs(-5564.51 - returnval.second[ModelType::GMINUS]) / 5564.51 < thresh);
	assert(abs(0.394624 - returnval.second[ModelType::NT]) / 0.394624 < thresh);
	assert(abs(2.e3 - returnval.second[ModelType::NV]) / 2.e3 < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	fitparam[ModelType::C] = 200.;
	fitparam[ModelType::D] = 4.;
	fitparam[ModelType::GMINUS] = 5.e-5;
	fitparam[ModelType::NT] = 3.;
	fitparam[ModelType::NV] = 0.25;
	fitparam[ModelType::NC] = 3.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(156.682 - returnval.first) / 156.682 < thresh);
	assert(abs(-0.402311 - returnval.second[ModelType::C]) / 0.402311 < thresh);
	assert(abs(6.95304 - returnval.second[ModelType::D]) / 6.95304 < thresh);
	assert(abs(-11.1376 - returnval.second[ModelType::GMINUS]) / 11.1376 < thresh);
	assert(abs(16.2273 - returnval.second[ModelType::NT]) / 16.2273 < thresh);
	assert(abs(500. - returnval.second[ModelType::NV]) / 500. < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	return 0;
}
