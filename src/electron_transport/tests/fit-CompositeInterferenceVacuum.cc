/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file tests/fit-CompositeInterferenceVacuum.cc
 * \brief Test suite for the composite model (interference and vacuum
 *    background) fit model.
 *
 * \test Test suite for the composite model (interference and vacuum
 *    background) fit model.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#include <cassert>
#include <array>
#include <list>
#include <utility>
#include <cmath>
#include <iostream>

#include <electron_transport/fitter_models/composite_interference_vacuum.h>

using namespace std;

/// Shortcut for the model type in this test.
using ModelType = molstat::transport::CompositeInterferenceVacuumFitModel;

/// Alias for the type of data values.
using DataType = std::array<double, 1>;

/**
 * \brief Main function for testing the symmetric-coupling, one-site model.
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
	fitparam[ModelType::F] = 1.;
	fitparam[ModelType::GMINUS] = 1.e-4;
	fitparam[ModelType::NORM] = 1.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(8.73755 - returnval.first) / 8.73755 < thresh);
	assert(abs(-1.07601e-2 - returnval.second[ModelType::F]) / 1.07601e-2 < thresh);
	assert(abs(-22919.8 - returnval.second[ModelType::GMINUS]) / 22919.8 < thresh);
	assert(abs(9.73755 - returnval.second[ModelType::NORM]) / 9.73755 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(44.0665 - returnval.first) / 44.0665 < thresh);
	assert(abs(-2.42792e-2 - returnval.second[ModelType::F]) / 2.42792e-2 < thresh);
	assert(abs(-126428. - returnval.second[ModelType::GMINUS]) / 126428. < thresh);
	assert(abs(45.0665 - returnval.second[ModelType::NORM]) / 45.0665 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(63.5689 - returnval.first) / 63.5689 < thresh);
	assert(abs(-1.22860e-2 - returnval.second[ModelType::F]) / 1.22860e-2 < thresh);
	assert(abs(-249957. - returnval.second[ModelType::GMINUS]) / 249957. < thresh);
	assert(abs(64.5689 - returnval.second[ModelType::NORM]) / 64.5689 < thresh);

	fitparam[ModelType::F] = 0.05;
	fitparam[ModelType::GMINUS] = 9.e-6;
	fitparam[ModelType::NORM] = 0.1;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(0.517917 - returnval.first) / 0.517917 < thresh);
	assert(abs(-1.07171e-4 - returnval.second[ModelType::F]) / 1.07171e-4 < thresh);
	assert(abs(-24901.2 - returnval.second[ModelType::GMINUS]) / 24901.2 < thresh);
	assert(abs(15.1792 - returnval.second[ModelType::NORM]) / 15.1792 < thresh);

	return 0;
}
