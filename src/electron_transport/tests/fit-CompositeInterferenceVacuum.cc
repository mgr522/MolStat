/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
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
	assert(abs(3.86807 - returnval.first) / 3.86807 < thresh);
	assert(abs(-5.37863e-3 - returnval.second[ModelType::F]) / 5.37863e-3 < thresh);
	assert(abs(-11459.9 - returnval.second[ModelType::GMINUS]) / 11459.9 < thresh);
	assert(abs(4.86807 - returnval.second[ModelType::NORM]) / 4.86807 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(21.5295 - returnval.first) / 21.5295 < thresh);
	assert(abs(-1.21355e-2 - returnval.second[ModelType::F]) / 1.21355e-2 < thresh);
	assert(abs(-63213.9 - returnval.second[ModelType::GMINUS]) / 63213.9 < thresh);
	assert(abs(22.5295 - returnval.second[ModelType::NORM]) / 22.5295 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(31.2776 - returnval.first) / 31.2776 < thresh);
	assert(abs(-6.13947e-3 - returnval.second[ModelType::F]) / 6.13947e-3 < thresh);
	assert(abs(-124975. - returnval.second[ModelType::GMINUS]) / 124975. < thresh);
	assert(abs(32.2776 - returnval.second[ModelType::NORM]) / 32.2776 < thresh);

	fitparam[ModelType::F] = 0.05;
	fitparam[ModelType::GMINUS] = 6.e-6;
	fitparam[ModelType::NORM] = 0.1;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(-0.195695 - returnval.first) / 0.195695 < thresh);
	assert(abs(-5.81034e-5 - returnval.second[ModelType::F]) / 5.81034e-5 < thresh);
	assert(abs(-18661.9 - returnval.second[ModelType::GMINUS]) / 18661.9 < thresh);
	assert(abs(8.04305 - returnval.second[ModelType::NORM]) / 8.04305 < thresh);

	return 0;
}
