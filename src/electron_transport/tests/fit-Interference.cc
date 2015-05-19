/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2015 Stony Brook University. */

/**
 * \file tests/fit-Interference.cc
 * \brief Test suite for transport near an interference feature.
 *
 * \test Test suite for transport near an interference feature.
 *
 * \author Matthew G.\ Reuter
 * \date February 2015
 */

#include <cassert>
#include <array>
#include <list>
#include <utility>
#include <cmath>

#include <electron_transport/fitter_models/interference.h>

using namespace std;

/// Shortcut for the model type in this test.
using ModelType = molstat::transport::InterferenceFitModel;

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
	fitparam[ModelType::NORM] = 1.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(0.116917 - returnval.first) / 0.116917 < thresh);
	assert(abs(-2.23383e-3 - returnval.second[ModelType::F]) / 2.23383e-3 < thresh);
	assert(abs(1.11692 - returnval.second[ModelType::NORM]) / 1.11692 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(5.02691 - returnval.first) / 5.02691 < thresh);
	assert(abs(-6.62960e-3 - returnval.second[ModelType::F]) / 6.62960e-3 < thresh);
	assert(abs(6.02691 - returnval.second[ModelType::NORM]) / 6.02691 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(10.1775 - returnval.first) / 10.1775 < thresh);
	assert(abs(-5.58877e-3 - returnval.second[ModelType::F]) / 5.58877e-3 < thresh);
	assert(abs(11.1775 - returnval.second[ModelType::NORM]) / 11.1775 < thresh);
	
	fitparam[ModelType::F] = 50.;
	fitparam[ModelType::NORM] = 0.1;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(-0.990823 - returnval.first) / 0.990823 < thresh);
	assert(abs(-9.17738e-4 - returnval.second[ModelType::F]) / 9.17738e-4 < thresh);
	assert(abs(9.17738e-2 - returnval.second[ModelType::NORM]) / 9.17738e-2 < thresh);

	return 0;
}
