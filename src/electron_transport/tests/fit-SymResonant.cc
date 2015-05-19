/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file tests/fit-SymResonant.cc
 * \brief Test suite for the symmetric-coupling, resonant-tunneling fit model.
 *
 * \test Test suite for the symmetric-coupling, resonant-tunneling fit model.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include <cassert>
#include <array>
#include <list>
#include <utility>
#include <cmath>

#include <electron_transport/fitter_models/symmetric_resonant.h>

using namespace std;

/// Shortcut for the model type in this test.
using ModelType = molstat::transport::SymmetricResonantFitModel;

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
	constexpr double thresh{ 1.e-6 };
	constexpr double data1val{ 0.8 }, data2val{ 0.9 }, data3val{ 0.95 };
	constexpr double data1obs{ 0.1 }, data2obs{ 0.1 }, data3obs{ 0.5 };

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
	fitparam[ModelType::GAMMA] = 10.;
	fitparam[ModelType::NORM] = 1.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(-0.999884 - returnval.first) < thresh);
	assert(abs(-2.91145e-4 - returnval.second[ModelType::GAMMA]) < thresh);
	assert(abs(1.16458e-4 - returnval.second[ModelType::NORM]) < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(-0.856818 - returnval.first) < thresh);
	assert(abs(-0.159091 - returnval.second[ModelType::GAMMA]) < thresh);
	assert(abs(0.143182 - returnval.second[ModelType::NORM]) < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(-0.304849 - returnval.first) < thresh);
	assert(abs(-0.365869 - returnval.second[ModelType::GAMMA]) < thresh);
	assert(abs(0.695151 - returnval.second[ModelType::NORM]) < thresh);
	

	fitparam[ModelType::GAMMA] = 5.;
	fitparam[ModelType::NORM] = 5.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(5.86515 - returnval.first) / 5.86515 < thresh);
	assert(abs(-8.58143 - returnval.second[ModelType::GAMMA]) / 8.58143 < thresh);
	assert(abs(1.37303 - returnval.second[ModelType::NORM]) / 1.37303 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(45.1763 - returnval.first) / 45.1763 < thresh);
	assert(abs(-25.6535 - returnval.second[ModelType::GAMMA]) / 25.6535 < thresh);
	assert(abs(9.23527 - returnval.second[ModelType::NORM]) / 9.23527 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(24.0155 - returnval.first) / 24.0155 < thresh);
	assert(abs(-6.58303 - returnval.second[ModelType::GAMMA]) / 6.58303 < thresh);
	assert(abs(5.00310 - returnval.second[ModelType::NORM]) / 5.00310 < thresh);

	return 0;
}
