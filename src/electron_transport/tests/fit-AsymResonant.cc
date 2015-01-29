/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/fit-AsymResonant.cc
 * \brief Test suite for the asymmetric-coupling, resonant-tunneling fit model.
 *
 * \test Test suite for the asymmetric-coupling, resonant-tunneling fit model.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include <cassert>
#include <array>
#include <list>
#include <utility>
#include <cmath>

#include <electron_transport/fitter_models/asymmetric_resonant.h>

using namespace std;

/// Shortcut for the model type in this test.
using ModelType = molstat::transport::AsymmetricResonantFitModel;

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
	constexpr double data1val{ 0.8 }, data2val{ 0.9 }, data3val{ 0.96 };
	constexpr double data1obs{ 0.1 }, data2obs{ 0.1 }, data3obs{ 0.3 };

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
	fitparam[ModelType::GAMMAL] = 10.;
	fitparam[ModelType::GAMMAR] = 10.;
	fitparam[ModelType::R] = 1.;
	fitparam[ModelType::NORM] = 1.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(-9.96992e-2 - returnval.first) / 9.96992e-2 < thresh);
	assert(abs(-3.62109e-4 - returnval.second[ModelType::GAMMAL]) / 3.62109e-4 < thresh);
	assert(abs(-3.62109e-4 - returnval.second[ModelType::GAMMAR]) / 3.62109e-4 < thresh);
	assert(abs(-6.55443e-3 - returnval.second[ModelType::R]) / 6.55443e-3 < thresh);
	assert(abs(3.00764e-4 - returnval.second[ModelType::NORM]) / 3.00764e-4 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(0.208156 - returnval.first) / 0.208156 < thresh);
	assert(abs(-0.158045 - returnval.second[ModelType::GAMMAL]) / 0.158045 < thresh);
	assert(abs(-0.158045 - returnval.second[ModelType::GAMMAR]) / 0.158045 < thresh);
	assert(abs(-2.87374 - returnval.second[ModelType::R]) / 2.87374 < thresh);
	assert(abs(0.308156 - returnval.second[ModelType::NORM]) / 0.308156 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(13.7398 - returnval.first) / 13.7398 < thresh);
	assert(abs(-2.24408 - returnval.second[ModelType::GAMMAL]) / 2.24408 < thresh);
	assert(abs(-2.24408 - returnval.second[ModelType::GAMMAR]) / 2.24408 < thresh);
	assert(abs(-41.5747 - returnval.second[ModelType::R]) / 41.5747 < thresh);
	assert(abs(14.0398 - returnval.second[ModelType::NORM]) / 14.0398 < thresh);


	fitparam[ModelType::GAMMAL] = 5.;
	fitparam[ModelType::GAMMAR] = 9.;
	fitparam[ModelType::R] = 1.4;
	fitparam[ModelType::NORM] = 3.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(10.4476 - returnval.first) / 10.4476 < thresh);
	assert(abs(-13.0172 - returnval.second[ModelType::GAMMAL]) / 13.0172 < thresh);
	assert(abs(5.80338 - returnval.second[ModelType::GAMMAR]) / 5.80338 < thresh);
	assert(abs(-10.4825 - returnval.second[ModelType::R]) / 10.4825 < thresh);
	assert(abs(3.51585 - returnval.second[ModelType::NORM]) / 3.51585 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(41.3249 - returnval.first) / 41.3249 < thresh);
	assert(abs(-2.96429 - returnval.second[ModelType::GAMMAL]) / 2.96429 < thresh);
	assert(abs(5.48431 - returnval.second[ModelType::GAMMAR]) / 5.48431 < thresh);
	assert(abs(-31.9057 - returnval.second[ModelType::R]) / 31.9057 < thresh);
	assert(abs(13.8083 - returnval.second[ModelType::NORM]) / 13.8083 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(40.9437 - returnval.first) / 40.9437 < thresh);
	assert(abs(37.4447 - returnval.second[ModelType::GAMMAL]) / 37.4447 < thresh);
	assert(abs(-19.7985 - returnval.second[ModelType::GAMMAR]) / 19.7985 < thresh);
	assert(abs(-19.6519 - returnval.second[ModelType::R]) / 19.6159 < thresh);
	assert(abs(13.7479 - returnval.second[ModelType::NORM]) / 13.7479 < thresh);

	return 0;
}
