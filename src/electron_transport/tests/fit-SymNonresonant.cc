/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2014 Northwestern University. */

/**
 * \file tests/fit-SymNonresonant.cc
 * \brief Test suite for the symmetric-coupling, nonresonant-tunneling fit
 *    model.
 *
 * \test Test suite for the symmetric-coupling, nonresonant-tunneling fit model.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include <cassert>
#include <array>
#include <list>
#include <utility>
#include <cmath>

#include <electron_transport/fitter_models/symmetric_nonresonant.h>

using namespace std;

/// Shortcut for the model type in this test.
using ModelType = molstat::transport::SymmetricNonresonantFitModel;

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
	constexpr double data1val{ 0.002 }, data2val{ 0.001 }, data3val{ 0.0005 };
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
	fitparam[ModelType::CEPSILON] = 100.;
	fitparam[ModelType::CGAMMA] = 5.;
	fitparam[ModelType::NORM] = 1.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(-0.442849 - returnval.first) / 0.442849  < thresh);
	assert(abs(0.458223 - returnval.second[ModelType::CEPSILON]) / 0.458223 < thresh);
	assert(abs(-10.2359 - returnval.second[ModelType::CGAMMA]) / 10.2359 < thresh);
	assert(abs(19.5572 - returnval.second[ModelType::NORM]) / 19.5572 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(0.868936 - returnval.first) / 0.868936  < thresh);
	assert(abs(0.340943 - returnval.second[ModelType::CEPSILON]) / 0.340943 < thresh);
	assert(abs(-10.7762 - returnval.second[ModelType::CGAMMA]) / 10.7762 < thresh);
	assert(abs(5.86894 - returnval.second[ModelType::NORM]) / 5.86894 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(-3.01677 - returnval.first) / 3.01677  < thresh);
	assert(abs(0.0607698 - returnval.second[ModelType::CEPSILON]) / 0.0607698 < thresh);
	assert(abs(-2.71703 - returnval.second[ModelType::CGAMMA]) / 2.71703 < thresh);
	assert(abs(0.983229 - returnval.second[ModelType::NORM]) / 0.983229 < thresh);
	

	fitparam[ModelType::CEPSILON] = 85.;
	fitparam[ModelType::CGAMMA] = 8.;
	fitparam[ModelType::NORM] = 6.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(-19.9797 - returnval.first) / 19.9797  < thresh);
	assert(abs(3.81479e-3 - returnval.second[ModelType::CEPSILON]) / 3.81479e-3 < thresh);
	assert(abs(-8.52159e-2 - returnval.second[ModelType::CGAMMA]) / 8.52159e-2 < thresh);
	assert(abs(3.38571e-3 - returnval.second[ModelType::NORM]) / 3.38571e-3 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(-4.99986 - returnval.first) / 4.99986  < thresh);
	assert(abs(2.39778e-5 - returnval.second[ModelType::CEPSILON]) / 2.39778e-5 < thresh);
	assert(abs(-7.57867e-4 - returnval.second[ModelType::CGAMMA]) / 7.57867e-4 < thresh);
	assert(abs(2.37842e-5 - returnval.second[ModelType::NORM]) / 2.37842e-5 < thresh);

	return 0;
}
