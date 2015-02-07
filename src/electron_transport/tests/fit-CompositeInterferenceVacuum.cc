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
	assert(abs(77.3615 - returnval.first) / 77.3615 < thresh);
	assert(abs(-0.107573 - returnval.second[ModelType::F]) / 0.107573 < thresh);
	assert(abs(-229198. - returnval.second[ModelType::GMINUS]) / 229198. < thresh);
	assert(abs(97.3615 - returnval.second[ModelType::NORM]) / 97.3615 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(107.648 - returnval.first) / 107.648 < thresh);
	assert(abs(-6.06773e-2 - returnval.second[ModelType::F]) / 6.06773e-2 < thresh);
	assert(abs(-316070. - returnval.second[ModelType::GMINUS]) / 316070. < thresh);
	assert(abs(112.648 - returnval.second[ModelType::NORM]) / 112.648 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(125.11 - returnval.first) / 125.11 < thresh);
	assert(abs(-2.45579e-2 - returnval.second[ModelType::F]) / 2.45579e-2 < thresh);
	assert(abs(-499900. - returnval.second[ModelType::GMINUS]) / 499900. < thresh);
	assert(abs(129.110 - returnval.second[ModelType::NORM]) / 129.110 < thresh);

	fitparam[ModelType::F] = 0.05;
	fitparam[ModelType::GMINUS] = 6.e-6;
	fitparam[ModelType::NORM] = 0.1;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(-3.91390 - returnval.first) / 3.91390 < thresh);
	assert(abs(-1.16207e-3 - returnval.second[ModelType::F]) / 1.16207e-3 < thresh);
	assert(abs(-373237. - returnval.second[ModelType::GMINUS]) / 373237. < thresh);
	assert(abs(160.861 - returnval.second[ModelType::NORM]) / 160.861 < thresh);

	return 0;
}
