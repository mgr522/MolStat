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
	assert(abs(552.090 - returnval.first) / 552.090 < thresh);
	assert(abs(2.71832 - returnval.second[ModelType::C]) / 2.71832 < thresh);
	assert(abs(-69.8812 - returnval.second[ModelType::D]) / 69.8812 < thresh);
	assert(abs(-375656. - returnval.second[ModelType::GMINUS]) / 375656. < thresh);
	assert(abs(72.0899 - returnval.second[ModelType::NT]) / 72.0899 < thresh);
	assert(abs(500. - returnval.second[ModelType::NV]) / 500. < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(918.022 - returnval.first) / 918.022 < thresh);
	assert(abs(0.835953 - returnval.second[ModelType::C]) / 0.835953 < thresh);
	assert(abs(-28.9171 - returnval.second[ModelType::D]) / 28.9171 < thresh);
	assert(abs(-117380. - returnval.second[ModelType::GMINUS]) / 117380. < thresh);
	assert(abs(13.9307 - returnval.second[ModelType::NT]) / 13.9307 < thresh);
	assert(abs(909.091 - returnval.second[ModelType::NV]) / 909.091 < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(1996.79 - returnval.first) / 1996.79 < thresh);
	assert(abs(4.51766e-2 - returnval.second[ModelType::C]) / 4.51766e-2 < thresh);
	assert(abs(-2.52890 - returnval.second[ModelType::D]) / 2.52890 < thresh);
	assert(abs(-11130.9 - returnval.second[ModelType::GMINUS]) / 11130.9 < thresh);
	assert(abs(0.790382 - returnval.second[ModelType::NT]) / 0.790382 < thresh);
	assert(abs(2.e3 - returnval.second[ModelType::NV]) / 2.e3 < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	fitparam[ModelType::C] = 200.;
	fitparam[ModelType::D] = 4.;
	fitparam[ModelType::GMINUS] = 5.e-5;
	fitparam[ModelType::NT] = 3.;
	fitparam[ModelType::NV] = 0.25;
	fitparam[ModelType::NC] = 3.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(205.364 - returnval.first) / 205.364 < thresh);
	assert(abs(-0.804622 - returnval.second[ModelType::C]) / 0.804622 < thresh);
	assert(abs(13.9061 - returnval.second[ModelType::D]) / 13.9061 < thresh);
	assert(abs(-22.2796 - returnval.second[ModelType::GMINUS]) / 22.2796 < thresh);
	assert(abs(32.4545 - returnval.second[ModelType::NT]) / 32.4545 < thresh);
	assert(abs(500. - returnval.second[ModelType::NV]) / 500. < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	// check a parameter set where g < gminus
	fitparam[ModelType::C] = 150.;
	fitparam[ModelType::D] = 2.6;
	fitparam[ModelType::GMINUS] = 6.3e-5;
	fitparam[ModelType::NT] = 0.89;
	fitparam[ModelType::NV] = 0.18;
	fitparam[ModelType::NC] = 510.;
	returnval = model.resid_j(fitparam, {{6.2e-5}}, 4450.);
	assert(abs(-1035.84 - returnval.first) / 1035.84 < thresh);
	assert(abs(1.71252e-3 - returnval.second[ModelType::C]) / 1.71252e-3 < thresh);
	assert(abs(-2.31626 - returnval.second[ModelType::D]) / 2.31626 < thresh);
	assert(abs(-562045. - returnval.second[ModelType::GMINUS]) / 562045. < thresh);
	assert(abs(1.04651 - returnval.second[ModelType::NT]) / 1.04651 < thresh);
	assert(abs(16129. - returnval.second[ModelType::NV]) / 16129. < thresh);
	assert(abs(1. - returnval.second[ModelType::NC]) < thresh);

	return 0;
}
