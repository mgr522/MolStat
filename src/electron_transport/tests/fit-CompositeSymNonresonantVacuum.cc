/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file tests/fit-CompositeSymNonresonantVacuum.cc
 * \brief Test suite for the composite model (symmetric-coupling,
 *    nonresonant-tunneling and vacuum background) fit model.
 *
 * \test Test suite for the composite model (symmetric-coupling,
 *    nonresonant-tunneling and vacuum background) fit model.
 *
 * \author Matthew G.\ Reuter
 * \date January 2015
 */

#include <cassert>
#include <array>
#include <list>
#include <utility>
#include <cmath>

#include <electron_transport/fitter_models/composite_symmetric_nonresonant_vacuum.h>

using namespace std;

/// Shortcut for the model type in this test.
using ModelType = molstat::transport::CompositeSymmetricNonresonantVacuumFitModel;

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
	fitparam[ModelType::C] = 100.;
	fitparam[ModelType::D] = 5.;
	fitparam[ModelType::GMINUS] = 1.e-4;
	fitparam[ModelType::NORM] = 1.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(52.0899 - returnval.first) / 52.0899 < thresh);
	assert(abs(2.71832 - returnval.second[ModelType::C]) / 2.71832 < thresh);
	assert(abs(-69.8812 - returnval.second[ModelType::D]) / 69.8812 < thresh);
	assert(abs(-375656. - returnval.second[ModelType::GMINUS]) / 375656. < thresh);
	assert(abs(72.0899 - returnval.second[ModelType::NORM]) / 72.0899 < thresh);

	returnval = model.resid_j(fitparam, data2, data2obs);
	assert(abs(8.93066 - returnval.first) / 8.93066 < thresh);
	assert(abs(0.835953 - returnval.second[ModelType::C]) / 0.835953 < thresh);
	assert(abs(-28.9171 - returnval.second[ModelType::D]) / 28.9171 < thresh);
	assert(abs(-117380. - returnval.second[ModelType::GMINUS]) / 117380. < thresh);
	assert(abs(13.9307 - returnval.second[ModelType::NORM]) / 13.9307 < thresh);

	returnval = model.resid_j(fitparam, data3, data3obs);
	assert(abs(-3.20962 - returnval.first) / 3.20962 < thresh);
	assert(abs(4.51766e-2 - returnval.second[ModelType::C]) / 4.51766e-2 < thresh);
	assert(abs(-2.52890 - returnval.second[ModelType::D]) / 2.52890 < thresh);
	assert(abs(-11130.9 - returnval.second[ModelType::GMINUS]) / 11130.9 < thresh);
	assert(abs(0.790382 - returnval.second[ModelType::NORM]) / 0.790382 < thresh);

	fitparam[ModelType::C] = 200.;
	fitparam[ModelType::D] = 4.;
	fitparam[ModelType::GMINUS] = 5.e-5;
	fitparam[ModelType::NORM] = 2.;
	returnval = model.resid_j(fitparam, data1, data1obs);
	assert(abs(44.9091 - returnval.first) / 44.9091 < thresh);
	assert(abs(-0.536414 - returnval.second[ModelType::C]) / 0.536414 < thresh);
	assert(abs(9.27073 - returnval.second[ModelType::D]) / 9.27073 < thresh);
	assert(abs(-14.8531 - returnval.second[ModelType::GMINUS]) / 14.8531 < thresh);
	assert(abs(32.4545 - returnval.second[ModelType::NORM]) / 32.4545 < thresh);

	// check a parameter set where g < gminus
	fitparam[ModelType::C] = 150.;
	fitparam[ModelType::D] = 2.6;
	fitparam[ModelType::GMINUS] = 6.3e-5;
	fitparam[ModelType::NORM] = 0.89;
	returnval = model.resid_j(fitparam, {{6.2e-5}}, 4450.);
	assert(abs(-4449.07 - returnval.first) / 4449.07 < thresh);
	assert(abs(1.71252e-3 - returnval.second[ModelType::C]) / 1.71252e-3 < thresh);
	assert(abs(-2.31626 - returnval.second[ModelType::D]) / 2.31626 < thresh);
	assert(abs(-562045. - returnval.second[ModelType::GMINUS]) / 562045. < thresh);
	assert(abs(1.04651 - returnval.second[ModelType::NORM]) / 1.04651 < thresh);

	return 0;
}
