/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file histogram.cc
 * \brief Provides tools for constructing histograms.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "histogram.h"
#include "bin_style.h"
#include <limits>

namespace molstat {

Histogram::Histogram(std::size_t ndim_)
	: haveBinned(false), ndim(ndim_), data(),
	  extremes(ndim_, {{std::numeric_limits<double>::max(),
	                    std::numeric_limits<double>::lowest()}}),
	  nbin_dim(ndim_, 0), bin_value(0), binned_data(0)
{
}

void Histogram::add_data(std::valarray<double> v)
{
	if(haveBinned)
		throw std::runtime_error("Cannot add data after binning the histogram.");

	if(v.size() != ndim)
		throw std::invalid_argument("Data has incorrect dimensionality.");

	// check the limits
	for(std::size_t j = 0; j < ndim; ++j)
	{
		if(v[j] < extremes[j][0]) // are we smaller than the min?
			extremes[j][0] = v[j];
		if(v[j] > extremes[j][1]) // are we larger than the max?
			extremes[j][1] = v[j];
	}

	// move the data
	data.emplace_front(move(v));
}

void Histogram::bin_data(
	const std::vector<std::shared_ptr<const BinStyle>> &binstyles)
{
	if(haveBinned)
		throw std::runtime_error("Data has already been binned.");

	if(binstyles.size() != ndim)
		throw std::invalid_argument("Incorrect number of binning styles.");

	// make sure that, if more than 1 bin is specified in a dimension, there is
	// a range of data values
	std::size_t total_bins = 1;
	for(std::size_t j = 0; j < ndim; ++j)
	{
		if(binstyles[j] == nullptr)
			throw std::runtime_error(
				"No binning style specified for dimension " + std::to_string(j) +
					".");

		if(binstyles[j]->nbins == 0)
			throw std::runtime_error(
				"There must be at least 1 bin in every dimension.");

		if(extremes[j][0] == extremes[j][1] && binstyles[j]->nbins != 1)
			throw std::runtime_error("Unable to bin data with >1 bins in " \
				"dimension" + std::to_string(j) + ".");

		total_bins *= binstyles[j]->nbins;
	}

	// determine the bounds of each dimension (in masked coordinates), as well
	// as the width of each bin (in masked coordinates)
	std::vector<std::array<double, 3>> bounds(ndim);
	bin_value.resize(ndim);
	for(std::size_t j = 0; j < ndim; ++j)
	{
		bounds[j][0] = binstyles[j]->mask(extremes[j][0]); // minimum
		bounds[j][1] = binstyles[j]->mask(extremes[j][1]); // maximum
		bounds[j][2] = (bounds[j][1] - bounds[j][0]) / binstyles[j]->nbins;
	
		// set the values of the bins for this dimension
		bin_value[j] = bin_values(bounds[j][0], bounds[j][1], bounds[j][2],
			binstyles[j]);
	}

	// set the size of binned_data vector
	binned_data.resize(total_bins);
	for(std::size_t j = 0; j < total_bins; ++j)
		binned_data[j] = 0.;

	// set up the indexor for accesing the bins in binned_data
	// -> need to get the number of bins in each dimension
	for(std::size_t j = 0; j < ndim; ++j)
		nbin_dim[j] = binstyles[j]->nbins;
	CounterIndex ci{ nbin_dim };

	// go through the data
	while(!data.empty())
	{
		for(std::size_t j = 0; j < ndim; ++j)
		{
			// convert each value to the masked space
			double element = binstyles[j]->mask(data.front()[j]);

			// figure out which bin for this dimension
			if(binstyles[j]->nbins == 1)
				ci.setIndex(j, 0);
			else
				ci.setIndex(j, static_cast<std::size_t>(
					(element - bounds[j][0]) / bounds[j][2]
				));
		}

		// increase the bin count
		binned_data[ci.arrayOffset()] += 1.;

		// discard this data point
		data.pop_front();
	}

	// apply the weight function to account for the bin sizes
	for(ci.reset(); !ci.at_end(); ++ci)
	{
		for(std::size_t j = 0; j < ndim; ++j)
		{
			binned_data[ci.arrayOffset()] *= binstyles[j]->dmaskdx(
				bin_value[j][ci[j]]);
		}
	}
}

std::vector<double> Histogram::bin_values(double dmin, double dmax,
	double dwidth, std::shared_ptr<const BinStyle> bstyle)
{
	std::vector<double> ret(bstyle->nbins); // initialize the size of the vector

	// lower bound of bin j is (dmin + j * dwidth)
	// upper bound of bin j is (dmin + (j+1) * dwidth)
	// we unmask these values and get the average
	for(std::size_t j = 0; j < bstyle->nbins; ++j)
	{
		const double lower = bstyle->invmask(dmin + j * dwidth);
		const double upper = bstyle->invmask(dmin + (j+1) * dwidth);
		ret[j] = 0.5 * (lower + upper);
	}

	return ret;
}

CounterIndex Histogram::begin() const
{
	return CounterIndex(nbin_dim);
}

std::valarray<double> Histogram::getCoordinates(const CounterIndex &index)
	const
{
	std::valarray<double> ret(ndim);

	// set the coordinate for each dimension
	for(std::size_t j = 0; j < ndim; ++j)
		ret[j] = bin_value[j][index[j]];

	return ret;
}

double Histogram::getBinCount(const CounterIndex &index) const
{
	return binned_data[index.arrayOffset()];
}

} // namespace molstat
