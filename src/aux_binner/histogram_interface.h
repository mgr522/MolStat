/**
 * \file histogram_interface.h
 * \brief Provides an interface for the GSL histogram (1D and 2D) functions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __histogram_interface_h__
#define __histogram_interface_h__

#include <memory>
#include <array>
#include <functional>

/**
 * \brief Interface for making histograms using GSL.
 *
 * This class provides a consistent interface to the GSL histogram functions
 * (both 1D and 2D). At construction, the user specifies the number of bins in
 * the histogram, as well as the minimum and maximum values for all variables.
 * The data elements are then added (one at a time) to make the histogram.
 *
 * The user may also specify a pair of functions, gmask and invgmask. gmask,
 * the conductance mask, will be applied to data before it is binned. The
 * inverse conductance mask, invgmask, will be applied after binning but before
 * output.
 *
 * An iterator-style class is provided for iterating through the bins.
 */
template <std::size_t N>
class Histogram {
protected:
	/**
	 * \brief Iterator class for accessing the histogram bins.
	 *
	 * This class needs the inverse conductance function for calculating the
	 * middle of the bin.
	 */
	class const_iterator {
	protected:
		/**
		 * \brief The bin index.
		 */
		std::array<size_t, N> bin;

		/**
		 * \brief The value of the variable(s) in the middle of this bin.
		 */
		std::array<double, N> val;

		/**
		 * \brief The bin count of this bin.
		 */
		double bincount;

		/**
		 * \brief The inverse conductance mask.
		 */
		const std::function<double(double)> invgmask;

		/**
		 * \brief Advances the values of bin to the next bin.
		 */
		virtual void next_bin() = 0;

		/**
		 * \brief Sets the output values for this bin.
		 */
		virtual void set_output() = 0;

	public:
		const_iterator() = delete;

		/**
		 * \brief Constructor specifying the inverse conductance mask.
		 */
		const_iterator(const std::function<double(double)> invgmask_);

		/**
		 * \brief Default destructor.
		 */
		virtual ~const_iterator() = default;

		/**
		 * \brief Get the value of the variable in the middle of this bin.
		 *
		 * \return The value of the variable in the middle of this bin.
		 */
		const std::array<double, N> &get_variable() const;

		/**
		 * \brief Get the bin count of this bin.
		 *
		 * \return The bin count of this bin.
		 */
		double get_bin_count() const;
	};

	/**
	 * \brief The function to apply to conductance values before binning.
	 */
	const std::function<double(double)> gmask;

	/**
	 * \brief The inverse mask function for conductances.
	 */
	const std::function<double(double)> invgmask;

public:
	/**
	 * \brief Default constructor. Assumes the identity function for the masks.
	 */
	Histogram();

	/**
	 * \brief Constructor that processes the conductance mask functions.
	 *
	 * \param[in] gmask_ The conductance mask to use when binning data.
	 * \param[in] invgmask_ The inverse conductance mask.
	 */
	Histogram(const std::function<double(double)> gmask_,
		const std::function<double(double)> invgmask_);

	/**
	 * \brief Destructor.
	 */
	virtual ~Histogram() = default;

	/**
	 * \brief Adds a data element to the histogram.
	 *
	 * \param[in] v The variable.
	 */
	virtual void add_data(const std::array<double, N> &v) = 0;
};

// Implementations of non-virtual functions
template<std::size_t N>
Histogram<N>::Histogram()
	: gmask([] (double g) { return g; }),
	  invgmask([] (double g) { return g; }) {
}

template<std::size_t N>
Histogram<N>::Histogram(const std::function<double(double)> gmask_,
	const std::function<double(double)> invgmask_)
	: gmask(gmask_), invgmask(invgmask_) {
}

template<std::size_t N>
Histogram<N>::const_iterator::const_iterator(
	const std::function<double(double)> invgmask_)
	: bincount(0), invgmask(invgmask_) {

	for(std::size_t i = 0; i < N; ++i) {
		bin[i] = 0;
		val[i] = 0.;
	}
}

template<std::size_t N>
const std::array<double, N> &Histogram<N>::const_iterator::get_variable() const
{
	return val;
}

template<std::size_t N>
double Histogram<N>::const_iterator::get_bin_count() const {
	return bincount;
}

#endif
