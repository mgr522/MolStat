/**
 * \file histogram1d.h
 * \brief Provides an interface for the GSL histogram functions.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __histogram1d_h__
#define __histogram1d_h__

#include <memory>
#include <gsl/gsl_histogram.h>
#include <queue>

/**
 * \brief Interface for making histograms using GSL.
 *
 * In the beginning, this class simply collects (and stores) data that will
 * eventually be binned into a histogram. At some point, the user decides to
 * bin; the bounds of the bins are determined by the data. After the initial
 * binning, data can still be added. It is not stored, but rather just added
 * to the bin counts.
 *
 * An iterator-style class is also provided for iterating through the bins.
 */
class Histogram1D {
protected:
	/**
	 * \brief List of the data to be binned into the histogram.
	 */
	std::queue<double> data;

	/**
	 * \brief The minimum value of the variable.
	 */
	double minval;

	/**
	 * \brief The maximum value of the variable.
	 */
	double maxval;

	/**
	 * \brief The handle to the GSL histogram functions.
	 */
	std::shared_ptr<gsl_histogram> hist;

public:
	class const_iterator;

	/**
	 * \brief Default constructor.
	 */
	Histogram1D();

	/**
	 * \brief Adds a data element to the histogram.
	 *
	 * \param[in] v The variable.
	 */
	void add_data(const double v);

	/**
	 * \brief Bin the based on what's stored in the queue.
	 *
	 * \param[in] n Number of bins to use.
	 */
	void bin(const size_t n);

	/**
	 * \brief Get an iterator to the front of this histogram.
	 *
	 * \return An iterator to the beginning of this histogram.
	 */
	const_iterator begin() const;

	/**
	 * \brief Get an iterator to the end of this histogram.
	 *
	 * \return An iterator to the end of this histogram.
	 */
	const_iterator end() const;

	/**
	 * \brief Iterator class for accessing the histogram bins.
	 */
	class const_iterator {
	private:
		/**
		 * \brief The GSL histogram handle.
		 */
		const std::shared_ptr<const gsl_histogram> hist;

		/**
		 * \brief The bin index.
		 */
		size_t bin;

		/**
		 * \brief The value of the variable in the middle of this bin.
		 */
		double val_;

		/**
		 * \brief The bin count of this bin.
		 */
		double bincount_;

		/**
		 * \brief Advances the values of bin to the next bin.
		 */
		void next_bin();

		/**
		 * \brief Sets the values of val_ and bincount_ for this bin.
		 */
		void set_output();

	public:
		/**
		 * \brief Default constructor.
		 */
		const_iterator();

		/**
		 * \brief Constucts the iterator; based on the GSL histogram handle.
		 *
		 * \param[in] h The GSL histogram handle.
	 	 */
		const_iterator(const std::shared_ptr<const gsl_histogram> h);

		/**
		 * \brief Prefix forward iteration operator.
		 *
		 * \return State of the iterator before iteration.
		 */
		const_iterator operator++();

		/**
		 * \brief Postfix forward iteration operator.
		 *
		 * \return State of the iterator after iteration.
		 */
		const_iterator operator++(int);

		/**
		 * \brief Get the value of the variable in the middle of this bin.
		 *
		 * \return The value of the variable in the middle of this bin.
		 */
		double variable() const;

		/**
		 * \brief Get the bin count of this bin.
		 *
		 * \return The bin count of this bin.
		 */
		double bin_count() const;

		/**
		 * \brief Determine if these iterators refer to the same state.
		 *
		 * \param[in] rhs The iterator to compare to.
		 * \return True if they have the same state; false otherwise.
		 */
		bool operator== (const const_iterator &rhs) const;

		/**
		 * \brief Determine if these iterators refer to different states.
		 *
		 * \param[in] rhs The iterator to compare to.
		 * \return True if they have different states; false otherwise.
		 */
		bool operator!= (const const_iterator &rhs) const;

		/**
		 * \brief Histogram1D::end() const needs to set the iterator to the end.
	 	 */
		friend const_iterator Histogram1D::end() const;
	};
};

#endif
