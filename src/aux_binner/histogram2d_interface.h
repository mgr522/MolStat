/**
 * \file histogram2d_interface.h
 * \brief Provides a class for making 2D histograms from data.
 *
 * \author Matthew G.\ Reuter
 * \date May 2014
 */

#ifndef __histogram2d_interface_h__
#define __histogram2d_interface_h__

#include <memory>
#include <gsl/gsl_histogram2d.h>
#include <queue>

/**
 * \brief Class that manages, and provides options for, binning data into 2D
 *    histograms using the GSL interface.
 *
 * In the beginning, this class simply collects (and stores) data that will
 * eventually be binned into a histogram. The data is in the form of a pair,
 * consisting of two ordered variables. At some point, the user decides to
 * bin; the bounds of the bins are determined by the data. After the initial
 * binning, data can still be added. It is not stored, but rather just added
 * to the bin counts.
 *
 * An iterator-style class is also provided for iterating through the bins.
 */
class Histogram2D {
protected:
	/**
	 * \brief List of the data to be binned into the histogram.
	 */
	std::queue<std::pair<double, double>> data;

	/**
	 * \brief Stores the minimum values of the first and second variables.
	 */
	std::pair<double, double> mins;

	/**
	 * \brief Stores the maximum values of the first and second variables.
	 */
	std::pair<double, double> maxs;

	/**
	 * \brief The handle to the GSL histogram functions.
	 */
	std::shared_ptr<gsl_histogram2d> hist;

public:
	class const_iterator;

	/**
	 * \brief Default constructor.
	 */
	Histogram2D();

	/**
	 * \brief Adds a data element to the histogram.
	 *
	 * \param[in] v1 The first variable.
	 * \param[in] v2 The second variable.
	 */
	void add_data(const double v1, const double v2);

	/**
	 * \brief Bin the based on what's stored in the queue.
	 *
	 * \param[in] n1 Number of bins for the first variable.
	 * \param[in] n2 Number of bins for the second variable.
	 */
	void bin(const size_t n1, const size_t n2);

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
		const std::shared_ptr<const gsl_histogram2d> hist;

		/**
		 * \brief The bin index for the first variable.
		 */
		size_t bin1;

		/**
		 * \brief The bin index for the second variable.
		 */
		size_t bin2;

		/**
		 * \brief The value of the first variable in the middle of this bin.
		 */
		double val1_;

		/**
		 * \brief The value of the second variable in the middle of this bin.
		 */
		double val2_;

		/**
		 * \brief The bin count of this bin.
		 */
		double bincount_;

		/**
		 * \brief Advances the values of bin1 and bin2 to the next bin.
		 */
		void next_bin();

		/**
		 * \brief Sets the values of val1_, val2_, and bincount_ for this bin.
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
		const_iterator(const std::shared_ptr<const gsl_histogram2d> h);

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
		 * \brief Get the value of the first variable in the middle of this bin.
		 *
		 * \return The value of the first variable in the middle of this bin.
		 */
		double variable1() const;

		/**
		 * \brief Get the value of the second variable in the middle of this bin.
		 *
		 * \return The value of the second variable in the middle of this bin.
		 */
		double variable2() const;

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

		friend const_iterator Histogram2D::end() const;
	};
};

#endif
