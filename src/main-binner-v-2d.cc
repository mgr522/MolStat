/*
This work is licensed under the Creative Commons Attribution 3.0 United States
License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/3.0/us/ or send a letter to Creative
Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

Copyright (C) 2013 Oak Ridge National Laboratory
*/
/**
 * \file main-binner-v-2d.cc
 * \brief Main function for binning data into histograms.
 *
 * Reads a list of data from standard in and bins it into a 2D histogram with
 * the specified number of bins. For the conductance histogram application,
 * binning for conductance is done logarithmically. Each line of data should be
 * structured to give the voltage and then the conductance.
 *
 * \todo Allow option to bin linearly (in addition to logarithmically).
 *
 * Two command-line arguments are required:
 *    -# The number of conductance measurements to be read (from standard
 *       input) and binned.
 *    -# The number of bins to use (for both voltage and conductance).
 *
 * The bins are output to standard out. Note that the requested number of bins
 * is presently a maximum. In an effort to reduce noise, bins with zero (or
 * only a few) counts are suppressed. The actual number of bins used is
 * output to standard error.
 *
 * \todo Add options to turn on/off bin suppression.
 *
 * \author Matthew G.\ Reuter
 * \date July 2012, May 2013
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <gsl/gsl_histogram2d.h>

/**
 * \brief Main function for binning.
 *
 * Parses the input parameters and outputs the binned data.
 *
 * \param[in] argc The number of command-line arguments.
 * \param[in] argv The command-line arguments.
 * \return Exit status; 0 for normal.
 */
int main(int argc, char **argv) {
	int nbin, ntrials, i, j, usedbin;
	double t, mint, maxt, *logt, minv, maxv, *v;
	gsl_histogram2d *h;

	// get the command-line arguments
	if(argc != 3) {
		fprintf(stderr, "Usage: ./binner ntrials nbin\n" \
			"   ntrials is the number of trials in the input data\n" \
			"   nbin is the number of bins to use\n" \
			"NOTE: The data is expected through stdin\n");
		return 0;
	}

	ntrials = atoi(argv[1]);
	if(ntrials < 1) {
		fprintf(stderr, "Error: Use at least trial.\n");
		return 0;
	}
	nbin = atoi(argv[2]);
	if(nbin < 1) {
		fprintf(stderr, "Error: Use at least one bin.\n");
		return 0;
	}

	// read in the data
	// get the min and max values, store the logarithms
	logt = (double*)malloc(ntrials*sizeof(double));
	v = (double*)malloc(ntrials*sizeof(double));
	mint = 1.0;
	maxt = 0.0;
	minv = DBL_MAX;
	maxv = -DBL_MAX;
	for(i = 0; i < ntrials; ++i) {
		scanf("%le %le", v + i, &t);

		if(v[i] < minv)
			minv = v[i];
		if(v[i] > maxv)
			maxv = v[i];

		if(t < mint)
			mint = t;
		if(t > maxt)
			maxt = t;

		logt[i] = log10(t);
	}
	mint = log10(mint);
	maxt = log10(maxt);
	//maxt = log10(1.001*maxt); // the upper bound is exclusive in gsl

	// make the histogram (in logarithm space)
	h = gsl_histogram2d_alloc(nbin, nbin);
	gsl_histogram2d_set_ranges_uniform(h, minv, maxv, mint, maxt);
	for(i = 0; i < ntrials; ++i)
		gsl_histogram2d_increment(h, v[i], logt[i]);

	// scale the bins so that the biggest bin value is 1
	gsl_histogram2d_scale(h, 1.0 / gsl_histogram2d_max_val(h));

	// print it out
	// count the number of bins with meaningful population
	usedbin = 0;
	for(i = 0; i < nbin; ++i) {
		gsl_histogram2d_get_xrange(h, i, &minv, &maxv);

		for(j = 0; j < nbin; ++j) {
			t = gsl_histogram2d_get(h, i, j);
#if 0
			// don't print out an empty or sparsely filled bin
			if(t < 0.005*ntrials)
				continue;
#endif
			++usedbin;

			gsl_histogram2d_get_yrange(h, j, &mint, &maxt);
#if 0
			// convert the range back to regular 'g' (from 'log g')
			mint = pow(10.0, mint);
			maxt = pow(10.0, maxt);
#endif
// keep the conductances on log-scale

			// also have to scale the histogram count by the width of the bin and
			// the total number of trials
			printf("%.6e %.6e %.6e\n", 0.5*(maxv + minv), 0.5*(maxt + mint), t);
		}

		// need a newline for plotting with gnuplot; not necessary for
		// mathematica
		//printf("\n");
	}

	// clean up
	gsl_histogram2d_free(h);
	free(logt);
	fprintf(stderr, "%d\n", usedbin);

	return 0;
}
