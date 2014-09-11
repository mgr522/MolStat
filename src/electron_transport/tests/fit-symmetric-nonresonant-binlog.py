# This file is a part of MolStat, which is distributed under the Creative
# Commons Attribution-NonCommercial 4.0 International Public License.
# MolStat (c) 2014, Northwestern University.

##
 # @file tests/fit-symmetric-nonresonant-binlog.py
 # @brief Test suite for fitting to the symmetric, nonresonant-tunneling model
 #    with logarithmic binning.
 # 
 # @test Test suite for fitting to the symmetric, nonresonant-tunneling model
 #    with logarithmic binning. This test is really designed to test the use
 #    of logarithmic binning, not the actual fitting model.
 #
 # See the documentation in fit-symmetric-nonresonant.py for how the histogram
 # data was generated. Note: that process produces a histogram in \f$g\f$; this
 # program tests the use of alternate binning styles. That data was then run
 # through
 # @verbatim
 # awk '//{ print log($1)/log(10), $1*log(10)*$2 }'
 # @endverbatim
 # to produce the data file used by this test.
 #
 # @author Matthew G.\ Reuter
 # @date September 2014

import subprocess

# test the log-binning aspect of the fitter program

## @cond

process = subprocess.Popen('../../molstat-fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricNonresonant\n' \
'symmetric-nonresonant-binlog.dat\n' \
'bin log 10.')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 1.599844e+04\n' \
'c=5.9267e+01, d=9.8774e+00, norm=1.0206e+05\n')

## @endcond
