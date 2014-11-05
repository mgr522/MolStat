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
 # @date November 2014

import subprocess
import math

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
tokens = output[0].split()
assert(tokens[0] == 'Resid')
assert(tokens[1] == '=')
# check the residual -- this is empirical
assert(math.fabs(float(tokens[2]) - 287.) / 287. < 5.e-2)

# check c -- not empirical
cline = tokens[3].split('=')
assert(cline[0] == 'c')
# need to remove the last character (a comma) from the number for the float call
assert(math.fabs(float(cline[1][:-1]) - 60.) / 60. < 5.e-2) #5% relative error

# check d -- not empirical
dline = tokens[4].split('=')
assert(dline[0] == 'd')
# need to remove the last character (a comma) from the number for the float call
assert(math.fabs(float(dline[1][:-1]) - 10.) / 10. < 5.e-2) #5% relative error

# check norm -- this is empirical
normline = tokens[5].split('=')
assert(normline[0] == 'norm')
assert(math.fabs(float(normline[1]) - 567.) / 567. < 5.e-2)

## @endcond
