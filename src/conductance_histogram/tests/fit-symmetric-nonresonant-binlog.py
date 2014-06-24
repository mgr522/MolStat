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
 # @author Matthew G.\ Reuter
 # @date May 2014

import subprocess

# test the log-binning aspect of the fitter program

## @cond

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricNonresonant\n' \
'symmetric-nonresonant-binlog.dat\n' \
'bin log 10.')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 2.752118e+01\nc=5.8240e+01, d=9.6983e+00, norm=1.1800e+01\n')

## @endcond
