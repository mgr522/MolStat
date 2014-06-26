# This file is a part of MolStat, which is distributed under the Creative
# Commons Attribution-NonCommercial 4.0 International Public License.
# MolStat (c) 2014, Northwestern University.

##
 # @file tests/fit-symmetric-nonresonant.py
 # @brief Test suite for fitting to the symmetric, nonresonant-tunneling model.
 # 
 # @test Test suite for fitting to the symmetric, nonresonant-tunneling model.
 #
 # The data in symmetric-nonresonant.dat was generated by the simulator with
 # the following input
 # @verbatim
 # SymmetricOneSiteModel
 # ZeroBias
 # 1000000
 # 100 log 10.
 # 0.
 # gamma normal 0.5 0.05
 # epsilon normal 3. 0.05
 # @endverbatim
 #
 # @author Matthew G.\ Reuter
 # @date May 2014

import subprocess

## @cond

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricNonresonant\n' \
'symmetric-nonresonant.dat\n')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Resid = 4.519332e+03\n' \
'c=5.9138e+01, d=9.8513e+00, norm=1.0147e+04\n')

## @endcond
