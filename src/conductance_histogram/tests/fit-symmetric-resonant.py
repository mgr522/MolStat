# This file is a part of MolStat, which is distributed under the Creative
# Commons Attribution-NonCommercial 4.0 International Public License.
# MolStat (c) 2014, Northwestern University.

##
 # @file tests/fit-symmetric-resonant.py
 # @brief Test suite for fitting to the symmetric, resonant-tunneling model.
 # 
 # @test Test suite for fitting to the symmetric, resonant-tunneling model.
 #
 # @author Matthew G.\ Reuter
 # @date May 2014

import subprocess

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate('SymmetricResonant\nsymmetric-resonant.dat\nnoprint')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 4.310146e-01\ngamma=9.9462e+00, norm=3.9567e+00\n')
