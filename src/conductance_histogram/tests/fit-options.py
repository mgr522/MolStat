# This file is a part of MolStat, which is distributed under the Creative
# Commons Attribution-NonCommercial 4.0 International Public License.
# MolStat (c) 2014, Northwestern University.

##
 # @file tests/fit-options.py
 # @brief Test suite for input options (initial guess specification) in the
 #    conductance histogram fitter.
 # 
 # @test Test suite for input options (initial guess specification) in the
 #    conductance histogram fitter.
 #
 # The symmetric, resonant-tunneling model; the symmetric, non-resonant
 # tunneling model; and the asymmetric, resonant-tunneling model are all
 # tested with various initial guess combinations.
 #
 # @author Matthew G.\ Reuter
 # @date June 2014

import subprocess

## @cond

# test for the symmetric, non-resonant model
process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricNonresonant\n' \
'symmetric-nonresonant.dat\n' \
'print\n' \
'guess c 60. d 10.' \
)

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Iter=  0, c=6.0000e+01, d=1.0000e+01, norm=1.0000e+00\n' \
'Iter=  1, c=3.9143e+01, d=6.4286e+00, norm=1.1799e+01\n' \
'Iter=  2, c=5.9172e+01, d=9.8836e+00, norm=1.1650e+01\n' \
'Iter=  3, c=5.8225e+01, d=9.6957e+00, norm=1.1796e+01\n' \
'Iter=  4, c=5.8240e+01, d=9.6983e+00, norm=1.1800e+01\n' \
'Iter=  5, c=5.8240e+01, d=9.6984e+00, norm=1.1800e+01\n' \
'Residual = 2.752039e+01\n' \
'\n' \
'Resid = 2.752039e+01\n' \
'c=5.8240e+01, d=9.6984e+00, norm=1.1800e+01\n' \
)

# test for the symmetric, resonant tunneling model
process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricResonant\n' \
'symmetric-resonant.dat\n' \
'print\n' \
'guess gamma 8. norm 3.' \
)

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Iter=  0, gamma=8.0000e+00, norm=3.0000e+00\n' \
'Iter=  1, gamma=1.0225e+01, norm=3.8700e+00\n' \
'Iter=  2, gamma=9.9417e+00, norm=3.9589e+00\n' \
'Iter=  3, gamma=9.9460e+00, norm=3.9567e+00\n' \
'Iter=  4, gamma=9.9462e+00, norm=3.9567e+00\n' \
'Residual = 4.310146e-01\n' \
'\n' \
'Resid = 4.310146e-01\n' \
'gamma=9.9462e+00, norm=3.9567e+00\n' \
)

# test for the asymmetric, resonant tunneling model
process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'AsymmetricResonant\n' \
'asymmetric-resonant.dat\n' \
'print\n' \
'guess gammal 17. gammar 12. r 25. norm 8.' \
)

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Iter=  0, gammaL=1.7000e+01, gammaR=1.2000e+01, r=2.5000e+01, norm=8.0000e+00\n' \
'Iter=  1, gammaL=1.6531e+01, gammaR=1.1350e+01, r=2.5219e+01, norm=8.1379e+00\n' \
'Iter=  2, gammaL=1.7103e+01, gammaR=1.1748e+01, r=2.3902e+01, norm=7.7229e+00\n' \
'Iter=  3, gammaL=1.7102e+01, gammaR=1.1749e+01, r=2.3910e+01, norm=7.7203e+00\n' \
'Iter=  4, gammaL=1.7103e+01, gammaR=1.1751e+01, r=2.3909e+01, norm=7.7142e+00\n' \
'Iter=  5, gammaL=1.7103e+01, gammaR=1.1751e+01, r=2.3909e+01, norm=7.7141e+00\n' \
'Residual = 5.084054e+00\n' \
'\n' \
'Resid = 5.084054e+00\n' \
'gammaL=1.7103e+01, gammaR=1.1751e+01, r=2.3909e+01, norm=7.7141e+00\n' \
)

## @endcond
