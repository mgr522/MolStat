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
 # See the documentation in fit-symmetric-nonresonant.py,
 # fit-symmetric-resonant.py, and fit-asymmetric-nonresonant.py for
 # information on how the input data was generated.
 #
 # @author Matthew G.\ Reuter
 # @date June 2014

import subprocess

## @cond

# test for the symmetric, non-resonant model
process = subprocess.Popen('../../molstat-fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricNonresonant\n' \
'symmetric-nonresonant.dat\n' \
'print\n' \
'guess c 60. d 10. norm 1.e4' \
)

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Iter=  0, c=6.0000e+01, d=1.0000e+01, norm=1.0000e+04\n' \
'Iter=  1, c=5.9112e+01, d=9.8467e+00, norm=1.0146e+04\n' \
'Iter=  2, c=5.9138e+01, d=9.8512e+00, norm=1.0147e+04\n' \
'Iter=  3, c=5.9138e+01, d=9.8513e+00, norm=1.0147e+04\n' \
'Residual = 4.519332e+03\n' \
'\n' \
'Resid = 4.519332e+03\n' \
'c=5.9138e+01, d=9.8513e+00, norm=1.0147e+04\n' \
)

# test for the symmetric, resonant tunneling model
process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricResonant\n' \
'symmetric-resonant.dat\n' \
'print\n' \
'guess gamma 9.7 norm 570.' \
)

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Iter=  0, gamma=9.7000e+00, norm=5.7000e+02\n' \
'Iter=  1, gamma=9.7117e+00, norm=5.7172e+02\n' \
'Iter=  2, gamma=9.7136e+00, norm=5.7256e+02\n' \
'Iter=  3, gamma=9.7139e+00, norm=5.7269e+02\n' \
'Iter=  4, gamma=9.7140e+00, norm=5.7271e+02\n' \
'Residual = 3.516796e+00\n' \
'\n' \
'Resid = 3.516796e+00\n' \
'gamma=9.7140e+00, norm=5.7271e+02\n' \
)

# test for the asymmetric, resonant tunneling model
process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'AsymmetricResonant\n' \
'asymmetric-resonant.dat\n' \
'print\n' \
'guess gammal 17. gammar 12. r 0.8 norm 44.' \
)

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Iter=  0, gammaL=1.7000e+01, gammaR=1.2000e+01, r=8.0000e-01, norm=4.4000e+01\n' \
'Iter=  1, gammaL=1.7803e+01, gammaR=1.2826e+01, r=6.6801e-01, norm=4.0152e+01\n' \
'Iter=  2, gammaL=1.7716e+01, gammaR=1.2705e+01, r=7.2670e-01, norm=4.3077e+01\n' \
'Iter=  3, gammaL=1.7688e+01, gammaR=1.2675e+01, r=7.4131e-01, norm=4.3908e+01\n' \
'Iter=  4, gammaL=1.7688e+01, gammaR=1.2675e+01, r=7.4177e-01, norm=4.3937e+01\n' \
'Iter=  5, gammaL=1.7688e+01, gammaR=1.2675e+01, r=7.4176e-01, norm=4.3936e+01\n' \
'Residual = 2.759250e+02\n' \
'\n' \
'Resid = 2.759250e+02\n' \
'gammaL=1.7688e+01, gammaR=1.2675e+01, r=7.4176e-01, norm=4.3936e+01\n' \
)

## @endcond
