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
 # @date September 2014

import subprocess

## @cond

# test for the symmetric, non-resonant model
process = subprocess.Popen('../../molstat-fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricNonresonant\n' \
'symmetric-nonresonant.dat\n' \
'print\n' \
'guess c 60. d 10. norm 1.e5' \
)

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == \
'Iter=  0, c=6.0000e+01, d=1.0000e+01, norm=1.0000e+05\n' \
'Iter=  1, c=5.9241e+01, d=9.8731e+00, norm=1.0205e+05\n' \
'Iter=  2, c=5.9267e+01, d=9.8774e+00, norm=1.0206e+05\n' \
'Iter=  3, c=5.9267e+01, d=9.8774e+00, norm=1.0206e+05\n' \
'Residual = 1.599687e+04\n' \
'\n' \
'Resid = 1.599687e+04\n' \
'c=5.9267e+01, d=9.8774e+00, norm=1.0206e+05\n' \
)

# test for the symmetric, resonant tunneling model
process = subprocess.Popen('../../molstat-fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
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
'Iter=  1, gamma=9.9347e+00, norm=8.6442e+02\n' \
'Iter=  2, gamma=9.8875e+00, norm=8.7932e+02\n' \
'Iter=  3, gamma=9.8877e+00, norm=8.7813e+02\n' \
'Iter=  4, gamma=9.8877e+00, norm=8.7814e+02\n' \
'Residual = 3.922257e+00\n' \
'\n' \
'Resid = 3.922257e+00\n' \
'gamma=9.8877e+00, norm=8.7814e+02\n' \
)

# test for the asymmetric, resonant tunneling model
process = subprocess.Popen('../../molstat-fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
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
'Iter=  1, gammaL=1.7946e+01, gammaR=1.2943e+01, r=6.3351e-01, norm=3.8541e+01\n' \
'Iter=  2, gammaL=1.7851e+01, gammaR=1.2814e+01, r=7.0020e-01, norm=4.1773e+01\n' \
'Iter=  3, gammaL=1.7800e+01, gammaR=1.2762e+01, r=7.2074e-01, norm=4.2944e+01\n' \
'Iter=  4, gammaL=1.7799e+01, gammaR=1.2761e+01, r=7.2187e-01, norm=4.3013e+01\n' \
'Iter=  5, gammaL=1.7799e+01, gammaR=1.2761e+01, r=7.2185e-01, norm=4.3011e+01\n' \
'Residual = 2.788663e+02\n' \
'\n' \
'Resid = 2.788663e+02\n' \
'gammaL=1.7799e+01, gammaR=1.2761e+01, r=7.2185e-01, norm=4.3011e+01\n' \
)

## @endcond
