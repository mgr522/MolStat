import subprocess

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
