import subprocess

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate('SymmetricNonresonant\nsymmetric-nonresonant.dat\nnoprint')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 2.752039e+01\nc=5.8240e+01, d=9.6984e+00, norm=1.1800e+01\n')
